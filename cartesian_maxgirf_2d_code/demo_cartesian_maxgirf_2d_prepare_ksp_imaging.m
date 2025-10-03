% demo_cartesian_maxgirf_2d_prepare_ksp_imaging.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 01/13/2025, Last modified: 02/08/2025

%% Clean slate
close all; clearvars -except json_number nr_json_files json_files json_file grad_file_path; clc;

%% Start a stopwatch timer
start_time = tic;

%% Read a .json file
tstart = tic; fprintf('%s: Reading a .json file: %s... ', datetime, json_file);
fid = fopen(json_file); 
json_txt = fread(fid, [1 inf], 'char=>char'); 
fclose(fid); 
json = jsondecode(json_txt);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Define the full path of a filename
%--------------------------------------------------------------------------
if ispc
    siemens_twix_file  = strrep(json.siemens_twix_file, '/', '\');
    ismrmrd_data_file  = strrep(json.ismrmrd_data_file, '/', '\');
    ismrmrd_noise_file = strrep(json.ismrmrd_noise_file, '/', '\');
    output_path        = strrep(json.output_path, '/', '\');
else
    siemens_twix_file  = json.siemens_twix_file;
    ismrmrd_data_file  = json.ismrmrd_data_file;
    ismrmrd_noise_file = json.ismrmrd_noise_file;
    output_path        = json.output_path;
end

%--------------------------------------------------------------------------
% Define the BART path
%--------------------------------------------------------------------------
bart_path = json.bart_path;

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
Lmax              = json.recon_parameters.Lmax;              % maximum rank of the SVD approximation of a higher-order encoding matrix
cutoff_percentage = json.recon_parameters.cutoff_percentage; % cutoff percentage to calculate an approximate full rank of a higher-order encoding matrix
slice_type        = json.recon_parameters.slice_type;        % type of an excitation slice: "curved" vs "flat"
phc_flag          = json.recon_parameters.phc_flag;          % 1=yes, 0=no
gridding_flag     = json.recon_parameters.gridding_flag;     % 1=yes, 0=no
cfc_flag          = json.recon_parameters.cfc_flag;          % 1=yes, 0=no
sfc_flag          = json.recon_parameters.sfc_flag;          % 1=yes, 0=no
gnc_flag          = json.recon_parameters.gnc_flag;          % 1=yes, 0=no
topup_flag        = json.recon_parameters.topup_flag;        % 1=yes, 0=no

%% Make an output path
mkdir(output_path);

%% Set up BART commands
%--------------------------------------------------------------------------
% Define a BART command
%--------------------------------------------------------------------------
if ispc
    command_prefix = 'wsl';
else
    command_prefix = '';
end
bart_command = sprintf('%s %s/bart', command_prefix, bart_path);

%--------------------------------------------------------------------------
% Translate from a Windows path to a WSL path 
%--------------------------------------------------------------------------
if ispc
    bart_output_path = strrep(output_path, '\', '/');
    bart_output_path = sprintf('/mnt/%s/%s/', lower(bart_output_path(1)), bart_output_path(4:end));
else
    bart_output_path = sprintf('%s/', output_path);
end

%% Read an ISMRMRD file (k-space data)
tstart = tic; fprintf('%s: Reading an ISMRMRD file: %s... ', datetime, ismrmrd_data_file);
if exist(ismrmrd_data_file, 'file')
    dset = ismrmrd.Dataset(ismrmrd_data_file, 'dataset');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
else
    error('File %s does not exist.  Please generate it.' , ismrmrd_data_file);
end

%% Get imaging parameters from an XML header
header = ismrmrd.xml.deserialize(dset.readxml);

%--------------------------------------------------------------------------
% measurement information
%--------------------------------------------------------------------------
patient_position = header.measurementInformation.patientPosition;

%--------------------------------------------------------------------------
% experimental conditions
%--------------------------------------------------------------------------
gamma = 4257.59 * (1e4 * 2 * pi); % gyromagnetic ratio for 1H [rad/sec/T]
B0 = header.experimentalConditions.H1resonanceFrequency_Hz * (2 * pi / gamma); % [Hz] * [2pi rad/cycle] / [rad/sec/T] => [T]

%--------------------------------------------------------------------------
% sequence parameters
%--------------------------------------------------------------------------
TE = header.sequenceParameters.TE * 1e-3;

%--------------------------------------------------------------------------
% Encoding Space (Nkx, Nky, Nkz)
%--------------------------------------------------------------------------
encoded_fov(1) = header.encoding.encodedSpace.fieldOfView_mm.x * 1e-3; % [m] RO
encoded_fov(2) = header.encoding.encodedSpace.fieldOfView_mm.y * 1e-3; % [m] PE
encoded_fov(3) = header.encoding.encodedSpace.fieldOfView_mm.z * 1e-3; % [m] SL

Nkx = header.encoding.encodedSpace.matrixSize.x; % number of readout samples in k-space
Nky = header.encoding.encodedSpace.matrixSize.y; % number of phase-encoding steps in k-space
Nkz = header.encoding.encodedSpace.matrixSize.z; % number of slice-encoding steps in k-space

encoded_resolution = encoded_fov ./ [Nkx Nky Nkz]; % [m]

kspace_encoding_step_1_center = header.encoding.encodingLimits.kspace_encoding_step_1.center;
kspace_encoding_step_2_center = header.encoding.encodingLimits.kspace_encoding_step_2.center;

kspace_encoding_step_1_maximum = header.encoding.encodingLimits.kspace_encoding_step_1.maximum;
kspace_encoding_step_1_minimum = header.encoding.encodingLimits.kspace_encoding_step_1.minimum;

kspace_encoding_step_2_maximum = header.encoding.encodingLimits.kspace_encoding_step_2.maximum;
kspace_encoding_step_2_minimum = header.encoding.encodingLimits.kspace_encoding_step_2.minimum;

acceleration_factor1 = header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1;
acceleration_factor2 = header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_2;

%--------------------------------------------------------------------------
% Recon Space (Nx, Ny, Nz)
%--------------------------------------------------------------------------
recon_fov(1) = header.encoding.reconSpace.fieldOfView_mm.x * 1e-3; % [m] RO
recon_fov(2) = header.encoding.reconSpace.fieldOfView_mm.y * 1e-3; % [m] PE
recon_fov(3) = header.encoding.reconSpace.fieldOfView_mm.z * 1e-3; % [m] SL

Nx = header.encoding.reconSpace.matrixSize.x; % number of samples in image space (RO)
Ny = header.encoding.reconSpace.matrixSize.y; % number of samples in image space (PE)
Nz = header.encoding.reconSpace.matrixSize.z; % number of samples in image space (SL)

recon_resolution = recon_fov ./ [Nx Ny Nz]; % [m]

%--------------------------------------------------------------------------
% Readout oversampling factor
%--------------------------------------------------------------------------
readout_os_factor = encoded_fov(1) / recon_fov(1); % readout oversampling factor

%--------------------------------------------------------------------------
% Number of receive channels
%--------------------------------------------------------------------------
Nc = header.acquisitionSystemInformation.receiverChannels;

%% Read acquisitions
tstart = tic; fprintf('%s: Reading acquisitions... ', datetime);
raw_data = dset.readAcquisition(); % read all the acquisitions
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Get data type
%--------------------------------------------------------------------------
acq_is_noise_measurement                = raw_data.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
acq_is_parallel_calibration             = raw_data.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION');
acq_is_parallel_calibration_and_imaging = raw_data.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING');
acq_is_reverse                          = raw_data.head.flagIsSet('ACQ_IS_REVERSE');
acq_is_navigation_data                  = raw_data.head.flagIsSet('ACQ_IS_NAVIGATION_DATA');
acq_is_phasecorr_data                   = raw_data.head.flagIsSet('ACQ_IS_PHASECORR_DATA');
acq_is_hpfeedback_data                  = raw_data.head.flagIsSet('ACQ_IS_HPFEEDBACK_DATA');
acq_is_dummyscan_data                   = raw_data.head.flagIsSet('ACQ_IS_DUMMYSCAN_DATA');
acq_is_rtfeedback_data                  = raw_data.head.flagIsSet('ACQ_IS_RTFEEDBACK_DATA');
acq_is_surfacecoilcorrectionscan_data   = raw_data.head.flagIsSet('ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA');

%% Get navigator data for imaging and imaging data
% navigator data = nav for calibration + nav for imaging
nav_data = raw_data.select(find(acq_is_phasecorr_data));

%% Get imaging data
img_data = raw_data.select(find(~acq_is_noise_measurement & ~acq_is_parallel_calibration & ~acq_is_navigation_data & ~acq_is_phasecorr_data & ~acq_is_dummyscan_data));

%% Parse an ISMRMRD header
%--------------------------------------------------------------------------
% ISMRMRD header
%--------------------------------------------------------------------------
% uint16_t version;                                    /**< First unsigned int indicates the version */
% uint64_t flags;                                      /**< bit field with flags */
% uint32_t measurement_uid;                            /**< Unique ID for the measurement */
% uint32_t scan_counter;                               /**< Current acquisition number in the measurement */
% uint32_t acquisition_time_stamp;                     /**< Acquisition clock */
% uint32_t physiology_time_stamp[ISMRMRD_PHYS_STAMPS]; /**< Physiology time stamps, e.g. ecg, breating, etc. */
% uint16_t number_of_samples;                          /**< Number of samples acquired */
% uint16_t available_channels;                         /**< Available coils */
% uint16_t active_channels;                            /**< Active coils on current acquisiton */
% uint64_t channel_mask[ISMRMRD_CHANNEL_MASKS];        /**< Mask to indicate which channels are active. Support for 1024 channels */
% uint16_t discard_pre;                                /**< Samples to be discarded at the beginning of  acquisition */
% uint16_t discard_post;                               /**< Samples to be discarded at the end of acquisition */
% uint16_t center_sample;                              /**< Sample at the center of k-space */
% uint16_t encoding_space_ref;                         /**< Reference to an encoding space, typically only one per acquisition */
% uint16_t trajectory_dimensions;                      /**< Indicates the dimensionality of the trajectory vector (0 means no trajectory) */
% float sample_time_us;                                /**< Time between samples in micro seconds, sampling BW */
% float position[3];                                   /**< Three-dimensional spatial offsets from isocenter */
% float read_dir[3];                                   /**< Directional cosines of the readout/frequency encoding */
% float phase_dir[3];                                  /**< Directional cosines of the phase */
% float slice_dir[3];                                  /**< Directional cosines of the slice direction */
% float patient_table_position[3];                     /**< Patient table off-center */
% ISMRMRD_EncodingCounters idx;                        /**< Encoding loop counters, see above */
% int32_t user_int[ISMRMRD_USER_INTS];                 /**< Free user parameters */
% float user_float[ISMRMRD_USER_FLOATS];               /**< Free user parameters */
%--------------------------------------------------------------------------
% Where EncodingCounters are defined as:
% uint16_t kspace_encode_step_1;    /**< e.g. phase encoding line number */
% uint16_t kspace_encode_step_2;    /**< e.g. partition encoding number */
% uint16_t average;                 /**< e.g. signal average number */
% uint16_t slice;                   /**< e.g. imaging slice number */
% uint16_t contrast;                /**< e.g. echo number in multi-echo */
% uint16_t phase;                   /**< e.g. cardiac phase number */
% uint16_t repetition;              /**< e.g. dynamic number for dynamic scanning */
% uint16_t set;                     /**< e.g. flow encoding set */
% uint16_t segment;                 /**< e.g. segment number for segmented acquisition */
% uint16_t user[ISMRMRD_USER_INTS]; /**< Free user parameters */
%--------------------------------------------------------------------------
number_of_samples       = double(max(img_data.head.number_of_samples));
discard_pre             = double(max(img_data.head.discard_pre));
discard_post            = double(max(img_data.head.discard_post));
center_sample           = double(max(img_data.head.center_sample));
nr_channels             = double(max(img_data.head.active_channels));
nr_phase_encoding_steps = double(max(img_data.head.idx.kspace_encode_step_1)) + 1;
nr_slice_encoding_steps = double(max(img_data.head.idx.kspace_encode_step_2)) + 1;
nr_averages             = double(max(img_data.head.idx.average)) + 1;
nr_slices               = double(max(img_data.head.idx.slice)) + 1;
nr_contrasts            = double(max(img_data.head.idx.contrast)) + 1;
nr_phases               = double(max(img_data.head.idx.phase)) + 1;
nr_repetitions          = double(max(img_data.head.idx.repetition)) + 1;
nr_sets                 = double(max(img_data.head.idx.set)) + 1;
nr_segments             = double(max(img_data.head.idx.segment)) + 1;

%--------------------------------------------------------------------------
% Get the dimensionality of the trajectory vector (0 means no trajectory)
%--------------------------------------------------------------------------
trajectory_dimensions = double(max(img_data.head.trajectory_dimensions));

%% Display an ISMRMRD header
fprintf('========================= ISMRMRD header =========================\n');
fprintf('encoded_fov        = %8.4f %8.4f %8.4f [mm]\n', encoded_fov(1) * 1e3, encoded_fov(2) * 1e3, encoded_fov(3) * 1e3);
fprintf('Nkx Nky Nkz        = %d      %d        %d\n', Nkx, Nky, Nkz);
fprintf('encoded_resolution = %8.4f %8.4f %8.4f [mm]\n', encoded_resolution(1) * 1e3, encoded_resolution(2) * 1e3, encoded_resolution(3) * 1e3);
fprintf('------------------------------------------------------------------\n');
fprintf('recon_fov          = %8.4f %8.4f %8.4f [mm]\n', recon_fov(1) * 1e3, recon_fov(2) * 1e3, recon_fov(3) * 1e3);
fprintf('Nx Ny Nz           = %d      %d        %d\n', Nx, Ny, Nz);
fprintf('recon_resolution   = %8.4f %8.4f %8.4f [mm]\n', recon_resolution(1) * 1e3, recon_resolution(2) * 1e3, recon_resolution(3) * 1e3);
fprintf('------------------------------------------------------------------\n');
fprintf('trajectory              = %s\n', header.encoding.trajectory);
fprintf('number_of_samples       = %d\n', number_of_samples);
fprintf('discard_pre             = %d\n', discard_pre);
fprintf('discard_post            = %d\n', discard_post);
fprintf('center_sample           = %d\n', center_sample);
fprintf('nr_channels             = %d\n', nr_channels);
fprintf('nr_phase_encoding_steps = %d\n', nr_phase_encoding_steps);
fprintf('nr_slice_encoding_steps = %d\n', nr_slice_encoding_steps);
fprintf('nr_averages             = %d\n', nr_averages);
fprintf('nr_slices               = %d\n', nr_slices);
fprintf('nr_contrasts            = %d\n', nr_contrasts);
fprintf('nr_phases               = %d\n', nr_phases);
fprintf('nr_repetitions          = %d\n', nr_repetitions);
fprintf('nr_sets                 = %d\n', nr_sets); % flow encoding set
fprintf('nr_segments             = %d\n', nr_segments); % even vs odd lines (ACQ_IS_REVERSE)
fprintf('==================================================================\n');

%% Get a list of segments
%--------------------------------------------------------------------------
% Get a list of navigator acquisitions
%--------------------------------------------------------------------------
nav_acq_list = find((nav_data.head.idx.slice == 0) & (nav_data.head.idx.repetition == 0));

%--------------------------------------------------------------------------
% Get a list of segments (navigator)
%--------------------------------------------------------------------------
nav_segment_list = double(nav_data.head.idx.segment(nav_acq_list));

%--------------------------------------------------------------------------
% Get a list of imaging acquisitions
%--------------------------------------------------------------------------
img_acq_list = find((img_data.head.idx.slice == 0) & (img_data.head.idx.repetition == 0));

%--------------------------------------------------------------------------
% Get a list of segments (imaging)
%--------------------------------------------------------------------------
img_segment_list = double(img_data.head.idx.segment(img_acq_list));

%% Define parameters for image reconstruction
%--------------------------------------------------------------------------
% Set the total number of voxels in image space
%--------------------------------------------------------------------------
N = Nkx * Nky * Nkz;

%--------------------------------------------------------------------------
% Set the number of spatial basis functions
%--------------------------------------------------------------------------
Nl = 19;

%--------------------------------------------------------------------------
% Set the number of phase-encoding steps
%--------------------------------------------------------------------------
Ni = length(img_acq_list);

%--------------------------------------------------------------------------
% Set the total number of k-space samples
%--------------------------------------------------------------------------
Nk = Nkx * Ni;

%--------------------------------------------------------------------------
% Set the oversampling parameter for randomized SVD
%--------------------------------------------------------------------------
os = 5;

%% Calculate the receiver noise matrix (Nc x Nc)
[Psi,inv_L] = calculate_receiver_noise_matrix(ismrmrd_noise_file);

%% Read a Siemens .dat file
fprintf('%s: Reading a Siemens .dat file: %s\n', datetime, siemens_twix_file);
twix = mapVBVD(siemens_twix_file);

%% Get the name of a gradient set
if isfield(twix{1}.hdr.Meas, 'tGradientCoil')
    tGradientCoil = twix{1}.hdr.Meas.tGradientCoil;
elseif isfield(twix{2}.hdr.Meas, 'tGradientCoil')
    tGradientCoil = twix{2}.hdr.Meas.tGradientCoil;
end

%% Reduce the TWIX dataset
if length(twix) > 1
    twix = twix{end};
end

%% Get a slice normal vector from Siemens TWIX format
%--------------------------------------------------------------------------
% dNormalSag: Sagittal component of a slice normal vector (in the PCS)
%--------------------------------------------------------------------------
if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal, 'dSag')
    dNormalSag = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal.dSag;
else
    dNormalSag = 0;
end

%--------------------------------------------------------------------------
% dNormalCor: Coronal component of a slice normal vector (in the PCS)
%--------------------------------------------------------------------------
if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal, 'dCor')
    dNormalCor = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal.dCor;
else
    dNormalCor = 0;
end

%--------------------------------------------------------------------------
% dNormalTra: Transverse component of a slice normal vector (in the PCS)
%--------------------------------------------------------------------------
if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal, 'dTra')
    dNormalTra = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal.dTra;
else
    dNormalTra = 0;
end

%--------------------------------------------------------------------------
% dRotAngle: Slice rotation angle ("swap Fre/Pha")
%--------------------------------------------------------------------------
if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{1}, 'dInPlaneRot')
    dRotAngle = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.dInPlaneRot; % [rad]
else
    dRotAngle = 0; % [rad]
end

%% Calculate a rotation matrix from the GCS to the PCS
[R_gcs2pcs, phase_sign, read_sign, main_orientation] = siemens_calculate_transform_gcs_to_pcs(dNormalSag, dNormalCor, dNormalTra, dRotAngle);

%% Get a rotation matrix from the GCS to the PCS (ISMRMRD format)
phase_dir = double(img_data.head.phase_dir(:,1));
read_dir  = double(img_data.head.read_dir(:,1));
slice_dir = double(img_data.head.slice_dir(:,1));
R_gcs2pcs_ismrmrd = [phase_dir read_dir slice_dir];

%% Calculate a transformation matrix from the RCS to the GCS [r,c,s] <=> [PE,RO,SL]
R_rcs2gcs = [0    1    0 ; % [PE]   [0 1 0] * [r]
             1    0    0 ; % [RO] = [1 0 0] * [c]
             0    0    1]; % [SL]   [0 0 1] * [s]

%% Calculate a rotation matrix from the PCS to the DCS
R_pcs2dcs = siemens_calculate_transform_pcs_to_dcs(patient_position);

%% Calculate a rotation matrix from the GCS to the DCS
R_gcs2dcs = R_pcs2dcs * R_gcs2pcs_ismrmrd;

%% Calculate a rotation matrix from [RO,PE,SL] to [x,y,z]
%--------------------------------------------------------------------------
% [x]   [         ][0 1 0][RO]
% [y] = [R_gcs2dcs][1 0 0][PE]
% [z]   [         ][0 0 1][SL]
%--------------------------------------------------------------------------
R_rps2dcs = R_gcs2dcs * [0 1 0; 1 0 0; 0 0 1];

%% Get information about a conventional EPI trajectory
%--------------------------------------------------------------------------
%
%     rampUpTime   flatTopTime  rampDownTime
%       |      |_______________|      |
%       |     /|               |\     |
%       |    / |               | \    |
%       |   /  |               |  \   |
%       |  /   |               |   \  |
%       | /    |               |    \ |
%  _____|/     |               |     \|                             ______ 
%       |<---->|<------------->|<---->|\                           /
%       |  90        1050         90  | \                         /
%       |<--------------------------->|  \                       /
%       |     echoSpacing = 1230      |   \                     /
%       |<->|xxxxxxxxxxxxxxxxxxxxx|<->|    \                   /
%           numSamples * dwellTime          \_________________/
%        34 uses  = 1161.6 usec    34.4 usec
%
% (echoSpacing - (numSamples * dwellTime)) / 2 = 34.2 usec
% acqDelayTime = 34 usec
%--------------------------------------------------------------------------
etl            = header.encoding.trajectoryDescription.userParameterLong(1).value;        % etl
nr_navigators  = header.encoding.trajectoryDescription.userParameterLong(2).value;        % numberOfNavigators
ramp_up_time   = header.encoding.trajectoryDescription.userParameterLong(3).value * 1e-6; % rampUpTime [usec]
ramp_down_time = header.encoding.trajectoryDescription.userParameterLong(4).value * 1e-6; % rampDownTime [usec]
flat_top_time  = header.encoding.trajectoryDescription.userParameterLong(5).value * 1e-6; % flatTopTime [usec]
echo_spacing   = header.encoding.trajectoryDescription.userParameterLong(6).value * 1e-6; % echoSpacing [usec]
acq_delay_time = header.encoding.trajectoryDescription.userParameterLong(7).value * 1e-6; % acqDelayTime [usec]
num_samples    = header.encoding.trajectoryDescription.userParameterLong(8).value;        % numSamples
dwell_time     = header.encoding.trajectoryDescription.userParameterDouble.value * 1e-6;  % dwellTime [usec]

%--------------------------------------------------------------------------
% Correct the EPI factor (i.e., echo train length)
%--------------------------------------------------------------------------
etl = Ni;

%--------------------------------------------------------------------------
% Define the duration of an RF pulse
%--------------------------------------------------------------------------
rf_duration = 3e-3;

%% Calculate the magnitude of a readout gradient lobe [mT/m]
%--------------------------------------------------------------------------
%                         |<--------------------->| truncated trapezoid 
%                      |  |   |       |       |   |  |
%                  Gx  |  |   |_______|_______|   |  |
%                      |  |  /|       |       |\  |  |
%                      |  | / |       |       | \ |  |
%                      |  |/  |       |       |  \|  |
%                      |  /   |       |       |   \  |
%                      | /|   |       |       |   |\ |
%  ____________________|/_|___|_______|_______|___|_\|____________________
%                     -t2 |  -t1      0       t1  |  t2
%                      |<>| ta                    |<>| tb
%
% SR = Gx / (-t1 + t2)
%
% The area of a truncated trapezoid (from -t2 + ta to t2 - tb) gives the maximum
% k-space extent.
%
%               [ SR * (t + t2)      t in [-t2, -t1)
% G_{trap}(t) = [ Gx                 t in [-t1,  t1)
%               [ SR * (t2 - t)      t in [ t1,  t2)
%
% Area of a trapezoid = gamma * (2 * t2 + 2 * t1) / 2 * Gx
%                     = gamma * Gx * (t1 + t2)
%
% Area of a triangle1 = gamma * ta / 2 * SR * (-t2 + ta + t2)
%                     = gamma * SR * ta^2 / 2
%
% Area of a triangle2 = gamma * tb / 2 * SR * (t2 - (t2 - tb))
%                     = gamma * SR * tb^2 / 2
%
% Area of a truncated trapezoid
% = gamma * Gx * (t1 + t2) - gamma * SR * ta^2 / 2 - gamma * SR * tb^2 / 2
% = gamma * Gx * [(t1 + t2) - (ta^2 + tb^2) / (2 * (-t1 + t2))]
%--------------------------------------------------------------------------
t1 = flat_top_time / 2; % [sec]
t2 = echo_spacing / 2; % [sec]

%--------------------------------------------------------------------------
% Set ta
%--------------------------------------------------------------------------
ta = acq_delay_time;

%--------------------------------------------------------------------------
% Set tb
%--------------------------------------------------------------------------
tb = 2 * t2 - (ta + num_samples * dwell_time);

%--------------------------------------------------------------------------
% Calculate the readout gradient amplitude [mT/m]
%--------------------------------------------------------------------------
Gu = (2 * pi / encoded_resolution(1)) / ((gamma * 1e-3) * ((t1 + t2) - (ta^2 + tb^2) / (2 * (-t1 + t2)))); % [mT/m]

%--------------------------------------------------------------------------
% Set the slew rate
%--------------------------------------------------------------------------
SR = Gu / (-t1 + t2); % [mT/m/sec]

%% Calculate the amplitude of a readout gradient lobe in the RPS [mT/m]
g_rps = cat(1, Gu, 0, 0); % [RO,PE,SL]

%% Calculate the amplitude of a readout gradient lobe in the DCS [mT/m]
g_dcs = R_rps2dcs * g_rps;

%% Calculate a vertex-based readout gradient lobe [mT/m]
t_vertex_readout_gradient_lobe = [0 ramp_up_time ramp_up_time + flat_top_time ramp_up_time + flat_top_time + ramp_down_time].'; % [sec]
g_vertex_readout_gradient_lobe = cat(1, 0, Gu, Gu, 0); % [mT/m]

%% Create a list of phase-encoding gradient areas
deltak_phase = (2 * pi) / encoded_fov(2); % [rad/m]
phase_areas = (-floor(Nky/2):ceil(Nky/2)-1).' * deltak_phase; % [rad/m]

%% Create a list of partition-encoding gradient areas
deltak_partition = (2 * pi) / encoded_fov(3); % [rad/m]
partition_areas = (-floor(Nkz/2):ceil(Nkz/2)-1).' * deltak_partition; % [rad/m]

%% Calculate a list of phase encoding line numbers (starting from 1)
phase_encoding_line_number_list = double(img_data.head.idx.kspace_encode_step_1(img_acq_list)) - (kspace_encoding_step_1_maximum - kspace_encoding_step_1_center + 1) + floor(Nky/2) + 1;

%% Calculate a list of partition encoding numbers (starting from 1)
partition_encoding_number_list = double(img_data.head.idx.kspace_encode_step_2(img_acq_list)) - (kspace_encoding_step_2_maximum - kspace_encoding_step_2_center + 1) + floor(Nkz/2) + 1;

%% Calculate the number of shots
nr_shots = length(phase_encoding_line_number_list) / etl;

%% Calculate the number of readout gradient lobes played to acquire the center of k-space
for shot_number = 1:nr_shots
    shot_range = (1:etl).' + (shot_number - 1) * etl;
    phase_encoding_line_number_list_per_shot = phase_encoding_line_number_list(shot_range);
    N_center = find(phase_encoding_line_number_list_per_shot == (floor(Nky/2) + 1));
    if ~isempty(N_center)
        break;
    end
end

%% Calculate the start time of a train of alternating polarity readout gradient lobes
t_s = TE - N_center * echo_spacing + echo_spacing / 2;

%% Calculate a train of alternating polarity readout gradient lobes [mT/m]
%--------------------------------------------------------------------------
%                                        TE
%                                        |
%               eco = 1   eco = 2    eco = 3     eco = 4
%               ______                ___|__
%              /|    |\              /|  | |\
%             / |    | \            / |  | | \
%            *--*----*--*--*----*--*--*--|-*--*--*----*--*-----> t
%           ts           \ |    | /      |     \ |    | /
%                         \|____|/       |      \|____|/
%            1  2    3  1  2    3  1  2    3  1  2    3  4
%--------------------------------------------------------------------------
nr_vertices_train = 3 * etl + 1;

t_vertex_train = zeros(nr_vertices_train, 1, 'double');
g_vertex_train = zeros(nr_vertices_train, 1, 'double');

%--------------------------------------------------------------------------
% Calculate the gradient lobes of a EPI readout waveform
%--------------------------------------------------------------------------
for eco = 1:etl

    %----------------------------------------------------------------------
    % Set the sign of a gradient lobe
    %----------------------------------------------------------------------
    if img_segment_list(eco) == 0 % forward (segment = 0)
        echo_sign = +1;
    else % reverse (segment == 1)
        echo_sign = -1;
    end

    %----------------------------------------------------------------------
    %               eco = 1   eco = 2    eco = 3     eco = 4
    %               ______                ______
    %              /|    |\              /|    |\
    %             / |    | \            / |    | \
    %            *--*----*--*--*----*--*--*----*--*--*----*--*-----
    %                        \ |    | /            \ |    | /
    %                         \|____|/              \|____|/
    %            1  2    3  1  2    3  1  2    3  1  2    3  4 
    %----------------------------------------------------------------------
    if eco < etl
        index_range = (1:3).';
    else
        index_range = (1:4).';
    end
    index_offset = (eco - 1) * 3;

    %----------------------------------------------------------------------
    % Calculate the time samples of a gradient lobe
    %----------------------------------------------------------------------
    t_vertex_train(index_range + index_offset) = t_s + t_vertex_readout_gradient_lobe(index_range) + (eco - 1) * echo_spacing;

    %----------------------------------------------------------------------
    % Calculate the gradient samples of a gradient lobe
    %----------------------------------------------------------------------
    g_vertex_train(index_range + index_offset) = echo_sign * g_vertex_readout_gradient_lobe(index_range);
end

%% Calculate a vertex-based prephasing gradient lobe [mT/m]
if img_segment_list(1) == 0 % forward (segment = 0)
    echo_sign = +1;
else % reverse (segment == 1)
    echo_sign = -1;
end

t_vertex_prephasing_gradient_lobe = cat(1, -echo_spacing + t_s, -echo_spacing + t_s + ramp_up_time, t_s - ramp_up_time, t_s);
g_vertex_prephasing_gradient_lobe = cat(1, 0, -echo_sign * Gu / 2, -echo_sign * Gu / 2, 0); % [mT/m]

%% Calculate a list of phase-encoding dephasers [mT/m]
%--------------------------------------------------------------------------
% Calculate the duration of a phase-encoding dephaser
%--------------------------------------------------------------------------
phase_encoding_dephaser_duration = ramp_up_time + flat_top_time + ramp_down_time; % [sec]

%--------------------------------------------------------------------------
% Calculate the area of a phase-encoding dephaser [mT/m*sec]
%--------------------------------------------------------------------------
phase_encoding_dephaser_area = zeros(nr_shots, 1, 'double');
for shot_number = 1:nr_shots
    shot_range = (1:etl).' + (shot_number - 1) * etl;
    phase_encoding_line_number_list_per_shot = phase_encoding_line_number_list(shot_range);
    % [rad/m] / ([rad/sec/T] * [T/1e3mT]) => [mT/m*sec]
    phase_encoding_dephaser_area(shot_number) = phase_areas(phase_encoding_line_number_list_per_shot(1)) / (gamma * 1e-3);
end

%--------------------------------------------------------------------------
% Calculate the amplitude of a phase-encoding dephaser [mT/m]
% area = amplitude * (dephaser_duration + flat_top_time) / 2
%--------------------------------------------------------------------------
G_phase_encoding_dephaser = 2 * phase_encoding_dephaser_area / (phase_encoding_dephaser_duration + flat_top_time);

%--------------------------------------------------------------------------
% Calculate the time samples of a phase-encoding dephaser [sec]
%--------------------------------------------------------------------------
t_vertex_phase_encoding_dephaser = cat(1, -echo_spacing + t_s, -echo_spacing + t_s + ramp_up_time, t_s - ramp_down_time, t_s);

%--------------------------------------------------------------------------
% Calculate the gradient samples of a phase-encoding dephaser
%--------------------------------------------------------------------------
g_vertex_phase_encoding_dephaser = zeros(4, nr_shots, 'double'); % [mT/m]
for shot_number = 1:nr_shots
    g_vertex_phase_encoding_dephaser(:,shot_number) = cat(1, 0, G_phase_encoding_dephaser(shot_number), G_phase_encoding_dephaser(shot_number), 0);
end

%% Calculate phase blips [mT/m]
%--------------------------------------------------------------------------
% Calculate the duration of a phase blip
%--------------------------------------------------------------------------
blip_duration = echo_spacing - (num_samples * dwell_time); % [sec]

%--------------------------------------------------------------------------
% Calculate the magnitude of a phase blip [mT/m]
% [rad/m] = [rad/sec/T] * [T/1e3mT] * [mT/m] * [sec]
% dkv = gamma * amplitude * duration / 2
% => amplitude = 2 * dkv / gamma / duration
%--------------------------------------------------------------------------
blip_amplitude = 2 * (2 * pi / encoded_fov(2)) / (gamma * 1e-3) / blip_duration; % [mT/m]

%--------------------------------------------------------------------------
% Calculate phase blips
%--------------------------------------------------------------------------
nr_vertices_blips = 3 * (etl - 1);

t_vertex_blips = zeros(nr_vertices_blips, nr_shots, 'double');
g_vertex_blips = zeros(nr_vertices_blips, nr_shots, 'double');

for shot_number = 1:nr_shots
    %----------------------------------------------------------------------
    % Calculate the range of indices per shot
    %----------------------------------------------------------------------
    shot_range = (1:etl).' + (shot_number - 1) * etl;

    %----------------------------------------------------------------------
    % Get phase encoding line numbers per shot
    %----------------------------------------------------------------------
    phase_encoding_line_number_list_per_shot = phase_encoding_line_number_list(shot_range);

    for eco = 1:etl-1
        %------------------------------------------------------------------
        %    ______                ______
        %   /|    |\              /|    |\
        %  / |    | \            / |    | \
        % *--*----*--*--*----*--*--*----*--*--*----*--*-----
        %          |  \ |    | /         |  \ |    | /
        %          |   \|____|/          |   \|____|/
        %          | _        | _        | _
        %          |/ \       |/ \       |/ \
        % *--------|***-------|***-------|***---------------
        %
        %------------------------------------------------------------------
        index_offset = (eco - 1) * 3;

        %------------------------------------------------------------------
        % Calculate the scale factor
        %------------------------------------------------------------------
        scale_factor = phase_encoding_line_number_list_per_shot(eco + 1) - phase_encoding_line_number_list_per_shot(eco);

        %------------------------------------------------------------------
        % Calculate the time points of a phase blip
        %------------------------------------------------------------------
        t_vertex_blips((1:3) + index_offset, shot_number) = t_s + eco * echo_spacing - blip_duration / 2 + cat(1, 0, blip_duration / 2, blip_duration);

        %------------------------------------------------------------------
        % Calculate the gradient values of a phase blip
        %------------------------------------------------------------------
        g_vertex_blips((1:3) + index_offset, shot_number) = scale_factor * cat(1, 0, blip_amplitude, 0);
    end
end

%% Calculate a list of partition-encoding dephasers [mT/m]
%--------------------------------------------------------------------------
% Calculate the duration of a partition-encoding dephaser
%--------------------------------------------------------------------------
partition_encoding_dephaser_duration = ramp_up_time + flat_top_time + ramp_down_time; % [sec]

%--------------------------------------------------------------------------
% Calculate the area of a partition-encoding dephaser [mT/m*sec]
%--------------------------------------------------------------------------
partition_encoding_dephaser_area = zeros(nr_shots, 1, 'double');
for shot_number = 1:nr_shots
    shot_range = (1:etl).' + (shot_number - 1) * etl;
    partition_encoding_number_list_per_shot = partition_encoding_number_list(shot_range);

    if partition_encoding_number_list_per_shot(1) == 0
        partition_encoding_dephaser_area(shot_number) = 0;
    else
        % [rad/m] / ([rad/sec/T] * [T/1e3mT]) => [mT/m*sec]
        partition_encoding_dephaser_area(shot_number) = partition_areas(partition_encoding_number_list_per_shot(1)) / (gamma * 1e-3);
    end
end

%--------------------------------------------------------------------------
% Calculate the amplitude of a partition-encoding dephaser [mT/m]
% area = amplitude * (dephaser_duration + flat_top_time) / 2
%--------------------------------------------------------------------------
G_partition_encoding_dephaser = 2 * partition_encoding_dephaser_area / (partition_encoding_dephaser_duration + flat_top_time);

%--------------------------------------------------------------------------
% Calculate the time samples of a partition-encoding dephaser [sec]
%--------------------------------------------------------------------------
t_vertex_partition_encoding_dephaser = cat(1, -echo_spacing + t_s, -echo_spacing + t_s + ramp_up_time, t_s - ramp_down_time, t_s);

%--------------------------------------------------------------------------
% Calculate the gradient samples of a partition-encoding dephaser
%--------------------------------------------------------------------------
g_vertex_partition_encoding_dephaser = zeros(4, nr_shots, 'double'); % [mT/m]
for shot_number = 1:nr_shots
    g_vertex_partition_encoding_dephaser(:,shot_number) = cat(1, 0, G_partition_encoding_dephaser(shot_number), G_partition_encoding_dephaser(shot_number), 0);
end

%% Calculata a vertex-based RO gradient waveform [mT/m]
t_vertex_readout_waveform = cat(1, t_vertex_prephasing_gradient_lobe, t_vertex_train);
g_vertex_readout_waveform = cat(1, g_vertex_prephasing_gradient_lobe, g_vertex_train);

[c,ia,ic] = unique(t_vertex_readout_waveform);
t_vertex_readout_waveform = t_vertex_readout_waveform(ia); % 307 x 1
g_vertex_readout_waveform = g_vertex_readout_waveform(ia); % 307 x 1

%% Calculate a vertex-based PE gradient waveform [mT/m]
t_vertex_phase_waveform = zeros(4 + nr_vertices_blips + 1, nr_shots, 'double');
g_vertex_phase_waveform = zeros(4 + nr_vertices_blips + 1, nr_shots, 'double');

for shot_number = 1:nr_shots
    t_vertex_phase_waveform(:,shot_number) = cat(1, t_vertex_phase_encoding_dephaser, t_vertex_blips(:,shot_number), t_vertex_readout_waveform(end));
    g_vertex_phase_waveform(:,shot_number) = cat(1, g_vertex_phase_encoding_dephaser(:,shot_number), g_vertex_blips(:,shot_number), 0);
end

%% Calculate a vertex-based SL gradient waveform [mT/m]
t_vertex_slice_waveform = zeros(4 + 1, nr_shots, 'double');
g_vertex_slice_waveform = zeros(4 + 1, nr_shots, 'double');

for shot_number = 1:nr_shots
    t_vertex_slice_waveform(:,shot_number) = cat(1, t_vertex_partition_encoding_dephaser, t_vertex_readout_waveform(end));
    g_vertex_slice_waveform(:,shot_number) = cat(1, g_vertex_partition_encoding_dephaser(:,shot_number), 0);
end

%% Calculate the total number of gradient samples with a gradient raster time of 1 us (fine sampling!)
grad_raster_time = 1e-6; % [sec]
grad_samples = round(echo_spacing * (1 + etl) / grad_raster_time);

%% Calculate a time axis (GRT) [sec]
t_grt = t_vertex_readout_waveform(1) + (0:grad_samples-1).' * grad_raster_time;

%% Interpolate nominal gradient waveforms in the GCS (GRT) [mT/m] [PE,RO,SL]
g_gcs_grt_nominal = zeros(3, grad_samples, nr_shots, 'double'); % [PE,RO,SL]
for shot_number = 1:nr_shots
    g_gcs_grt_nominal(1,:,shot_number) = interp1(t_vertex_phase_waveform(:,shot_number), g_vertex_phase_waveform(:,shot_number), t_grt, 'linear', 'extrap'); % PE (gv)
    g_gcs_grt_nominal(2,:,shot_number) = interp1(t_vertex_readout_waveform, g_vertex_readout_waveform, t_grt, 'linear', 'extrap'); % RO (gu)
    g_gcs_grt_nominal(3,:,shot_number) = interp1(t_vertex_slice_waveform(:,shot_number), g_vertex_slice_waveform(:,shot_number), t_grt, 'linear', 'extrap'); % SL (gw)
end

%% Calculate nominal k-space trajectories in the GCS (GRT) [rad/m] [PE,RO,SL]
%--------------------------------------------------------------------------
% Numerically integrate the coefficients
% [rad/sec/T] * [mT/m] * [T/1e3 mT] * [sec] => * 1e-3 [rad/m]
%--------------------------------------------------------------------------
k_gcs_grt_nominal = cumsum(gamma * g_gcs_grt_nominal * 1e-3 * grad_raster_time, 2); % 3 x grad_samples x nr_shots [rad/m]

%% Calculate the total number of ADC samples per shot
adc_samples = num_samples * etl;

%% Calculate a time axis for ADC samples (ART) [sec]
t_adc = reshape(bsxfun(@plus, acq_delay_time + ((0:num_samples-1).' + 0.5) * dwell_time, t_s + echo_spacing * (0:etl-1)), [adc_samples 1]);

if 0
figure('Color', 'w', 'Position', [1 640 1239 338]);
hold on;
plot(t_vertex_readout_waveform * 1e3, g_vertex_readout_waveform);
plot(t_adc * 1e3, t_adc * 0, '.');
end

%% Calculate a time axis for static field correction
% ACR phantom: 2D EPI => 'SK\SP'
if strfind(twix.hdr.Dicom.tSequenceVariant, 'SP') % spin echo
    t_sfc = t_adc - TE / 2;
else
    t_sfc = t_adc;
end

%% Interpolate nominal k-space trajectories in the GCS (ART) [rad/m] [PE,RO,SL]
tstart = tic; fprintf('%s: Interpolating nominal k-space trajectories in the GCS (ART)... ', datetime);
k_gcs_adc_nominal = zeros(3, adc_samples, nr_shots, 'double'); % [PE,RO,SL]
for shot_number = 1:nr_shots
    k_gcs_adc_nominal(1,:,shot_number) = interp1(t_grt, k_gcs_grt_nominal(1,:,shot_number), t_adc, 'linear', 'extrap'); % PE
    k_gcs_adc_nominal(2,:,shot_number) = interp1(t_grt, k_gcs_grt_nominal(2,:,shot_number), t_adc, 'linear', 'extrap'); % RO
    k_gcs_adc_nominal(3,:,shot_number) = interp1(t_grt, k_gcs_grt_nominal(3,:,shot_number), t_adc, 'linear', 'extrap'); % SL
end
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Calculate nominal gradient waveforms in the DCS (GRT) [mT/m] [x,y,z]
g_dcs_grt_nominal = pagemtimes(R_gcs2dcs, g_gcs_grt_nominal);

%% Calculate nominal k-space trajectories in the DCS (GRT) [rad/m] [x,y,z]
k_dcs_grt_nominal = pagemtimes(R_gcs2dcs, k_gcs_grt_nominal);

%% Calculate nominal k-space trajectories in the DCS (ADC) [rad/m] [x,y,z]
k_dcs_adc_nominal = pagemtimes(R_gcs2dcs, k_gcs_adc_nominal);

%% Calculate the time courses of phase coefficients from nominal gradient waveforms (grad_samples x Nl x nr_shots) (GRT) [rad/m], [rad/m^2], [rad/m^3]
tstart = tic; fprintf('%s: Calculating the time courses of phase coefficients (nominal)... ', datetime);
k_grt_nominal = zeros(grad_samples, Nl, nr_shots, 'double');
for shot_number = 1:nr_shots
    k_grt_nominal(:,:,shot_number) = calculate_concomitant_field_coefficients(g_dcs_grt_nominal(1,:,shot_number).', g_dcs_grt_nominal(2,:,shot_number).', g_dcs_grt_nominal(3,:,shot_number).', Nl, B0, gamma, grad_raster_time);
end
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Interpolate the time courses of phase coefficients from nominal gradient waveforms (Nk x Nl x nr_shots) (ADC) [rad/m], [rad/m^2], [rad/m^3]
k_adc_nominal = interp1(t_grt, k_grt_nominal, t_adc, 'linear', 'extrap'); % Nk x Nl x nr_shots

%% Calculate the time courses of phase coefficients (grad_samples x Nl x nr_shots) (GRT) [rad/m], [rad/m^2], [rad/m^3]
tstart = tic; fprintf('%s: Calculating the time courses of phase coefficients... ', datetime);
k_grt = zeros(grad_samples, Nl, nr_shots, 'double');
for shot_number = 1:nr_shots
    k_grt(:,:,shot_number) = calculate_concomitant_field_coefficients(g_dcs_grt_nominal(1,:,shot_number).', g_dcs_grt_nominal(2,:,shot_number).', g_dcs_grt_nominal(3,:,shot_number).', Nl, B0, gamma, grad_raster_time);
end
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Interpolate the time courses of phase coefficients (grad_samples x Nl x nr_shots) (ADC) [rad/m], [rad/m^2], [rad/m^3]
k_adc = interp1(t_grt, k_grt, t_adc, 'linear', 'extrap'); % grad_samples x Nl x nr_shots

%% Set the time courses of phase coefficients
if cfc_flag
    k = k_adc;
else
    k = zeros(Nk, Nl, nr_shots, 'double');
end

%% Write a .cfl file
%--------------------------------------------------------------------------
% t_grt (grad_samples x 1)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 't_img_grt');
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, t_grt);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% g_gcs_grt_nominal (3 x grad_samples)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'g_img_gcs_grt_nominal');
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, g_gcs_grt_nominal);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% k_gcs_grt_nominal (3 x grad_samples)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'k_img_gcs_grt_nominal');
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, k_gcs_grt_nominal);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% k_grt (grad_samples x Nl)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'k_img_grt');
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, k_grt);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% t_adc (adc_samples x 1)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 't_img_adc');
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, t_adc);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% k_gcs_adc_nominal (3 x adc_samples)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'k_img_gcs_adc_nominal');
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, k_gcs_adc_nominal);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% k_adc (adc_samples x Nl)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'k_img_adc');
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, k_adc);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Calculate a nonuniform readout k-space trajectory [rad/m]
%--------------------------------------------------------------------------
% Calculate a time axis [sec]
%--------------------------------------------------------------------------
t = acq_delay_time + ((0:num_samples-1).' + 0.5) * dwell_time - echo_spacing / 2; % [sec]

%--------------------------------------------------------------------------
% t in [-t2, -t1)
%--------------------------------------------------------------------------
index_range1 = (t >= -t2) & (t < -t1);
% [rad/sec/T] * [T/1e3mT] * ([mT/m/sec] * [sec]^2 - [mT/m] * [sec]) => [rad/m]
k1 = (gamma * 1e-3) / 2 * (SR * (t(index_range1) + t2).^2 - Gu * (t1 + t2));

%--------------------------------------------------------------------------
% t in [-t1, t1)
%--------------------------------------------------------------------------
index_range2 = (t >= -t1) & (t < t1);
% [rad/sec/T] * [T/1e3mT] * [mT/m] * [sec] => [rad/m]
k2 = (gamma * 1e-3) * Gu * t(index_range2);

%--------------------------------------------------------------------------
% t in [t1, t2)
%--------------------------------------------------------------------------
index_range3 = (t >= t1) & (t < t2);
% [rad/sec/T] * [T/1e3mT] * ([mT/m/sec] * [sec]^2 - [mT/m] * [sec]) => [rad/m]
k3 = (gamma * 1e-3) / 2 * (2 * Gu * t1 + 2 * SR * t2 * (t(index_range3) - t1) - SR * (t(index_range3).^2 - t1^2));

%--------------------------------------------------------------------------
% Combine all segments
%--------------------------------------------------------------------------
ku_nonuniform = cat(1, k1, k2, k3); % [rad/m]

%% Calculate the maximum k-space value [rad/m]
kumax = Gu * (gamma * 1e-3) * ((t1 + t2)  - (ta^2 + tb^2) / (2 * (-t1 + t2))) / 2;

%% Calculate a uniform readout k-space trajectory [rad/m]
dku = (2 * pi) / encoded_fov(1); % [rad/m]
ku_uniform = (-floor(Nkx/2):ceil(Nkx/2)-1).' * dku;

%% Calculate the bandwidth per pixel [Hz/Px]
%--------------------------------------------------------------------------
% [Hz] = number of pixels * readout_os_factor * [Hz/Px]
% dwellTime = 1 / (Nkx * bandwith); % [sec]
%--------------------------------------------------------------------------
bandwidth = 1 / (Nkx * dwell_time); % [Hz/Px]

%% Display a readout gradient lobe
baby_blue = [193 220 243] / 255;
blue      = [0   173 236] / 255;
orange    = [239 173 127] / 255;
green     = [205 235 188] / 255;
yellow    = [253 234 155] / 255;

orange_siemens = [236 99 0] / 255;
green_siemens = [3 153 153] / 255;

red_color = [201 37 31] / 255;
blue_color = [86 120 191] / 255;

color_order{1} = '#1f77b4';
color_order{2} = '#ff7f0e';
color_order{3} = '#2ca02c';
color_order{4} = '#d62728';
color_order{5} = '#9467bd';
color_order{6} = '#8c564b';
color_order{7} = '#e377c2';
color_order{8} = '#7f7f7f';
color_order{9} = '#bcbd22';
color_order{10} = '#17becf';

figure('Color', 'w', 'Position', [2 2 1173 593]);
ax1 = subplot(1,2,1);
hold on;
plot(ax1, (t_vertex_readout_gradient_lobe - echo_spacing / 2) * 1e6, g_vertex_readout_gradient_lobe, '.-', 'Color', color_order{1}, 'LineWidth', 1, 'MarkerSize', 10);
plot(ax1, t * 1e6, t * 0, '.', 'Color', 'r', 'LineWidth', 1, 'MarkerSize', 10);
axis square;
grid on; grid minor;
set(ax1, 'XAxisLocation', 'origin', 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontSize', 10);
ylim([-25 25]);
aa = xlabel(ax1, 'Time [$\mu$sec]', 'Interpreter', 'latex', 'FontSize', 12, 'VerticalAlignment', 'baseline');
ylabel(ax1, 'Amplitude [mT/m]', 'Interpreter', 'latex', 'FontSize', 12);
subtitle('Readout gradient lobe', 'Interpreter', 'latex', 'FontSize', 12);
legend(ax1, 'Readout gradient lobe', 'ADC sample locations', 'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', 10);

text(ax1, ax1.XLim(2) * 1.15, 29.5 + 6, 'Ramp sampling in EPI', 'Color', blue, 'Interpreter', 'latex', 'FontSize', 20, 'HorizontalAlignment', 'center');
text(ax1, ax1.XLim(2) * 1.15, 26.5 + 5.7, sprintf('$$G_{\\mathrm{u}} = %4.2f$$ [mT/m], $$t_{a}/t_{\\mathrm{ramp-up}}/t_{\\mathrm{plateau}}/t_{\\mathrm{ramp-down}}$$ = %3.0f/%2.0f/%3.0f/%2.0f [$$\\mu$$sec]', Gu, acq_delay_time * 1e6, ramp_up_time * 1e6, flat_top_time * 1e6, ramp_down_time * 1e6), 'Color', green_siemens, 'Interpreter', 'latex', 'FontSize', 14, 'HorizontalAlignment', 'center');
text(ax1, ax1.XLim(2) * 1.15, 23.5 + 5.4, sprintf('bandwidth = %3.0f [Hz/Px], ADC samples = %d, dwell time = %3.1f [$\\mu$sec]', bandwidth, num_samples, dwell_time * 1e6), 'Color', green_siemens, 'Interpreter', 'latex', 'FontSize', 14, 'HorizontalAlignment', 'center');

ax2 = subplot(1,2,2);
hold on;
plot(ax2, t * 1e6, ku_nonuniform / (2 * kumax), '-', 'Color', color_order{1}, 'LineWidth', 1);
plot(ax2, t * 1e6, ku_uniform / (2 * kumax), '-', 'Color', color_order{2}, 'LineWidth', 1);
axis square;
grid on; grid minor;
set(ax2, 'XAxisLocation', 'origin', 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontSize', 10);
xlabel('Time [$\mu$sec]', 'Interpreter', 'latex', 'FontSize', 12, 'VerticalAlignment', 'baseline');
ylabel('$k_{u} / k_{u,\mathrm{max}}$ [unitless]', 'Interpreter', 'latex', 'FontSize', 12);
subtitle({'Nonuniform vs uniform readout k-space trajectories'}, 'Interpreter', 'latex', 'FontSize', 12);
legend(ax2, 'Nonuniform readout k-space traj.', 'Uniform readout k-space traj.', 'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 10);
drawnow;

set(ax1, 'Position', [0.1317 0.0223 0.3347 0.8150]);
set(ax2, 'Position', [0.5303 0.0223 0.3347 0.8150]);

[filepath,json_filename,ext] = fileparts(json_file);
export_fig(fullfile(output_path, sprintf('ramp_sampling_epi_%s', json_filename)), '-r300', '-tif', '-c[100, 440, 150, 330]'); % [top,right,bottom,left]
close gcf;

%% Read a Siemens .grad file
grad_filename = sprintf('coeff_%s.grad', tGradientCoil); % Gradient Coil daVinci (r0=0.25m, s.u.), Siemens 0.55T Aera
grad_file = fullfile(grad_file_path, grad_filename);
grad = read_siemens_grad_file(grad_file);

%% Parse a Siemens .grad file
grad_info = grad.info;
R0        = grad.R0;      % radius of a gradient coil [m]
alpha_z   = grad.alpha_z; % spherical harmonic coefficients (cosine) for the z gradient coil
alpha_x   = grad.alpha_x; % spherical harmonic coefficients (cosine) for the x gradient coil
alpha_y   = grad.alpha_y; % spherical harmonic coefficients (cosine) for the y gradient coil
beta_z    = grad.beta_z;  % spherical harmonic coefficients (sine)   for the z gradient coil
beta_x    = grad.beta_x;  % spherical harmonic coefficients (sine)   for the x gradient coil
beta_y    = grad.beta_y;  % spherical harmonic coefficients (sine)   for the y gradient coil

%% Calculate the number of total slices
nr_recons = nr_slices * nr_repetitions;

%% Read navigator data (1 x Nkx x nr_navigators x nr_channels x ones(1,6) x (11)nr_repetitions x ones(1,2) x (14)nr_slices)
%--------------------------------------------------------------------------
% ( 1)not used   x ( 2)readout dimension x ( 3)number of TRs  x ( 4)COIL_DIM
% ( 5)MAPS_DIM   x ( 6)TE_DIM            x ( 7)COEFF_DIM      x ( 8)COEFF2_DIM
% ( 9)ITER_DIM   x (10)CSHIFT_DIM        x (11)TIME_DIM       x (12)TIME2_DIM
% (13)LEVEL_DIM  x (14)SLICE_DIM         x (15)AVG_DIM        x (16)BATCH_DIM
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Calculate the number of navigator data per slice
%--------------------------------------------------------------------------
nav_acq_list = find((nav_data.head.idx.slice == 0) & (nav_data.head.idx.repetition == 0));
nr_navigators = length(nav_acq_list);

%--------------------------------------------------------------------------
% Preallocate navigator data
%--------------------------------------------------------------------------
ksp_img_nav = complex(zeros([Nkx nr_navigators Nc nr_repetitions nr_slices], 'single'));

for idx = 1:nr_recons
    %----------------------------------------------------------------------
    % Get information about the current slice
    %----------------------------------------------------------------------
    [slice_number, repetition_number] = ind2sub([nr_slices nr_repetitions], idx);

    %----------------------------------------------------------------------
    % Get a list of navigator acquisitions for the current slice
    %----------------------------------------------------------------------
    nav_acq_list = find((nav_data.head.idx.slice      == (slice_number - 1)) & ...
                        (nav_data.head.idx.repetition == (repetition_number - 1)));

    %----------------------------------------------------------------------
    % Sort one acquisition at a time
    %----------------------------------------------------------------------
    nr_acquisitions = length(nav_acq_list);
    for i = 1:nr_acquisitions
        tstart = tic; fprintf('%s:(SLC=%2d/%2d)(REP=%2d/%2d)(i=%d/%d) Reading "navigator" k-space data... ', datetime, slice_number, nr_slices, repetition_number, nr_repetitions, i, nr_acquisitions);
        ksp_img_nav(:,i,:,repetition_number,slice_number) = nav_data.data{nav_acq_list(i)};
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end
end

%--------------------------------------------------------------------------
% Reshape navigator data
%--------------------------------------------------------------------------
ksp_img_nav = reshape(ksp_img_nav, [1 Nkx nr_navigators Nc ones(1,6) nr_repetitions ones(1,2) nr_slices]);

%% Write a .cfl file
%--------------------------------------------------------------------------
% ksp_img_nav (1 x Nkx x nr_navigators x nr_channels x ones(1,6) x nr_repetitions x ones(1,2) x nr_slices)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'ksp_img_nav');
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
writecfl(cfl_file, ksp_img_nav);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Perform coil compression
%--------------------------------------------------------------------------
% Usage: cc [-p d] [-M] [-r d:d:d] [-A] [-S] [-G] [-E] <kspace> <coeff|proj_kspace>
%
% Performs coil compression.
%
% -p N    perform compression to N virtual channels
% -M      output compression matrix
% -r S    size of calibration region
% -A      use all data to compute coefficients
% -S      type: SVD
% -G      type: Geometric
% -E      type: ESPIRiT
% -h      help
%--------------------------------------------------------------------------
ksp_nav_file    = strcat(bart_output_path, 'ksp_img_nav');
ksp_nav_cc_file = strcat(bart_output_path, 'ksp_img_nav_cc');
command = sprintf('%s cc -p 1 -S %s %s', bart_command, ksp_nav_file, ksp_nav_cc_file);
tstart = tic; fprintf('%s:[BART] Performing coil compression:\n%s\n', datetime, command);
[status_cc,result_cc] = system(command);
fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

%% Read a .cfl file
cfl_file = fullfile(output_path, 'ksp_img_nav_cc');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
ksp_img_nav_cc = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
ksp_img_nav_cc = reshape(ksp_img_nav_cc, [Nkx nr_navigators 1 nr_repetitions nr_slices]);

%% Create a sampling mask
mask = false(Nkx, Nky, Nkz, 'logical');

for shot_number = 1:nr_shots
    %----------------------------------------------------------------------
    % Calculate the range of indices per shot
    %----------------------------------------------------------------------
    shot_range = (1:etl).' + (shot_number - 1) * etl;

    %----------------------------------------------------------------------
    % Get phase encoding line numbers per shot
    %----------------------------------------------------------------------
    phase_encoding_line_number_list_per_shot = phase_encoding_line_number_list(shot_range);

    %----------------------------------------------------------------------
    % Get partition encoding numbers per shot
    %----------------------------------------------------------------------
    partition_encoding_number_list_per_shot = partition_encoding_number_list(shot_range);

    %----------------------------------------------------------------------
    % Assign 1 to acquired locations
    %----------------------------------------------------------------------
    if sum(partition_encoding_number_list_per_shot(:)) == 0
        mask(:,phase_encoding_line_number_list_per_shot,:) = true;
    else
        mask(:,phase_encoding_line_number_list_per_shot,partition_encoding_number_list_per_shot) = true;
    end
end

%% Calculate the number of total slices
nr_recons = nr_slices * nr_contrasts * nr_phases * nr_repetitions * nr_sets;

%% Prepare k-space data per slice
for idx = 1:nr_recons

    %% Get information about the current slice
    [slice_number, contrast_number, phase_number, repetition_number, set_number] = ind2sub([nr_slices nr_contrasts nr_phases nr_repetitions nr_sets], idx);

%     if ~(slice_number == 15 || slice_number == 16)
%         continue;
%     end

%     if repetition_number ~= 1
%         continue;
%     end

    %% Set the description of a dataset
    description = sprintf('slc%d', slice_number);

    if nr_contrasts > 1
        description = sprintf('%s_eco%d', description, contrast_number);
    end

    if nr_phases > 1
        description = sprintf('%s_pha%d', description, phase_number);
    end

    if nr_repetitions > 1
        description = sprintf('%s_rep%d', description, repetition_number);
    end

    if nr_sets > 1
        description = sprintf('%s_set%d', description, set_number);
    end

    %% Define variables
    %----------------------------------------------------------------------
    % Get a list of imaging acquisitions for the current slice
    %----------------------------------------------------------------------
    img_acq_list = find((img_data.head.idx.slice      == (slice_number - 1))      & ...
                        (img_data.head.idx.contrast   == (contrast_number - 1))   & ...
                        (img_data.head.idx.phase      == (phase_number - 1))      & ...
                        (img_data.head.idx.repetition == (repetition_number - 1)) & ...
                        (img_data.head.idx.set        == (set_number - 1)));

    %----------------------------------------------------------------------
    % Get a list of segments (imaging)
    %----------------------------------------------------------------------
    img_segment_list = double(img_data.head.idx.segment(img_acq_list));

    %----------------------------------------------------------------------
    % Get a list of navigator acquisitions for the current slice
    %----------------------------------------------------------------------
    nav_acq_list = find((nav_data.head.idx.slice      == (slice_number - 1))      & ...
                        (nav_data.head.idx.contrast   == (contrast_number - 1))   & ...
                        (nav_data.head.idx.phase      == (phase_number - 1))      & ...
                        (nav_data.head.idx.repetition == (repetition_number - 1)) & ...
                        (nav_data.head.idx.set        == (set_number - 1)));

    %----------------------------------------------------------------------
    % Get a list of segments (navigator)
    %----------------------------------------------------------------------
    nav_segment_list = double(nav_data.head.idx.segment(nav_acq_list));

    %% Process navigator data
    if nav_segment_list(1) == 0 % forward (segment = 0)
        %------------------------------------------------------------------
        % Get the first and third echo signals, S1+ and S3+
        %------------------------------------------------------------------
        S1_forward = ksp_img_nav_cc(:,1,1,repetition_number,slice_number);
        S3_forward = ksp_img_nav_cc(:,3,1,repetition_number,slice_number);

        %------------------------------------------------------------------
        % Calculate a synthetic positive readout echo S2+
        %------------------------------------------------------------------
        S2_forward = (S1_forward + S3_forward) / 2;

        %------------------------------------------------------------------
        % Get the second echo signal, S2-
        %------------------------------------------------------------------
        S2_reverse = ksp_img_nav_cc(:,2,1,repetition_number,slice_number);

    else % reverse (segment = 1)
        %------------------------------------------------------------------
        % Get the first and third echo signals, S1+ and S3+
        %------------------------------------------------------------------
        S1_reverse = ksp_img_nav_cc(:,1,1,repetition_number,slice_number);
        S3_reverse = ksp_img_nav_cc(:,3,1,repetition_number,slice_number);

        %------------------------------------------------------------------
        % Calculate a synthetic positive readout echo S2+
        %------------------------------------------------------------------
        S2_reverse = (S1_reverse + S3_reverse) / 2;

        %------------------------------------------------------------------
        % Get the second echo signal, S2-
        %------------------------------------------------------------------
        S2_forward = ksp_img_nav_cc(:,2,1,repetition_number,slice_number);
    end

    %----------------------------------------------------------------------
    % Perform spline interpolation
    % forward (segment = 0) vs reverse (segment = 1)
    %----------------------------------------------------------------------
    S2_forward_gridded = interp1(+ku_nonuniform, S2_forward, ku_uniform, 'spline', 0);
    S2_reverse_gridded = interp1(-ku_nonuniform, S2_reverse, ku_uniform, 'spline', 0);

    %----------------------------------------------------------------------
    % Apply forward FFT to move from k-space to image space
    % Siemens: k-space <=> image space
    %----------------------------------------------------------------------
    S2_forward_gridded_proj = 1 / sqrt(Nkx) * fftshift(fft(ifftshift(S2_forward_gridded, 1), [], 1), 1);
    S2_reverse_gridded_proj = 1 / sqrt(Nkx) * fftshift(fft(ifftshift(S2_reverse_gridded, 1), [], 1), 1);

    %----------------------------------------------------------------------
    % Calculate the phase difference between S2+ and S2-
    %----------------------------------------------------------------------
    phase_difference = angle(S2_forward_gridded_proj .* conj(S2_reverse_gridded_proj));

    %----------------------------------------------------------------------
    % Calculate a voxel mask
    %----------------------------------------------------------------------
    voxel_list = find(abs(S2_forward_gridded_proj) > max(abs(S2_forward_gridded_proj)) * 0.2);
    nr_voxels = length(voxel_list);

    %----------------------------------------------------------------------
    % Perform a linear least-squares fit
    % y = dphi1 * x + dphi0 = [x 1][dphi1 dphi0].'
    %----------------------------------------------------------------------
    index_range = (-floor(Nkx/2):ceil(Nkx/2)-1).';

    A_index_only = cat(2, index_range(voxel_list), ones(nr_voxels,1));

    c = A_index_only \ phase_difference(voxel_list,:);

    A = cat(2, index_range, ones(Nkx,1));

    phase_difference_fit = A * c;

    %% Display a phase difference signal
    figure('Color', 'w');
    hold on;
    plot(phase_difference);
    plot(phase_difference_fit);
    grid on; grid minor;
    ylim([-pi pi]);
    xlim([1 Nkx]);
    set(gca, 'Box', 'on', 'TickLabelInterpreter', 'latex');
    xlabel('Sample index', 'Interpreter', 'latex');
    ylabel('Phase [rad]', 'Interpreter', 'latex');
    title({sprintf('Phase difference (imaging data)'), sprintf('SLC/REP = %d/%d, Gridding = %d, %d voxels', slice_number, repetition_number, gridding_flag, nr_voxels)}, 'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 12);
    legend('Raw', 'Linear least-squares fit', 'Interpreter', 'latex', 'FontSize', 10);
    export_fig(fullfile(output_path, sprintf('phase_difference_img_slc%d_rep%d_gridding%d', slice_number, repetition_number, gridding_flag)), '-r300', '-tif', '-c[0, 140, 20, 100]'); % [top,right,bottom,left]
    close gcf;

    %% Get "Cartesian-format" "imaging" k-space data (Nkx x Nky x Nkz x Nc)
    %----------------------------------------------------------------------
    % ( 1)not used   x ( 2)readout dimension x ( 3)number of TRs  x ( 4)COIL_DIM
    % ( 5)MAPS_DIM   x ( 6)TE_DIM            x ( 7)COEFF_DIM      x ( 8)COEFF2_DIM
    % ( 9)ITER_DIM   x (10)CSHIFT_DIM        x (11)TIME_DIM       x (12)TIME2_DIM
    % (13)LEVEL_DIM  x (14)SLICE_DIM         x (15)AVG_DIM        x (16)BATCH_DIM
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Calculate the index range (1st index) of a readout acquisition
    %----------------------------------------------------------------------
    kx_range = (discard_pre+1:number_of_samples-discard_post).' - (center_sample + 1) + floor(Nkx/2) + 1;

    %----------------------------------------------------------------------
    % Preallocate k-space data (Nkx x Nky x Nkz x Nc)
    %----------------------------------------------------------------------
    ksp_img_odd_cartesian  = complex(zeros(Nkx, Nky, Nkz, Nc, 'single'));
    ksp_img_even_cartesian = complex(zeros(Nkx, Nky, Nkz, Nc, 'single'));

    %----------------------------------------------------------------------
    % Sort one acquisition at a time
    %----------------------------------------------------------------------
    nr_acquisitions = length(img_acq_list);

    for i = 1:nr_acquisitions
        tstart = tic; fprintf('%s:(SLC=%2d/%2d)(REP=%2d/%2d)(i=%3d/%3d) Reading "imaging" k-space data... ', datetime, slice_number, nr_slices, repetition_number, nr_repetitions, i, nr_acquisitions);

        %------------------------------------------------------------------
        % Get the phase encoding line number
        %------------------------------------------------------------------
        kspace_encode_step_1 = double(img_data.head.idx.kspace_encode_step_1(img_acq_list(i)));
        phase_encoding_line_number = kspace_encode_step_1 - (kspace_encoding_step_1_maximum - kspace_encoding_step_1_center + 1) + floor(Nky/2) + 1;

        %------------------------------------------------------------------
        % Get the partition encoding number
        %------------------------------------------------------------------
        kspace_encode_step_2 = double(img_data.head.idx.kspace_encode_step_2(img_acq_list(i)));
        partition_encoding_number = kspace_encode_step_2 - (kspace_encoding_step_2_maximum - kspace_encoding_step_2_center + 1) + floor(Nkz/2) + 1;
        if partition_encoding_number == 0, partition_encoding_number = partition_encoding_number + 1; end % For 2D imaging

        %------------------------------------------------------------------
        % Prewhiten k-space data
        %------------------------------------------------------------------
        readout = img_data.data{img_acq_list(i)}; % number_of_samples x Nc
        readout = (inv_L * readout.').';

        %------------------------------------------------------------------
        % Perform spline interpolation
        % forward (segment = 0) vs reverse (segment = 1)
        %------------------------------------------------------------------
        if gridding_flag
            if img_segment_list(i) == 0 % forward (segment = 0)
                readout = interp1(+ku_nonuniform, readout, ku_uniform, 'spline', 0);
            else % reverse (segment = 1)
                readout = interp1(-ku_nonuniform, readout, ku_uniform, 'spline', 0);
            end
        else
            %--------------------------------------------------------------
            % Perform row flipping
            %--------------------------------------------------------------
            if img_segment_list(i) == 1
                readout = flip(readout,1);
            end
        end

        %------------------------------------------------------------------
        % Perform phase correction
        %------------------------------------------------------------------
        if phc_flag
            if img_segment_list(i) == 1 % reverse (segment = 1)
                %----------------------------------------------------------
                % Apply forward FFT to move from k-space to image space
                % Siemens: k-space <=> image space
                %----------------------------------------------------------
                proj = 1 / sqrt(Nkx) * fftshift(fft(ifftshift(readout, 1), [], 1), 1);

                %----------------------------------------------------------
                % Perform odd-even echo phase correction
                %----------------------------------------------------------
                proj = proj .* exp(1j * phase_difference_fit);

                %----------------------------------------------------------
                % Apply inverse FFT to move from image space to k-space
                % Siemens: k-space <=> image space
                %----------------------------------------------------------
                readout = sqrt(Nkx) * fftshift(ifft(ifftshift(proj, 1), [], 1), 1);
            end
        end

        %------------------------------------------------------------------
        % Accumulate k-space
        %------------------------------------------------------------------
        if img_segment_list(i) == img_segment_list(1) % odd
            ksp_img_odd_cartesian(kx_range,phase_encoding_line_number,partition_encoding_number,:) = ksp_img_odd_cartesian(kx_range,phase_encoding_line_number,partition_encoding_number,:) + reshape(readout, [Nkx 1 1 Nc]);
        elseif img_segment_list(i) == img_segment_list(2) % even
            ksp_img_even_cartesian(kx_range,phase_encoding_line_number,partition_encoding_number,:) = ksp_img_even_cartesian(kx_range,phase_encoding_line_number,partition_encoding_number,:) + reshape(readout, [Nkx 1 1 Nc]);
        end
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %% Combine odd and even k-space data
    ksp_img_cartesian = ksp_img_odd_cartesian + ksp_img_even_cartesian;

    %% Calculate the actual slice number for Siemens interleaved multislice imaging
    %----------------------------------------------------------------------
    % slice_number: acquisition slice number
    % actual_slice_number is used to retrieve slice information from TWIX format
    % Q: For coronal? acq_slice_nr = slice_nr;
    %----------------------------------------------------------------------
    if nr_slices > 1 % multi-slice
        if mod(nr_slices,2) == 0 % even
            offset1 = 0;
            offset2 = 1;
        else % odd
            offset1 = 1;
            offset2 = 0;
        end
        if slice_number <= ceil(nr_slices / 2)
            actual_slice_number = 2 * slice_number - offset1;
        else
            actual_slice_number = 2 * (slice_number - ceil(nr_slices / 2)) - offset2;
        end
    else
        actual_slice_number = slice_number;
    end

    %% Get a slice offset in the PCS from Siemens TWIX format
    if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_number}, 'sPosition')
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_number}.sPosition, 'dSag')
            sag_offset_twix = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_number}.sPosition.dSag; % [mm]
        else
            sag_offset_twix = 0; % [mm]
        end
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_number}.sPosition, 'dCor')
            cor_offset_twix = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_number}.sPosition.dCor; % [mm]
        else
            cor_offset_twix = 0; % [mm]
        end
        if isfield(twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_number}.sPosition, 'dTra')
            tra_offset_twix = twix.hdr.MeasYaps.sSliceArray.asSlice{actual_slice_number}.sPosition.dTra; % [mm]
        else
            tra_offset_twix = 0; % [mm]
        end
    else
        sag_offset_twix = 0; % [mm]
        cor_offset_twix = 0; % [mm]
        tra_offset_twix = 0; % [mm]
    end

    %% Get a slice offset of a stack in the PCS from Siemens TWIX format
    pcs_offset = [sag_offset_twix; cor_offset_twix; tra_offset_twix] * 1e-3; % [mm] * [m/1e3mm] => [m]

    %% Add the relative table position
    %pcs_offset = pcs_offset - [relative_table_position.x relative_table_position.y relative_table_position.z].' * 1e-3;

    %% Get a slice offset in the PCS from ISMRMRD format
    sag_offset_ismrmrd = double(img_data.head.position(1,img_acq_list(1))); % [mm]
    cor_offset_ismrmrd = double(img_data.head.position(2,img_acq_list(1))); % [mm]
    tra_offset_ismrmrd = double(img_data.head.position(3,img_acq_list(1))); % [mm]
    pcs_offset_ismrmrd = [sag_offset_ismrmrd; cor_offset_ismrmrd; tra_offset_ismrmrd] * 1e-3; % [mm] * [m/1e3mm] => [m]

    %% Get a rotation matrix from the GCS to the PCS (ISMRMRD format)
    phase_dir = double(img_data.head.phase_dir(:,img_acq_list(1)));
    read_dir  = double(img_data.head.read_dir(:,img_acq_list(1)));
    slice_dir = double(img_data.head.slice_dir(:,img_acq_list(1)));
    R_gcs2pcs_ismrmrd = [phase_dir read_dir slice_dir];

    %% Calculate a slice offset in the DCS [m]
    dcs_offset = R_pcs2dcs * pcs_offset; % 3 x 1

    %% Display slice information
    fprintf('======================= SLICE INFORMATION ========================\n');
    fprintf('main_orientation = %d (SAGITTAL/CORONAL/TRANSVERSAL = 0/1/2)\n', main_orientation);
    fprintf('slice_nr = %d, actual_slice_nr = %d\n', slice_number, actual_slice_number);
    fprintf('dNormalSag = %+g \ndNormalCor = %+g \ndNormalTra = %+g \ndRotAngle = %g [rad]\n', dNormalSag, dNormalCor, dNormalTra, dRotAngle);
    fprintf('phase_sign = %+g, read_sign = %+g\n', phase_sign, read_sign);
    fprintf('---------------------- From Siemens TWIX format ------------------\n');
    fprintf('                   [sag]   %10.5f [mm]\n', sag_offset_twix);
    fprintf('slice offset(PCS): [cor] = %10.5f [mm]\n', cor_offset_twix);
    fprintf('                   [tra]   %10.5f [mm]\n', tra_offset_twix);
    fprintf('---------------------- From ISMRMRD format -----------------------\n');
    fprintf('                   [sag]   %10.5f [mm]\n', sag_offset_ismrmrd);
    fprintf('slice offset(PCS): [cor] = %10.5f [mm]\n', cor_offset_ismrmrd);
    fprintf('                   [tra]   %10.5f [mm]\n', tra_offset_ismrmrd);
    fprintf('-------------- From Siemens TWIX format (w/ gradient flips) ------\n');
    fprintf('                   [sag]   [%10.5f %10.5f %10.5f][PE]\n', R_gcs2pcs(1,1), R_gcs2pcs(1,2), R_gcs2pcs(1,3));
    fprintf('R_gcs2pcs        : [cor] = [%10.5f %10.5f %10.5f][RO]\n', R_gcs2pcs(2,1), R_gcs2pcs(2,2), R_gcs2pcs(2,3));
    fprintf('                   [tra]   [%10.5f %10.5f %10.5f][SL]\n', R_gcs2pcs(3,1), R_gcs2pcs(3,2), R_gcs2pcs(3,3));
    fprintf('-------------- From ISMRMRD format (w/o gradient flips) ----------\n');
    fprintf('                   [sag]   [%10.5f %10.5f %10.5f][PE]\n', R_gcs2pcs_ismrmrd(1,1), R_gcs2pcs_ismrmrd(1,2), R_gcs2pcs_ismrmrd(1,3));
    fprintf('R_gcs2pcs        : [cor] = [%10.5f %10.5f %10.5f][RO]\n', R_gcs2pcs_ismrmrd(2,1), R_gcs2pcs_ismrmrd(2,2), R_gcs2pcs_ismrmrd(2,3));
    fprintf('                   [tra]   [%10.5f %10.5f %10.5f][SL]\n', R_gcs2pcs_ismrmrd(3,1), R_gcs2pcs_ismrmrd(3,2), R_gcs2pcs_ismrmrd(3,3));
    fprintf('------------------------------------------------------------------\n');
    fprintf('                   [ x ]   [%10.5f %10.5f %10.5f][sag]\n', R_pcs2dcs(1,1), R_pcs2dcs(1,2), R_pcs2dcs(1,3));
    fprintf('R_pcs2dcs        : [ y ] = [%10.5f %10.5f %10.5f][cor]\n', R_pcs2dcs(2,1), R_pcs2dcs(2,2), R_pcs2dcs(2,3));
    fprintf('                   [ z ]   [%10.5f %10.5f %10.5f][tra]\n', R_pcs2dcs(3,1), R_pcs2dcs(3,2), R_pcs2dcs(3,3));
    fprintf('------------------------------------------------------------------\n');
    fprintf('                   [ x ]   [%10.5f %10.5f %10.5f][PE]\n', R_gcs2dcs(1,1), R_gcs2dcs(1,2), R_gcs2dcs(1,3));
    fprintf('R_gcs2dcs        : [ y ] = [%10.5f %10.5f %10.5f][RO]\n', R_gcs2dcs(2,1), R_gcs2dcs(2,2), R_gcs2dcs(2,3));
    fprintf('                   [ z ]   [%10.5f %10.5f %10.5f][SL]\n', R_gcs2dcs(3,1), R_gcs2dcs(3,2), R_gcs2dcs(3,3));
    fprintf('==================================================================\n');

    %% Read a .cfl file
    %----------------------------------------------------------------------
    % x (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('x_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    x = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % y (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('y_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    y = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % z (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('z_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    z = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % u (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('u_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    u = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % v (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('v_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    v = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % w (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('w_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    w = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Calculate the concomitant field phase
    %----------------------------------------------------------------------
    % Calculate the amplitude of a trapezoid
    %----------------------------------------------------------------------
    Gx = g_dcs(1);
    Gy = g_dcs(2);
    Gz = g_dcs(3);

    %----------------------------------------------------------------------
    % Calculate the rise time of a trapezoid
    %----------------------------------------------------------------------
    r = ramp_up_time;

    %----------------------------------------------------------------------
    % Calculate the duration of a trapezoid
    %----------------------------------------------------------------------
    w = echo_spacing;

    %----------------------------------------------------------------------
    % Calculate the concomitant field phase
    %----------------------------------------------------------------------
    % [rad/sec/T] / [T] [m]^2 * [mT/m]^2 * [T/1e3mT]^2 * [sec] => [rad]
    conc_phase1 = gamma / (2 * B0) * ((Gx * 1e-3) * z - (Gz * 1e-3) * x / 2).^2 * ((w - 2 * r) + 2 * r / 3);
    conc_phase2 = gamma / (2 * B0) * ((Gy * 1e-3) * z - (Gz * 1e-3) * y / 2).^2 * ((w - 2 * r) + 2 * r / 3);
    conc_phase = conc_phase1 + conc_phase2;

    %% Calculate a parabolic shift [m]
    cfl_file = fullfile(output_path, sprintf('parabolic_shift_slc%d_%s', slice_number, slice_type));
    if ~exist(strcat(cfl_file, '.cfl'), 'file')
        tstart = tic; fprintf('%s:(SLC=%2d/%2d) Calculating a parabolic shift... ', datetime, slice_number, nr_slices);
        %------------------------------------------------------------------
        % Calculate a sampling interval in the phase-encoding direction [rad/m]
        %------------------------------------------------------------------
        dkv = (2 * pi) / encoded_fov(2); % [rad/m]

        %------------------------------------------------------------------
        % Calculate the number of interleaves per kx-ky plane
        %------------------------------------------------------------------
        nr_interleaves = 1; % single-shot EPI

        %------------------------------------------------------------------
        % Calculate a parabolic shift [m]
        % dk_blip = nr_interleaves * acceleration_factor1 * dky [rad/m]
        %------------------------------------------------------------------
        parabolic_shift = conc_phase / (nr_interleaves * acceleration_factor1 * dkv); % [rad] / [rad/m] => [m]
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %% Write a .cfl file
    %----------------------------------------------------------------------
    % parabolic_shift (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('parabolic_shift_slc%d_%s', slice_number, slice_type));
    if ~exist(strcat(cfl_file, '.cfl'), 'file')
        tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
        writecfl(cfl_file, parabolic_shift);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %% Calculate a circular mask or a cube mask
    %----------------------------------------------------------------------
    % The correction is limited to a 2 * r0 cube enclosing the specified FOV.
    % Pixels outside the 2 * r0 FOV are set to zero for safety reasons.
    %----------------------------------------------------------------------
    circle_mask = zeros(Nkx, Nky, Nkz, 'single');
    %circle_mask((x.^2 + y.^2 + z.^2) < R0^2) = 1;
    circle_mask((abs(x) < R0) & (abs(y) < R0) & (abs(z) < R0)) = 1;

    %% Convert Cartesian coordinates to spherical coordinates
    cfl_file = fullfile(output_path, sprintf('dx_slc%d_%s', slice_number, slice_type));
    if ~exist(strcat(cfl_file, '.cfl'), 'file')
        tstart = tic; fprintf('%s:(SLC=%2d/%2d) Converting Cartesian coordinates (x,y,z) to spherical coordinates (r,theta,phi)... ', datetime, slice_number, nr_slices);
        radius = reshape(sqrt(x.^2 + y.^2 + z.^2), [N 1]); % N x 1 [m]
        theta = reshape(atan2(sqrt(x.^2 + y.^2), z), [N 1]); % same as acos(z./r)
        phi = reshape(atan2(y, x), [N 1]); % azimuth
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %% Calculate a displacement field along the x-axis in the DCS [m]
    cfl_file = fullfile(output_path, sprintf('dx_slc%d_%s', slice_number, slice_type));
    if ~exist(strcat(cfl_file, '.cfl'), 'file')
        tstart = tic; fprintf('%s:(SLC=%2d/%2d) Calculating a displacement field along the x-axis... ', datetime, slice_number, nr_slices);
        dx = reshape(siemens_B(radius, theta, phi, R0, alpha_x, beta_x), [Nkx Nky Nkz]); % N x 1 [m]
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %% Calculate a displacement field along the y-axis in the DCS [m]
    cfl_file = fullfile(output_path, sprintf('dx_slc%d_%s', slice_number, slice_type));
    if ~exist(strcat(cfl_file, '.cfl'), 'file')
        tstart = tic; fprintf('%s:(SLC=%2d/%2d) Calculating a displacement field along the y-axis... ', datetime, slice_number, nr_slices);
        dy = reshape(siemens_B(radius, theta, phi, R0, alpha_y, beta_y), [Nkx Nky Nkz]); % N x 1 [m]
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %% Calculate a displacement field along the z-axis in the DCS [m]
    cfl_file = fullfile(output_path, sprintf('dx_slc%d_%s', slice_number, slice_type));
    if ~exist(strcat(cfl_file, '.cfl'), 'file')
        tstart = tic; fprintf('%s:(SLC=%2d/%2d) Calculating a displacement field along the z-axis... ', datetime, slice_number, nr_slices);
        dz = reshape(siemens_B(radius, theta, phi, R0, alpha_z, beta_z), [Nkx Nky Nkz]); % N x 1 [m]
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %% Calculate displacement fields in the GCS [m] (3 x N) [v,u,w] = [PE,RO,SL]
    cfl_file = fullfile(output_path, sprintf('dx_slc%d_%s', slice_number, slice_type));
    if ~exist(strcat(cfl_file, '.cfl'), 'file')
        dr_gcs = R_gcs2dcs.' * cat(2, dx(:), dy(:), dz(:)).'; % 3 x N
        du = reshape(dr_gcs(2,:), [Nkx Nky Nkz]); % RO [m]
        dv = reshape(dr_gcs(1,:), [Nkx Nky Nkz]); % PE [m]
        dw = reshape(dr_gcs(3,:), [Nkx Nky Nkz]); % SL [m]
    end

    %% Read a .cfl file
    if sfc_flag
        %------------------------------------------------------------------
        % fieldmap (Nkx x Nky x Nkz)
        %------------------------------------------------------------------
        cfl_file = fullfile(output_path, sprintf('fieldmap_smooth_slc%d_gnc%d_%s', slice_number, gnc_flag, slice_type));
        tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
        fieldmap = real(readcfl(cfl_file));
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    else
        fieldmap = zeros(Nkx, Nky, Nkz, 'single');
    end

    %% Calculate the accumulated phase due to static off-resonance (phase right after a refocusing pulse)
    if strfind(twix.hdr.Dicom.tSequenceVariant, 'SP') % 'SK\SP' (spin echo ?)
        initial_phase = exp(-1j * 2 * pi * fieldmap * TE / 2);
    else
        initial_phase = complex(zeros(Nkx, Nky, Nkz, 'single'));
    end

    %% Calculate concomitant field basis functions (N x Nl) [m], [m^2], [m^3]
    cfl_file = fullfile(output_path, sprintf('U_img_slc%d_%s_cfc%d_sfc%d', slice_number, slice_type, cfc_flag, sfc_flag));
    if ~exist(strcat(cfl_file, '.cfl'), 'file')
        tstart = tic; fprintf('%s:(SLC=%2d/%2d) Calculating concomitant field basis functions... ', datetime, slice_number, nr_slices);
        p = calculate_concomitant_field_basis(x(:), y(:), z(:), Nl);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %% Calculate the SVD of a higher-order encoding matrix (Nk x N)
    if cfc_flag
        cfl_file = fullfile(output_path, sprintf('U_img_slc%d_%s_cfc%d_sfc%d', slice_number, slice_type, cfc_flag, sfc_flag));
        if ~exist(strcat(cfl_file, '.cfl'), 'file')
            tstart = tic; fprintf('%s:(SLC=%2d/%2d) Calculating randomized SVD... ', datetime, slice_number, nr_slices);
            [u_tilde,s_tilde,v_tilde] = calculate_rsvd_higher_order_encoding_matrix(k(:,4:Nl), p(:,4:Nl), Lmax, os, fieldmap, t_sfc, sfc_flag);
            U = single(u_tilde(:,1:Lmax)); % Nk x Lmax
            S = diag(single(s_tilde(1:Lmax,1:Lmax))); % Lmax x Lmax => Lmax x 1
            V = reshape(single(v_tilde(:,1:Lmax) * s_tilde(1:Lmax,1:Lmax)'), [Nkx Nky Nkz Lmax]); % Nkx x Nky x Nkz x Lmax
            fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
        else
            cfl_file = fullfile(output_path, sprintf('S_img_slc%d_%s_cfc%d_sfc%d', slice_number, slice_type, cfc_flag, sfc_flag));
            S = readcfl(cfl_file);
        end
    else
        U = ones(Nk, Lmax, 'single');
        S = cat(1, 1, zeros(Lmax - 1, 1, 'single'));
        V = ones(Nkx, Nky, Nkz, Lmax, 'single');
    end

    S_scaled = S / S(1);
    cumulative_sum = cumsum(S_scaled) / sum(S_scaled);
    indices = find(cumulative_sum > cutoff_percentage * 1e-2);
    L = indices(1);

    %% Display the singular values and cumulative sum
    if repetition_number == 1
        FontSize = 12;

        figure('Color', 'w', 'Position', [3 533 1048 455]);
        subplot(1,2,1);
        plot(S_scaled, '.-', 'LineWidth', 1, 'MarkerSize', 12, 'Color', color_order{1});
        grid on;
        grid minor;
        set(gca, 'YScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', FontSize);
        axis square;
        xlabel('$$\ell$$', 'Interpreter', 'latex', 'FontSize', FontSize + 2);
        ylabel('Singular values, $$\sigma_{\ell}$$', 'Interpreter', 'latex', 'FontSize', FontSize + 2);

        subplot(1,2,2);
        hold on;
        plot(cumulative_sum, '.-', 'LineWidth', 1, 'MarkerSize', 12, 'Color', color_order{1});
        plot([L L], [0 1.1], '--', 'Color', color_order{2}, 'LineWidth', 1);
        axis square;
        grid on;
        grid minor;
        set(gca, 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontSize', FontSize);
        ylim([0 1.1]);
        legend('', sprintf('%5.2f$$\\%%$$ cutoff, L = %d', cutoff_percentage, L), 'Location', 'Southwest', 'Interpreter', 'latex', 'FontSize', FontSize);
        xlabel('$$\ell$$', 'Interpreter', 'latex', 'FontSize', FontSize + 2);
        ylabel('Cumulative sum, $$(\Sigma_{k=1}^{\ell}\sigma_{k}) / (\Sigma_{k=1}^{L_{\mathrm{max}}}\sigma_{k})$$', 'Interpreter', 'latex', 'FontSize', FontSize + 2);
        export_fig(fullfile(output_path, sprintf('singular_values_and_cumulative_sum_slc%d_%s', slice_number, slice_type)), '-r300', '-tif'); % [top,right,bottom,left]
        close gcf;
    end

    %% Write a .cfl file
    %----------------------------------------------------------------------
    % L (1 x 1)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('L_slc%d_%s', slice_number, slice_type));
    if ~exist(strcat(cfl_file, '.cfl'), 'file')
        tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
        writecfl(cfl_file, L);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %% Calculate time-reversed U
    cfl_file = fullfile(output_path, sprintf('U_img_slc%d_%s_cfc%d_sfc%d', slice_number, slice_type, cfc_flag, sfc_flag));
    if ~exist(strcat(cfl_file, '.cfl'), 'file')
        tstart = tic; fprintf('%s:(SLC=%2d/%2d) Calculating time-reversed U... ', datetime, slice_number, nr_slices);
        U = reshape(U, [Nkx etl Lmax]);
        U_time_reversed = complex(zeros(Nkx, Nky, Nkz, Lmax, 'single'));

        for shot_number = 1:nr_shots
            %--------------------------------------------------------------
            % Calculate the range of indices per shot
            %--------------------------------------------------------------
            shot_range = (1:etl).' + (shot_number - 1) * etl;

            %--------------------------------------------------------------
            % Get phase encoding line numbers per shot
            %--------------------------------------------------------------
            phase_encoding_line_number_list_per_shot = phase_encoding_line_number_list(shot_range);

            %--------------------------------------------------------------
            % Get partition encoding numbers per shot
            %--------------------------------------------------------------
            partition_encoding_number_list_per_shot = partition_encoding_number_list(shot_range);

            %--------------------------------------------------------------
            % Get a list of segments per shot
            %--------------------------------------------------------------
            img_segment_list_per_shot = img_segment_list(shot_range);

            for eco = 1:etl
                %----------------------------------------------------------
                % Get the phase encoding line number
                %----------------------------------------------------------
                phase_encoding_line_number = phase_encoding_line_number_list_per_shot(eco);

                %----------------------------------------------------------
                % Get the partition encoding number
                %----------------------------------------------------------
                partition_encoding_number = partition_encoding_number_list_per_shot(eco);
                if partition_encoding_number == 0, partition_encoding_number = partition_encoding_number + 1; end

                %----------------------------------------------------------
                % Flip a "reverse" line in the left singular vectors
                %----------------------------------------------------------
                if img_segment_list_per_shot(eco) == 0 % forward (segment == 0)
                    U_time_reversed(:,phase_encoding_line_number,partition_encoding_number,:) = U(:,eco,:);
                else % reverse (segment == 1)
                    U_time_reversed(:,phase_encoding_line_number,partition_encoding_number,:) = flip(U(:,eco,:),1);
                end
            end
        end

        %------------------------------------------------------------------
        % Reorder U
        %------------------------------------------------------------------
        U = reshape(U, [Nk Lmax]);
        U_time_reversed = reshape(U_time_reversed, [Nkx * Nky * Nkz Lmax]);
        U_time_reversed = U_time_reversed((mask > 0),:);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %% Write a .cfl file
    %----------------------------------------------------------------------
    % ksp_img_cartesian (Nkx x Nky x Nkz x Nc)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('ksp_img_cartesian_slc%d_rep%d_gridding%d_phc%d', slice_number, repetition_number, gridding_flag, phc_flag));
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, ksp_img_cartesian);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % mask (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('mask_img_slc%d', slice_number));
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, mask);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % circle_mask (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('circle_mask_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, circle_mask);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % dx (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('dx_slc%d_%s', slice_number, slice_type));
    if ~exist(strcat(cfl_file, '.cfl'), 'file')
        tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
        writecfl(cfl_file, dx);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %----------------------------------------------------------------------
    % dy (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('dy_slc%d_%s', slice_number, slice_type));
    if ~exist(strcat(cfl_file, '.cfl'), 'file')
        tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
        writecfl(cfl_file, dy);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %----------------------------------------------------------------------
    % dz (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('dz_slc%d_%s', slice_number, slice_type));
    if ~exist(strcat(cfl_file, '.cfl'), 'file')
        tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
        writecfl(cfl_file, dz);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %----------------------------------------------------------------------
    % du (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('du_slc%d_%s', slice_number, slice_type));
    if ~exist(strcat(cfl_file, '.cfl'), 'file')
        tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
        writecfl(cfl_file, du);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %----------------------------------------------------------------------
    % dv (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('dv_slc%d_%s', slice_number, slice_type));
    if ~exist(strcat(cfl_file, '.cfl'), 'file')
        tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
        writecfl(cfl_file, dv);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %----------------------------------------------------------------------
    % dw (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('dw_slc%d_%s', slice_number, slice_type));
    if ~exist(strcat(cfl_file, '.cfl'), 'file')
        tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
        writecfl(cfl_file, dw);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %----------------------------------------------------------------------
    % initial_phase (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('initial_phase_slc%d_%s', slice_number, slice_type));
    if ~exist(strcat(cfl_file, '.cfl'), 'file')
        tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
        writecfl(cfl_file, initial_phase);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %----------------------------------------------------------------------
    % U (Nk * nr_shots x Lmax)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('U_img_slc%d_%s_cfc%d_sfc%d', slice_number, slice_type, cfc_flag, sfc_flag));
    if ~exist(strcat(cfl_file, '.cfl'), 'file')
        tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
        writecfl(cfl_file, U_time_reversed);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %----------------------------------------------------------------------
    % S (Lmax x 1)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('S_img_slc%d_%s_cfc%d_sfc%d', slice_number, slice_type, cfc_flag, sfc_flag));
    if ~exist(strcat(cfl_file, '.cfl'), 'file')
        tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
        writecfl(cfl_file, S);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end

    %----------------------------------------------------------------------
    % V (Nkx x Nky x Nkz x Lmax)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('V_img_slc%d_%s_cfc%d_sfc%d', slice_number, slice_type, cfc_flag, sfc_flag));
    if ~exist(strcat(cfl_file, '.cfl'), 'file')
        tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
        writecfl(cfl_file, V);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end
end
