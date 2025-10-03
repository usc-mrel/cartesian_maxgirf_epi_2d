% demo_figure3.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 02/10/2025, Last modified: 06/29/2025

%% Clean slate
close all; clear all; clc;

%% Start a stopwatch timer
start_time = tic;

%% Set source directories
package_path = 'D:\cartesian_maxgirf_epi_2d';
ismrmrd_path = 'D:\ismrmrd';

%% Add source directories to search path
addpath(genpath(package_path));
addpath(genpath(ismrmrd_path));

%% Define the full path of an output directory
% PE direction: A >> P
json_file1 = 'D:\cartesian_maxgirf_epi_2d\data\acr_phantom_20250209\json_files\meas_MID00125_FID08928_ep2d_se_bw976_tra_AP_avg16_fov256_256_pf_R2_gridding1_phc1_cfc1_sfc0_gnc0_topup0.json';
json_file2 = 'D:\cartesian_maxgirf_epi_2d\data\acr_phantom_20250209\json_files\meas_MID00125_FID08928_ep2d_se_bw976_tra_AP_avg16_fov256_256_pf_R2_gridding1_phc1_cfc1_sfc0_gnc1_topup0.json';
json_file3 = 'D:\cartesian_maxgirf_epi_2d\data\acr_phantom_20250209\json_files\meas_MID00125_FID08928_ep2d_se_bw976_tra_AP_avg16_fov256_256_pf_R2_gridding1_phc1_cfc1_sfc0_gnc1_topup1.json';

dicom_path1 = 'D:\cartesian_maxgirf_epi_2d\data\acr_phantom_20250209\dicom\S6_ep2d_se_bw976_tra_AP_avg16_fov256_256_pf_R2'; % ND
dicom_path2 = 'D:\cartesian_maxgirf_epi_2d\data\acr_phantom_20250209\dicom\S24_ep2d_se_bw976_tra_AP_avg16_fov256_256_pf_R2_S6_DIS2D'; % DIS2D

nii_applytopup_img1_cor_file = 'D:\cartesian_maxgirf_epi_2d\data\acr_phantom_20250209\topup_AP_PA\applytopup_img1_cor.nii.gz';

%% Get directory information
dir_info = dir(fullfile(dicom_path1, '*IMA'));
nr_files = length(dir_info);

%% Get a DICOM header
dicom_info = dicominfo(fullfile(dir_info(1).folder, dir_info(1).name));

%% Parse the DICOM header
%--------------------------------------------------------------------------
% Rows Attribute
% Number of rows in the image
%--------------------------------------------------------------------------
N1 = double(dicom_info.Rows);

%--------------------------------------------------------------------------
% Columns  Attribute
% Number of columns in the image
%--------------------------------------------------------------------------
N2 = double(dicom_info.Columns);

%% Read a .dicom file
img1_dicom = zeros(N1, N2, nr_files, 'single');

x1_dicom = zeros(N1, N2, nr_files, 'single');
y1_dicom = zeros(N1, N2, nr_files, 'single');
z1_dicom = zeros(N1, N2, nr_files, 'single');

slice_location = zeros(nr_files, 1, 'double');

for idx = 1:nr_files

    %% Get a DICOM header
    dicom_info = dicominfo(fullfile(dir_info(idx).folder, dir_info(idx).name));

    %% Read a DICOM image
    if isfield(dicom_info, 'RescaleSlope')
        RescaleSlope = dicom_info.RescaleSlope;
    else
        RescaleSlope = 1;
    end

    if isfield(dicom_info, 'RescaleIntercept')
        RescaleIntercept = dicom_info.RescaleIntercept;
    else
        RescaleIntercept = 0;
    end

    img1_dicom(:,:,idx) = RescaleSlope * double(dicomread(dicom_info)).' + RescaleIntercept; % transpose it!

    %% Parse the DICOM header
    %----------------------------------------------------------------------
    % Patient Position Attribute
    % Patient position descriptor relative to the equipment
    %----------------------------------------------------------------------
    patient_position = dicom_info.PatientPosition;

    %----------------------------------------------------------------------
    % Slice Thickness Attribute
    % Nominal slice thickness, in mm
    %----------------------------------------------------------------------
    slice_thickness = dicom_info.SliceThickness; % [mm]

    %----------------------------------------------------------------------
    % Image Position (Patient) Attribute
    % The x, y, and z coordinates of the upper left hand corner
    % (center of the first voxel transmitted) of the image, in mm
    %----------------------------------------------------------------------
    ipp = dicom_info.ImagePositionPatient; % [mm]

    %----------------------------------------------------------------------
    % Image Orientation (Patient) Attribute
    % The direction cosines of the first row and the first column with respect
    % to the patient
    %----------------------------------------------------------------------
    iop = dicom_info.ImageOrientationPatient;

    %----------------------------------------------------------------------
    % Pixel Spacing Attribute
    % Physical distance in the Patient between the center of each pixel, specified
    % by a numeric pair - adjacent row spacing, adjacent column spacing in mm
    %----------------------------------------------------------------------
    pixel_spacing = dicom_info.PixelSpacing; % [mm]

    %----------------------------------------------------------------------
    % Number of slices
    %----------------------------------------------------------------------
    N3 = 1;

    %----------------------------------------------------------------------
    % Slice number
    %----------------------------------------------------------------------
    instance_number = dicom_info.InstanceNumber;

    %----------------------------------------------------------------------
    % Acquisition number
    %----------------------------------------------------------------------
    acquisition_number = dicom_info.AcquisitionNumber;

    %----------------------------------------------------------------------
    % Slice location [mm]
    %----------------------------------------------------------------------
    slice_location(idx) = dicom_info.SliceLocation;
    fprintf('idx=%2d, instance number = %2d, acquisition number = %2d, SliceLocation = %6.2f [mm]\n', idx, instance_number, acquisition_number, dicom_info.SliceLocation);

    %% Calculate the total number of voxels
    N = N1 * N2 * N3;

    %% Calculate a rotation matrix from the PCS to the DCS
    R_pcs2dcs = siemens_calculate_transform_pcs_to_dcs(patient_position);

    %% Calculate a scaling matrix
    scaling_matrix_dicom = [pixel_spacing(1) 0 0; 0 pixel_spacing(2) 0; 0 0 slice_thickness] * 1e-3; % [mm] * [m/1e3mm] => [m]

    %% Calculate a trannsformation matrix from the RCS to the PCS [r,c,s] <=> [L,P,S]
    R_rcs2pcs = cat(2, iop(1:3), iop(4:6), cross(iop(1:3), iop(4:6)));

    %% Calculate spatial coordinates in the RCS [m]
    [I1,I2,I3] = ndgrid((0:N1-1).', (0:N2-1).', (0:N3-1).');
    r_rcs = scaling_matrix_dicom * cat(1, I1(:).', I2(:).', I3(:).'); % 3 x N

    %% Calculate spatial coordinates in the LPH [m] (R => L, A => P, I => S)
    % The DICOM LPH coordinate system is identical to the Siemens Patient Coordinate System (PCS)
    r_pcs = repmat(ipp * 1e-3, [1 N]) + R_rcs2pcs * r_rcs;

    %% Add a table position
    r_pcs = r_pcs + repmat(-dicom_info.Private_0019_1014 * 1e-3, [1 N]);

    %% Calculate spatial coordinates in the DCS [m]
    r_dcs = R_pcs2dcs * r_pcs;

    %% Save arrays
    x1_dicom(:,:,idx) = reshape(r_dcs(1,:), [N1 N2]); % N x 1 [m]
    y1_dicom(:,:,idx) = reshape(r_dcs(2,:), [N1 N2]); % N x 1 [m]
    z1_dicom(:,:,idx) = reshape(r_dcs(3,:), [N1 N2]); % N x 1 [m]
end

%% Sort the dicom images
[sorted_values, sorted_indices] = sort(slice_location);
img1_dicom = img1_dicom(:,:,sorted_indices);

%% Get directory information
dir_info = dir(fullfile(dicom_path2, '*IMA'));
nr_files = length(dir_info);

%% Get a DICOM header
dicom_info = dicominfo(fullfile(dir_info(1).folder, dir_info(1).name));

%% Parse the DICOM header
%--------------------------------------------------------------------------
% Rows Attribute
% Number of rows in the image
%--------------------------------------------------------------------------
N1 = double(dicom_info.Rows);

%--------------------------------------------------------------------------
% Columns  Attribute
% Number of columns in the image
%--------------------------------------------------------------------------
N2 = double(dicom_info.Columns);

%% Read a .dicom file
img2_dicom = zeros(N1, N2, nr_files, 'single');

x2_dicom = zeros(N1, N2, nr_files, 'single');
y2_dicom = zeros(N1, N2, nr_files, 'single');
z2_dicom = zeros(N1, N2, nr_files, 'single');

slice_location = zeros(nr_files, 1, 'double');

for idx = 1:nr_files

    %% Get a DICOM header
    dicom_info = dicominfo(fullfile(dir_info(idx).folder, dir_info(idx).name));

    %% Read a DICOM image
    if isfield(dicom_info, 'RescaleSlope')
        RescaleSlope = dicom_info.RescaleSlope;
    else
        RescaleSlope = 1;
    end

    if isfield(dicom_info, 'RescaleIntercept')
        RescaleIntercept = dicom_info.RescaleIntercept;
    else
        RescaleIntercept = 0;
    end

    img2_dicom(:,:,idx) = RescaleSlope * double(dicomread(dicom_info)).' + RescaleIntercept; % transpose it!

    %% Parse the DICOM header
    %----------------------------------------------------------------------
    % Patient Position Attribute
    % Patient position descriptor relative to the equipment
    %----------------------------------------------------------------------
    patient_position = dicom_info.PatientPosition;

    %----------------------------------------------------------------------
    % Slice Thickness Attribute
    % Nominal slice thickness, in mm
    %----------------------------------------------------------------------
    slice_thickness = dicom_info.SliceThickness; % [mm]

    %----------------------------------------------------------------------
    % Image Position (Patient) Attribute
    % The x, y, and z coordinates of the upper left hand corner
    % (center of the first voxel transmitted) of the image, in mm
    %----------------------------------------------------------------------
    ipp = dicom_info.ImagePositionPatient; % [mm]

    %----------------------------------------------------------------------
    % Image Orientation (Patient) Attribute
    % The direction cosines of the first row and the first column with respect
    % to the patient
    %----------------------------------------------------------------------
    iop = dicom_info.ImageOrientationPatient;

    %----------------------------------------------------------------------
    % Pixel Spacing Attribute
    % Physical distance in the Patient between the center of each pixel, specified
    % by a numeric pair - adjacent row spacing, adjacent column spacing in mm
    %----------------------------------------------------------------------
    pixel_spacing = dicom_info.PixelSpacing; % [mm]

    %----------------------------------------------------------------------
    % Number of slices
    %----------------------------------------------------------------------
    N3 = 1;

    %----------------------------------------------------------------------
    % Slice number
    %----------------------------------------------------------------------
    instance_number = dicom_info.InstanceNumber;

    %----------------------------------------------------------------------
    % Acquisition number
    %----------------------------------------------------------------------
    acquisition_number = dicom_info.AcquisitionNumber;

    %----------------------------------------------------------------------
    % Slice location [mm]
    %----------------------------------------------------------------------
    slice_location(idx) = dicom_info.SliceLocation;
    fprintf('idx=%2d, instance number = %2d, acquisition number = %2d, SliceLocation = %6.2f [mm]\n', idx, instance_number, acquisition_number, dicom_info.SliceLocation);

    %% Calculate the total number of voxels
    N = N1 * N2 * N3;

    %% Calculate a rotation matrix from the PCS to the DCS
    R_pcs2dcs = siemens_calculate_transform_pcs_to_dcs(patient_position);

    %% Calculate a scaling matrix
    scaling_matrix_dicom = [pixel_spacing(1) 0 0; 0 pixel_spacing(2) 0; 0 0 slice_thickness] * 1e-3; % [mm] * [m/1e3mm] => [m]

    %% Calculate a trannsformation matrix from the RCS to the PCS [r,c,s] <=> [L,P,S]
    R_rcs2pcs = cat(2, iop(1:3), iop(4:6), cross(iop(1:3), iop(4:6)));

    %% Calculate spatial coordinates in the RCS [m]
    [I1,I2,I3] = ndgrid((0:N1-1).', (0:N2-1).', (0:N3-1).');
    r_rcs = scaling_matrix_dicom * cat(1, I1(:).', I2(:).', I3(:).'); % 3 x N

    %% Calculate spatial coordinates in the LPH [m] (R => L, A => P, I => S)
    % The DICOM LPH coordinate system is identical to the Siemens Patient Coordinate System (PCS)
    r_pcs = repmat(ipp * 1e-3, [1 N]) + R_rcs2pcs * r_rcs;

    %% Add a table position
    r_pcs = r_pcs + repmat(-dicom_info.Private_0019_1014 * 1e-3, [1 N]);

    %% Calculate spatial coordinates in the DCS [m]
    r_dcs = R_pcs2dcs * r_pcs;

    %% Save arrays
    x2_dicom(:,:,idx) = reshape(r_dcs(1,:), [N1 N2]); % N x 1 [m]
    y2_dicom(:,:,idx) = reshape(r_dcs(2,:), [N1 N2]); % N x 1 [m]
    z2_dicom(:,:,idx) = reshape(r_dcs(3,:), [N1 N2]); % N x 1 [m]
end

%% Sort the dicom images
[sorted_values, sorted_indices] = sort(slice_location);
img2_dicom = img2_dicom(:,:,sorted_indices);

%% Read a .json file
tstart = tic; fprintf('%s: Reading a .json file: %s... ', datetime, json_file1);
fid = fopen(json_file1);
json_txt = fread(fid, [1 inf], 'char=>char');
fclose(fid);
json = jsondecode(json_txt);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Define the full path of a filename
%--------------------------------------------------------------------------
output_path1       = strrep(json.output_path, '/', '\');
ismrmrd_data_file1 = strrep(json.ismrmrd_data_file, '/', '\');

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
Lmax          = json.recon_parameters.Lmax;           % maximum rank of the SVD approximation of a higher-order encoding matrix
lambda        = json.recon_parameters.lambda;         % l2 regularization parameter
tol           = json.recon_parameters.tol;            % PCG tolerance
maxiter       = json.recon_parameters.maxiter;        % PCG maximum iteration 
slice_type    = json.recon_parameters.slice_type;     % type of an excitation slice: "curved" vs "flat"
phc_flag      = json.recon_parameters.phc_flag;       % 1=yes, 0=no
gridding_flag = json.recon_parameters.gridding_flag;  % 1=yes, 0=no
cfc_flag      = json.recon_parameters.cfc_flag;       % 1=yes, 0=no
sfc_flag      = json.recon_parameters.sfc_flag;       % 1=yes, 0=no
gnc_flag      = json.recon_parameters.gnc_flag;       % 1=yes, 0=no
topup_flag    = json.recon_parameters.topup_flag;     % 1=yes, 0=no

%--------------------------------------------------------------------------
% Number of slices
%--------------------------------------------------------------------------
if isfield(json, 'nr_slices')
    nr_slices = json.nr_slices;
else
    nr_slices = 1;
end

%--------------------------------------------------------------------------
% Number of repetitions
%--------------------------------------------------------------------------
if isfield(json, 'nr_repetitions')
    nr_repetitions = json.nr_repetitions;
else
    nr_repetitions = 1;
end

%--------------------------------------------------------------------------
% main_orientation (SAGITTAL/CORONAL/TRANSVERSAL = 0/1/2)
%--------------------------------------------------------------------------
if isfield(json, 'main_orientation')
    main_orientation = json.main_orientation;
else
    main_orientation = 2;
end

%% Make an output path
pclr_output_path1 = fullfile(output_path1, 'PCLR');

%% Read an ISMRMRD file (k-space data)
tstart = tic; fprintf('%s: Reading an ISMRMRD file: %s... ', datetime, ismrmrd_data_file1);
if exist(ismrmrd_data_file1, 'file')
    dset = ismrmrd.Dataset(ismrmrd_data_file1, 'dataset');
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
else
    error('File %s does not exist.  Please generate it.' , ismrmrd_data_file1);
end

%% Get imaging parameters from an XML header
header = ismrmrd.xml.deserialize(dset.readxml);

%--------------------------------------------------------------------------
% Encoding Space (Nkx, Nky, Nkz)
%--------------------------------------------------------------------------
encoded_fov(1) = header.encoding.encodedSpace.fieldOfView_mm.x * 1e-3; % [m] RO
encoded_fov(2) = header.encoding.encodedSpace.fieldOfView_mm.y * 1e-3; % [m] PE
encoded_fov(3) = header.encoding.encodedSpace.fieldOfView_mm.z * 1e-3; % [m] SL

Nkx = header.encoding.encodedSpace.matrixSize.x; % number of readout samples in k-space
Nky = header.encoding.encodedSpace.matrixSize.y; % number of phase-encoding steps in k-space
Nkz = header.encoding.encodedSpace.matrixSize.z; % number of slice-encoding steps in k-space

%% Read a .cfl file
img1 = complex(zeros(Nkx, Nky, nr_slices, 'single'));

x = complex(zeros(Nkx, Nky, nr_slices, 'single'));
y = complex(zeros(Nkx, Nky, nr_slices, 'single'));
z = complex(zeros(Nkx, Nky, nr_slices, 'single'));

dx = complex(zeros(Nkx, Nky, nr_slices, 'single'));
dy = complex(zeros(Nkx, Nky, nr_slices, 'single'));
dz = complex(zeros(Nkx, Nky, nr_slices, 'single'));

for slice_number = 13%1:nr_slices
    %----------------------------------------------------------------------
    % Calculate the actual slice number for Siemens interleaved multislice imaging
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

    %----------------------------------------------------------------------
    % img (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    img_filename = sprintf('img_maxgirf_pclr_slc%d_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_topup%d_i%d_l%4.2f', slice_number, slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, topup_flag, maxiter, lambda);
    cfl_file = fullfile(pclr_output_path1, img_filename);
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices, cfl_file);
    img1(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % x (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path1, sprintf('x_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices, cfl_file);
    x(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    
    %----------------------------------------------------------------------
    % y (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path1, sprintf('y_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices, cfl_file);
    y(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % z (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path1, sprintf('z_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices, cfl_file);
    z(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % dx (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path1, sprintf('dx_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices, cfl_file);
    dx(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % dy (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path1, sprintf('dy_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices, cfl_file);
    dy(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % dz (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path1, sprintf('dz_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices, cfl_file);
    dz(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%--------------------------------------------------------------------------
% read_sign
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path1, 'read_sign');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
read_sign = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% phase_sign
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path1, 'phase_sign');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
phase_sign = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Read a .json file
tstart = tic; fprintf('%s: Reading a .json file: %s... ', datetime, json_file2);
fid = fopen(json_file2);
json_txt = fread(fid, [1 inf], 'char=>char');
fclose(fid);
json = jsondecode(json_txt);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Define the full path of a filename
%--------------------------------------------------------------------------
output_path2 = strrep(json.output_path, '/', '\');

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
Lmax          = json.recon_parameters.Lmax;           % maximum rank of the SVD approximation of a higher-order encoding matrix
lambda        = json.recon_parameters.lambda;         % l2 regularization parameter
tol           = json.recon_parameters.tol;            % PCG tolerance
maxiter       = json.recon_parameters.maxiter;        % PCG maximum iteration 
slice_type    = json.recon_parameters.slice_type;     % type of an excitation slice: "curved" vs "flat"
phc_flag      = json.recon_parameters.phc_flag;       % 1=yes, 0=no
gridding_flag = json.recon_parameters.gridding_flag;  % 1=yes, 0=no
cfc_flag      = json.recon_parameters.cfc_flag;       % 1=yes, 0=no
sfc_flag      = json.recon_parameters.sfc_flag;       % 1=yes, 0=no
gnc_flag      = json.recon_parameters.gnc_flag;       % 1=yes, 0=no
topup_flag    = json.recon_parameters.topup_flag;     % 1=yes, 0=no

%% Make an output path
pclr_output_path2 = fullfile(output_path2, 'PCLR');

%% Read a .cfl file
img2 = complex(zeros(Nkx, Nky, nr_slices, 'single'));

for slice_number = 13%1:nr_slices
    %----------------------------------------------------------------------
    % Calculate the actual slice number for Siemens interleaved multislice imaging
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

    %----------------------------------------------------------------------
    % img (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    img_filename = sprintf('img_maxgirf_pclr_slc%d_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_topup%d_i%d_l%4.2f', slice_number, slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, topup_flag, maxiter, lambda);
    cfl_file = fullfile(pclr_output_path2, img_filename);
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices, cfl_file);
    img2(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Read a .json file
tstart = tic; fprintf('%s: Reading a .json file: %s... ', datetime, json_file3);
fid = fopen(json_file3);
json_txt = fread(fid, [1 inf], 'char=>char');
fclose(fid);
json = jsondecode(json_txt);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% Define the full path of a filename
%--------------------------------------------------------------------------
output_path3 = strrep(json.output_path, '/', '\');
topup_path3 = strrep(json.topup_path, '/', '\');

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
Lmax          = json.recon_parameters.Lmax;           % maximum rank of the SVD approximation of a higher-order encoding matrix
lambda        = json.recon_parameters.lambda;         % l2 regularization parameter
tol           = json.recon_parameters.tol;            % PCG tolerance
maxiter       = json.recon_parameters.maxiter;        % PCG maximum iteration 
slice_type    = json.recon_parameters.slice_type;     % type of an excitation slice: "curved" vs "flat"
phc_flag      = json.recon_parameters.phc_flag;       % 1=yes, 0=no
gridding_flag = json.recon_parameters.gridding_flag;  % 1=yes, 0=no
cfc_flag      = json.recon_parameters.cfc_flag;       % 1=yes, 0=no
sfc_flag      = json.recon_parameters.sfc_flag;       % 1=yes, 0=no
gnc_flag      = json.recon_parameters.gnc_flag;       % 1=yes, 0=no
topup_flag    = json.recon_parameters.topup_flag;     % 1=yes, 0=no

%% Make an output path
pclr_output_path3 = fullfile(output_path3, 'PCLR');

%% Read a .cfl file
img3 = complex(zeros(Nkx, Nky, nr_slices, 'single'));

for slice_number = 13%1:nr_slices
    %----------------------------------------------------------------------
    % Calculate the actual slice number for Siemens interleaved multislice imaging
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

    %----------------------------------------------------------------------
    % img (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    img_filename = sprintf('img_maxgirf_pclr_slc%d_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_topup%d_i%d_l%4.2f', slice_number, slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, topup_flag, maxiter, lambda);
    cfl_file = fullfile(pclr_output_path3, img_filename);
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices, cfl_file);
    img3(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Set variables
[Nkx,Nky,~] = size(img1);

Nx = Nkx / 2;
Ny = Nky;
Nz = 1;

N = Nx * Ny * Nz;

%% Read a .nii file
img3_dicom = single(niftiread(nii_applytopup_img1_cor_file));
img3_dicom = flip(img3_dicom,2);

%% Read a .cfl file
displacement = zeros(Nx, Ny, nr_slices, 'single');

for slice_number = 1:nr_slices
    %----------------------------------------------------------------------
    % Calculate the actual slice number for Siemens interleaved multislice imaging
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
    
    %----------------------------------------------------------------------
    % displacement (Nx x Ny)
    %----------------------------------------------------------------------
    img_filename = sprintf('displacement_slc%d_flat', slice_number);
    cfl_file = fullfile(topup_path3, img_filename);
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices, cfl_file);
    displacement(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

if read_sign < 0
    displacement = flip(displacement,1);
end

if phase_sign < 0
    displacement = flip(displacement,2);
end

%% Adjust the image size of a custom reconstruction
% -128:127 => flip => 127:-128 => crop => 
idx1_range = (-floor(Nx/2):ceil(Nx/2)-1).' + floor(Nkx/2) + 1;
idx2_range = (-floor(Ny/2):ceil(Ny/2)-1).' + floor(Nky/2) + 1;
idx3_range = (-floor(Nz/2):ceil(Nz/2)-1).' + floor(Nkz/2) + 1;

img1 = img1(idx1_range, idx2_range, :);
img2 = img2(idx1_range, idx2_range, :);
img3 = img3(idx1_range, idx2_range, :);

x = x(idx1_range, idx2_range, :);
y = y(idx1_range, idx2_range, :);
z = z(idx1_range, idx2_range, :);

dx = dx(idx1_range, idx2_range, :);
dy = dy(idx1_range, idx2_range, :);
dz = dz(idx1_range, idx2_range, :);

%% Flip variables
if read_sign < 0
    img1 = flip(img1,1);
    img2 = flip(img2,1);
    img3 = flip(img3,1);

    x = flip(x,1);
    y = flip(y,1);
    z = flip(z,1);

    dx = flip(dx,1);
    dy = flip(dy,1);
    dz = flip(dz,1);

    displacement = flip(displacement,1);
end

if phase_sign < 0
    img1 = flip(img1,2);
    img2 = flip(img2,2);
    img3 = flip(img3,2);

    x = flip(x,2);
    y = flip(y,2);
    z = flip(z,2);

    dx = flip(dx,2);
    dy = flip(dy,2);
    dz = flip(dz,2);

    displacement = flip(displacement,2);
end

%% Scale images
c1 = floor(Nx/2) + 1;
c2 = floor(Ny/2) + 1;

slice_number = 26;

scale_factor1 = abs(img1(c1,c2,slice_number));
scale_factor2 = abs(img2(c1,c2,slice_number));
scale_factor3 = abs(img3(c1,c2,slice_number));

img1_scaled = img1 / scale_factor1;
img2_scaled = img2 / scale_factor2;
img3_scaled = img3 / scale_factor3;

scale_factor1_dicom = abs(img1_dicom(c1,c2,slice_number));
scale_factor2_dicom = abs(img2_dicom(c1,c2,slice_number));
scale_factor3_dicom = abs(img3_dicom(c1,c2,slice_number));

img1_dicom_scaled = img1_dicom / scale_factor1_dicom;
img2_dicom_scaled = img2_dicom / scale_factor2_dicom;
img3_dicom_scaled = img3_dicom / scale_factor3_dicom;

%% Display images
close all;

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

cmap = flip(brewermap([],"RdBu"),1);

FontSize = 12;

idx1_range = (17:238).'; % left-right
idx2_range = (30:251).'; % up-down

N1 = length(idx1_range);
N2 = length(idx2_range);

idx1_range_zoom = (153:208).'; % left-right
idx2_range_zoom = (166:221).'; % up-down

N1_zoom = length(idx1_range_zoom);
N2_zoom = length(idx2_range_zoom);

slice_number = 26;

climits = [0 2.5];

%--------------------------------------------------------------------------
% Calculate a circle
%--------------------------------------------------------------------------
circle_radius = 95;
t = (0:0.01:1).';
x_circle = circle_radius * cos(2 * pi * t) + floor(N1/2) + 1 - 1;
y_circle = circle_radius * sin(2 * pi * t) + floor(N1/2) + 1;

figure('Color', 'w', 'Position', [0 31 1128 947]);

%--------------------------------------------------------------------------
% Cartesian MaxGIRF: Gridding/PHC/CFC/SFC/GNC/TOPUP = 1/1/1/0/0/0
%--------------------------------------------------------------------------
ax1 = subplot(2,4,1);
hold on;
imagesc(ax1, abs(img1_scaled(idx1_range,idx2_range,slice_number)).');
plot(x_circle, y_circle, '--', 'Color', 'r', 'LineWidth', 1);

% left, vertical (up-down)
plot(ax1, [idx1_range_zoom(1) - idx1_range(1) idx1_range_zoom(1) - idx1_range(1)], ...
          [idx2_range_zoom(1) - idx2_range(1) 194], 'Color', red_color, 'LineWidth', 1.5);
% [136 136], [136 194]

% right, vertical (up-down)
plot(ax1, [136 136] + 55, [136 194], 'Color', red_color, 'LineWidth', 1.5);

% top, horizotal (left-right)
plot(ax1, [136 136 + 55], [136 136], 'Color', red_color, 'LineWidth', 1.5);

% top, horizotal (left-right)
plot(ax1, [136 136 + 55], [194 194], 'Color', red_color, 'LineWidth', 1.5);

axis image ij off;
colormap(ax1, gray(256));
clim(ax1, climits);
title(ax1, '(Phase-const.) Cartesian MaxGIRF', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subtitle(ax1, {sprintf('Gridding/PHC/{\\color[rgb]{%f %f %f}CFC}', orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3)), 'complex averaging'}, 'Rotation', 0, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax1, 2, 0, '(A)', 'FontSize', 14, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
text(ax1, N1, N2 / 2, sprintf('PE direction (A >> P)'), 'FontSize', FontSize, 'Rotation', -90, 'Interpreter', 'tex', 'Color', 'r', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

%--------------------------------------------------------------------------
% Cartesian MaxGIRF: Gridding/PHC/CFC/SFC/GNC/TOPUP = 1/1/1/0/1/0
%--------------------------------------------------------------------------
ax2 = subplot(2,4,2);
hold on;
imagesc(ax2, abs(img2_scaled(idx1_range,idx2_range,slice_number)).');
plot(x_circle, y_circle, '--', 'Color', 'r', 'LineWidth', 1);
axis image ij off;
colormap(ax2, gray(256));
clim(ax2, climits);
title(ax2, '(Phase-const.) Cartesian MaxGIRF', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subtitle(ax2, {sprintf('Gridding/PHC/{\\color[rgb]{%f %f %f}CFC}/{\\color[rgb]{%f %f %f}GNC}', orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3)), 'complex averaging'}, 'Rotation', 0, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax2, 2, 0, '(B)', 'FontSize', 14, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% title
%--------------------------------------------------------------------------
text(ax2, N1, -90 - 16, sprintf('2D SE-EPI: axial, 1.0 x 1.0 mm^2, R = 2, PF = 6/8, ETL = 95, 16 NSA, z = %3.1f mm', z(1,1,slice_number) * 1e3), 'Color', blue, 'Interpreter', 'tex', 'FontSize', 20, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

text(ax2, N1, -38 - 16, {sprintf('Gridding/PHC: gridding for ramp sampling/odd-even echo phase correction'), ...
                    sprintf('{\\color[rgb]{%f %f %f}CFC}/{\\color[rgb]{%f %f %f}GNC} vs {\\color[rgb]{%f %f %f}CFC}/{\\color[rgb]{%f %f %f}GNC(2D)}: concomitant field correction/gradient nonlinearity correction {\\color[rgb]{%f %f %f}during recon} vs {\\color[rgb]{%f %f %f}after recon}', ...
                    orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), green_siemens(1,1), green_siemens(1,2), green_siemens(1,3), green_siemens(1,1), green_siemens(1,2), green_siemens(1,3), ...
                    orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), green_siemens(1,1), green_siemens(1,2), green_siemens(1,3)), ...
                    sprintf('{\\color[rgb]{%f %f %f}TOPUP} vs {\\color[rgb]{%f %f %f}TOPUP}: susceptibility-induced distortion correction using a displacement map from FSP TOPUP {\\color[rgb]{%f %f %f}during recon} vs {\\color[rgb]{%f %f %f}after recon}', ...
                    orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3),  green_siemens(1,1), green_siemens(1,2), green_siemens(1,3), ...
                    orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3),  green_siemens(1,1), green_siemens(1,2), green_siemens(1,3))}, ...
                    'Rotation', 0, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');


% Create line (top)
annotation(gcf, 'line', [0.0063 0.9920], [0.9110+0.021 0.9110+0.021], 'LineWidth', 2);

% Create line (bottom)
annotation(gcf, 'line', [0.0063 0.9920], [0.8478+0.016 0.8478+0.016], 'LineWidth', 2);

%--------------------------------------------------------------------------
% Cartesian MaxGIRF: Gridding/PHC/CFC/SFC/GNC/TOPUP = 1/1/1/0/1/1
%--------------------------------------------------------------------------
ax3 = subplot(2,4,3);
hold on;
imagesc(ax3, abs(img3_scaled(idx1_range,idx2_range,slice_number)).');
plot(x_circle, y_circle, '--', 'Color', 'r', 'LineWidth', 1);
axis image ij off;
colormap(ax3, gray(256));
clim(ax3, climits);
title(ax3, '(Phase-const.) Cartesian MaxGIRF', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subtitle(ax3, {sprintf('Gridding/PHC/{\\color[rgb]{%f %f %f}CFC}/{\\color[rgb]{%f %f %f}GNC}/{\\color[rgb]{%f %f %f}TOPUP}', ...
    orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3)), 'complex averaging'}, 'Rotation', 0, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax3, 2, 0, '(C)', 'FontSize', 14, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Displacement field along the x-axis
%--------------------------------------------------------------------------
ax4 = subplot(2,4,4);
hold on;
imagesc(ax4, dx(idx1_range,idx2_range,slice_number).' * 1e3);
contour(ax4, dx(idx1_range,idx2_range,slice_number).' * 1e3, (-3:0.5:3).', 'ShowText' ,'on', 'LevelStep', 4, 'LineWidth', 1, 'Color', 'k');

% left, vertical (up-down)
plot(ax4, [idx1_range_zoom(1) - idx1_range(1) idx1_range_zoom(1) - idx1_range(1)], ...
          [idx2_range_zoom(1) - idx2_range(1) 194], 'Color', red_color, 'LineWidth', 1.5);
% [136 136], [136 194]

% right, vertical (up-down)
plot(ax4, [136 136] + 55, [136 194], 'Color', red_color, 'LineWidth', 1.5);

% top, horizotal (left-right)
plot(ax4, [136 136 + 55], [136 136], 'Color', red_color, 'LineWidth', 1.5);

% top, horizotal (left-right)
plot(ax4, [136 136 + 55], [194 194], 'Color', red_color, 'LineWidth', 1.5);

axis image ij off;
colormap(ax4, cmap);
clim(ax4, [-3 3]);
title(ax4, 'Displacement field', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subtitle(ax4, {'along the x-axis (RO direction)'}, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax4, N1 / 2 - 3, 0, {'x-axis $$\longrightarrow$$'}, 'Rotation', 0, 'Color', 'k', 'Interpreter', 'latex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(ax4, 2, 0, '(D)', 'FontSize', 14, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
hc4 = colorbar;
set(hc4, 'Position', [0.9645 0.5321 0.0097 0.2378], 'FontSize', FontSize);
hTitle4 = title(hc4, '[mm]', 'FontSize', FontSize, 'Position', [12.4914 172.5567 0]);

%--------------------------------------------------------------------------
% DICOM: Gridding/PHC/CFC/SFC/GNC/TOPUP = 1/1/1/0/0/0
%--------------------------------------------------------------------------
ax5 = subplot(2,4,5);
hold on;
imagesc(ax5, abs(img1_dicom_scaled(idx1_range,idx2_range,slice_number)).');
plot(x_circle, y_circle, '--', 'Color', 'r', 'LineWidth', 1);
axis image ij off;
colormap(ax5, gray(256));
clim(ax5, climits);
title(ax5, 'Traditional reconstruction', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subtitle(ax5, {sprintf('Gridding/PHC/{\\color[rgb]{%f %f %f}CFC}', green_siemens(1,1), green_siemens(1,2), green_siemens(1,3)), 'magnitude averaging'}, 'Rotation', 0, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax5, 2, 0, '(E)', 'FontSize', 14, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% DICOM: Gridding/PHC/CFC/SFC/GNC/TOPUP = 1/1/1/0/1/0
%--------------------------------------------------------------------------
ax6 = subplot(2,4,6);
hold on;
imagesc(ax6, abs(img2_dicom_scaled(idx1_range,idx2_range,slice_number)).');
plot(x_circle, y_circle, '--', 'Color', 'r', 'LineWidth', 1);
axis image ij off;
colormap(ax6, gray(256));
clim(ax6, climits);
title(ax6, 'Traditional reconstruction', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subtitle(ax6, {sprintf('Gridding/PHC/{\\color[rgb]{%f %f %f}CFC}/{\\color[rgb]{%f %f %f}GNC}', green_siemens(1,1), green_siemens(1,2), green_siemens(1,3), green_siemens(1,1), green_siemens(1,2), green_siemens(1,3)), 'magnitude averaging'}, 'Rotation', 0, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax6, 2, 0, '(F)', 'FontSize', 14, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% DICOM: Gridding/PHC/CFC/SFC/GNC/TOPUP = 1/1/1/0/1/1
%--------------------------------------------------------------------------
ax7 = subplot(2,4,7);
hold on;
imagesc(ax7, abs(img3_dicom_scaled(idx1_range,idx2_range,slice_number)).');
plot(x_circle, y_circle, '--', 'Color', 'r', 'LineWidth', 1);
axis image ij off;
colormap(ax7, gray(256));
clim(ax7, climits);
title(ax7, 'Traditional reconstruction', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subtitle(ax7, {sprintf('Gridding/PHC/{\\color[rgb]{%f %f %f}CFC}/{\\color[rgb]{%f %f %f}GNC}/{\\color[rgb]{%f %f %f}TOPUP}', ...
    green_siemens(1,1), green_siemens(1,2), green_siemens(1,3), green_siemens(1,1), green_siemens(1,2), green_siemens(1,3), green_siemens(1,1), green_siemens(1,2), green_siemens(1,3)), 'magnitude averaging'}, 'Rotation', 0, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax7, 2, 0, '(G)', 'FontSize', 14, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Displacement field along the y-axis
%--------------------------------------------------------------------------
ax8 = subplot(2,4,8);
hold on;
imagesc(ax8, dy(idx1_range,idx2_range,slice_number).' * 1e3);
contour(ax8, dy(idx1_range,idx2_range,slice_number).' * 1e3, (-3:0.5:3).', 'ShowText' ,'on', 'LevelStep', 4, 'LineWidth', 1, 'Color', 'k');

% left, vertical (up-down)
plot(ax8, [idx1_range_zoom(1) - idx1_range(1) idx1_range_zoom(1) - idx1_range(1)], ...
          [idx2_range_zoom(1) - idx2_range(1) 194], 'Color', red_color, 'LineWidth', 1.5);
% [136 136], [136 194]

% right, vertical (up-down)
plot(ax8, [136 136] + 55, [136 194], 'Color', red_color, 'LineWidth', 1.5);

% top, horizotal (left-right)
plot(ax8, [136 136 + 55], [136 136], 'Color', red_color, 'LineWidth', 1.5);

% top, horizotal (left-right)
plot(ax8, [136 136 + 55], [194 194], 'Color', red_color, 'LineWidth', 1.5);

axis image ij off;
colormap(ax8, cmap);
clim(ax8, [-2.2 2.2]);
title(ax8, 'Displacement field', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
subtitle(ax8, {'along the y-axis (PE direction)'}, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax8, 0, N2 / 2 - 6, {'y-axis $$\longrightarrow$$'}, 'Rotation', 90, 'Color', 'k', 'Interpreter', 'latex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(ax8, 2, 0, '(H)', 'FontSize', 14, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
hc8 = colorbar;
set(hc8, 'Position', [0.9645 0.2047 - 0.021 0.0097 0.2378], 'FontSize', FontSize);
hTitle8 = title(hc8, '[mm]', 'FontSize', FontSize, 'Position', [12.4914 172.5567 0]);

%--------------------------------------------------------------------------
% Cartesian MaxGIRF: Zoom 1
%--------------------------------------------------------------------------
ax9 = axes;
imagesc(ax9, abs(img1_scaled(idx1_range_zoom,idx2_range_zoom,slice_number)).');
axis image;
colormap(ax9, gray(256));
clim(ax9, [0 1.5]);
set(ax9, 'XColor', red_color, 'YColor', red_color, 'XTickLabel', [], 'YTickLabel', [],  'XTick', [], 'YTick', [], 'LineWidth', 2);

%--------------------------------------------------------------------------
% Cartesian MaxGIRF: Zoom 2
%--------------------------------------------------------------------------
ax10 = axes;
imagesc(ax10, abs(img2_scaled(idx1_range_zoom,idx2_range_zoom,slice_number)).');
axis image;
colormap(ax10, gray(256));
clim(ax10, [0 1.5]);
set(ax10, 'XColor', red_color, 'YColor', red_color, 'XTickLabel', [], 'YTickLabel', [],  'XTick', [], 'YTick', [], 'LineWidth', 2);

%--------------------------------------------------------------------------
% Cartesian MaxGIRF: Zoom 3
%--------------------------------------------------------------------------
ax11 = axes;
imagesc(ax11, abs(img3_scaled(idx1_range_zoom,idx2_range_zoom,slice_number)).');
axis image;
colormap(ax11, gray(256));
clim(ax11, [0 1.5]);
set(ax11, 'XColor', red_color, 'YColor', red_color, 'XTickLabel', [], 'YTickLabel', [],  'XTick', [], 'YTick', [], 'LineWidth', 2);

%--------------------------------------------------------------------------
% DICOM: Zoom 1
%--------------------------------------------------------------------------
ax12 = axes;
imagesc(ax12, abs(img1_dicom_scaled(idx1_range_zoom,idx2_range_zoom,slice_number)).');
axis image;
colormap(ax12, gray(256));
clim(ax12, [0 1.5]);
set(ax12, 'XColor', red_color, 'YColor', red_color, 'XTickLabel', [], 'YTickLabel', [],  'XTick', [], 'YTick', [], 'LineWidth', 2);

%--------------------------------------------------------------------------
% DICOM: Zoom 2
%--------------------------------------------------------------------------
ax13 = axes;
imagesc(ax13, abs(img2_dicom_scaled(idx1_range_zoom,idx2_range_zoom,slice_number)).');
axis image;
colormap(ax13, gray(256));
clim(ax13, [0 1.5]);
set(ax13, 'XColor', red_color, 'YColor', red_color, 'XTickLabel', [], 'YTickLabel', [],  'XTick', [], 'YTick', [], 'LineWidth', 2);

%--------------------------------------------------------------------------
% DICOM: Zoom 3
%--------------------------------------------------------------------------
ax14 = axes;
imagesc(ax14, abs(img3_dicom_scaled(idx1_range_zoom,idx2_range_zoom,slice_number)).');
axis image;
colormap(ax14, gray(256));
clim(ax14, [0 1.5]);
set(ax14, 'XColor', red_color, 'YColor', red_color, 'XTickLabel', [], 'YTickLabel', [],  'XTick', [], 'YTick', [], 'LineWidth', 2);

set(ax1, 'Position', [0.0396-0.035 0.4604 0.2379 0.3896]);
set(ax2, 'Position', [0.2790-0.035 0.4604 0.2379 0.3896]);
set(ax3, 'Position', [0.5184-0.035 0.4604 0.2379 0.3896]);
set(ax4, 'Position', [0.7578-0.035 0.4604 0.2379 0.3896]);

set(ax5, 'Position', [0.0396-0.035 0.1325-0.021 0.2379 0.3896]);
set(ax6, 'Position', [0.2790-0.035 0.1325-0.021 0.2379 0.3896]);
set(ax7, 'Position', [0.5184-0.035 0.1325-0.021 0.2379 0.3896]);
set(ax8, 'Position', [0.7578-0.035 0.1325-0.021 0.2379 0.3896]);

set(ax9 , 'Position', [-0.0204 + 0.2394 * 0 0.5147 0.1594 0.1273]);
set(ax10, 'Position', [-0.0204 + 0.2394 * 1 0.5147 0.1594 0.1273]);
set(ax11, 'Position', [-0.0204 + 0.2394 * 2 0.5147 0.1594 0.1273]);

set(ax12, 'Position', [-0.0204 + 0.2394 * 0 0.5147 - 0.3279 - 0.021 0.1594 0.1273]);
set(ax13, 'Position', [-0.0204 + 0.2394 * 1 0.5147 - 0.3279 - 0.021 0.1594 0.1273]);
set(ax14, 'Position', [-0.0204 + 0.2394 * 2 0.5147 - 0.3279 - 0.021 0.1594 0.1273]);

export_fig(sprintf('figure3_slc%d', slice_number), '-r300', '-tif', '-c[70, 10, 470, 0]'); % [top,right,bottom,left]
