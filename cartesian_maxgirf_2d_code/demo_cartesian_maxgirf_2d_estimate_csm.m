% demo_cartesian_maxgirf_2d_estimate_csm.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 01/18/2025, Last modified: 02/08/2025

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
    ismrmrd_data_file = strrep(json.ismrmrd_data_file, '/', '\');
    output_path       = strrep(json.output_path, '/', '\');
    topup_path        = strrep(json.topup_path, '/', '\');
else
    ismrmrd_data_file = json.ismrmrd_data_file;
    output_path       = json.output_path;
    topup_path        = json.topup_path;
end

%--------------------------------------------------------------------------
% Define the BART directory
%--------------------------------------------------------------------------
bart_path = json.bart_path;

%--------------------------------------------------------------------------
% Reconstruction parameters
%--------------------------------------------------------------------------
Lmax          = json.recon_parameters.Lmax;          % maximum rank of the SVD approximation of a higher-order encoding matrix
lambda        = json.recon_parameters.lambda;        % l2 regularization parameter
tol           = json.recon_parameters.tol;           % PCG tolerance
maxiter       = json.recon_parameters.maxiter;       % PCG maximum iteration 
slice_type    = json.recon_parameters.slice_type;    % type of an excitation slice: "curved" vs "flat"
cal_size      = json.recon_parameters.cal_size.';    % size of calibration region
gridding_flag = json.recon_parameters.gridding_flag; % 1=yes, 0=no
phc_flag      = json.recon_parameters.phc_flag;      % 1=yes, 0=no
cfc_flag      = json.recon_parameters.cfc_flag;      % 1=yes, 0=no
sfc_flag      = json.recon_parameters.sfc_flag;      % 1=yes, 0=no
gnc_flag      = json.recon_parameters.gnc_flag;      % 1=yes, 0=no
topup_flag    = json.recon_parameters.topup_flag;    % 1=yes, 0=no

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
% Encoding Space (Nkx, Nky, Nkz)
%--------------------------------------------------------------------------
encoded_fov(1) = header.encoding.encodedSpace.fieldOfView_mm.x * 1e-3; % [m] RO
encoded_fov(2) = header.encoding.encodedSpace.fieldOfView_mm.y * 1e-3; % [m] PE
encoded_fov(3) = header.encoding.encodedSpace.fieldOfView_mm.z * 1e-3; % [m] SL

Nkx = header.encoding.encodedSpace.matrixSize.x; % number of readout samples in k-space
Nky = header.encoding.encodedSpace.matrixSize.y; % number of phase-encoding steps in k-space
Nkz = header.encoding.encodedSpace.matrixSize.z; % number of slice-encoding steps in k-space

encoded_resolution = encoded_fov ./ [Nkx Nky Nkz]; % [m]

%--------------------------------------------------------------------------
% Recon Space (Nx, Ny, Nz)
%--------------------------------------------------------------------------
Nx = header.encoding.reconSpace.matrixSize.x; % number of samples in image space (RO)
Ny = header.encoding.reconSpace.matrixSize.y; % number of samples in image space (PE)
Nz = header.encoding.reconSpace.matrixSize.z; % number of samples in image space (SL)

%--------------------------------------------------------------------------
% Number of receive channels
%--------------------------------------------------------------------------
Nc = header.acquisitionSystemInformation.receiverChannels;

%% Calculate the number of total slices
nr_recons = nr_slices;

%% Perform image reconstruction per slice
for idx = 1:nr_recons

    %% Get the slice number
    slice_number = ind2sub(nr_slices, idx);

%     if ~(slice_number == 15 || slice_number == 16)
%         continue;
%     end

    %% Read a .cfl file
    %----------------------------------------------------------------------
    % L (1 x 1)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('L_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    L = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % ksp_cal_cartesian (Nkx x Nky x Nkz x Nc)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('ksp_cal_cartesian_slc%d_gridding%d_phc%d', slice_number, gridding_flag, phc_flag));
    if ~exist(strcat(cfl_file, '.cfl'), 'file')
        cfl_file = fullfile(output_path, sprintf('ksp_img_cartesian_slc%d_rep1_gridding%d_phc%d', slice_number, gridding_flag, phc_flag));
    end
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    ksp = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % mask_img (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('mask_img_slc%d', slice_number));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    img_mask = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % mask_cal (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('mask_cal_slc%d', slice_number));
    if exist(strcat(cfl_file, '.cfl'), 'file') % R > 1
        tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
        cal_mask = readcfl(cfl_file);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    else % either R=1 or PF
        cal_mask = img_mask .* flip(img_mask,2);
    end

    %----------------------------------------------------------------------
    % circle_mask (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('circle_mask_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    circle_mask = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

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
    % dx (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('dx_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    dx = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % dy (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('dy_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    dy = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % dz (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('dz_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    dz = readcfl(cfl_file);
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

    %----------------------------------------------------------------------
    % du (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('du_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    du = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % dv (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('dv_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    dv = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % dw (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('dw_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    dw = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % read_sign
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, 'read_sign');
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    read_sign = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % phase_sign
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, 'phase_sign');
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    phase_sign = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % initial_phase (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('initial_phase_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    initial_phase = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % U (NkNs x Lmax)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('U_cal_slc%d_%s_cfc%d_sfc%d', slice_number, slice_type, cfc_flag, sfc_flag));
    if ~exist(strcat(cfl_file, '.cfl'), 'file')
        cfl_file = fullfile(output_path, sprintf('U_img_slc%d_%s_cfc%d_sfc%d', slice_number, slice_type, cfc_flag, sfc_flag));
    end
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    U = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % S (Lmax x 1)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('S_cal_slc%d_%s_cfc%d_sfc%d', slice_number, slice_type, cfc_flag, sfc_flag));
    if ~exist(strcat(cfl_file, '.cfl'), 'file')
        cfl_file = fullfile(output_path, sprintf('S_img_slc%d_%s_cfc%d_sfc%d', slice_number, slice_type, cfc_flag, sfc_flag));
    end
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    S = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % V (Nkx x Nky x Nkz x Lmax)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('V_cal_slc%d_%s_cfc%d_sfc%d', slice_number, slice_type, cfc_flag, sfc_flag));
    if ~exist(strcat(cfl_file, '.cfl'), 'file')
        cfl_file = fullfile(output_path, sprintf('V_img_slc%d_%s_cfc%d_sfc%d', slice_number, slice_type, cfc_flag, sfc_flag));
    end
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    V = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % displacement (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    if topup_flag
        cfl_file = fullfile(topup_path, sprintf('displacement_slc%d_%s', slice_number, slice_type));
        tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
        idx1_range = (-floor(Nx/2):ceil(Nx/2)-1).' + floor(Nkx/2) + 1;
        idx2_range = (-floor(Ny/2):ceil(Ny/2)-1).' + floor(Nky/2) + 1;
        displacement = zeros(Nkx, Nky, 'single');
        displacement(idx1_range,idx2_range) = readcfl(cfl_file);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    else
        displacement = zeros(size(dv), 'single');
    end

    if read_sign < 0
        displacement = flip(displacement,1);
    end

    if phase_sign < 0
        displacement = flip(displacement,2);
    end

    %% Subselect U
    if sum(img_mask(:)) == size(U,1) % without calibration
        img_list = find(img_mask);
        cal_list = find(cal_mask);
        U_cal = U(ismember(img_list, cal_list),:);
    else
        U_cal = U;
    end

    %% Calculate calibration only k-space data
    ksp_cal = bsxfun(@times, cal_mask, ksp);

    %% Calculate parameters for Type-1 and Type-2 NUFFTs
    p1 = u;
    p2 = v;

    if gnc_flag
        p1 = p1 + du;
        p2 = p2 + dv;
    end

    if topup_flag
        p2 = p2 + displacement;
    end

    p1 = p1 / encoded_fov(1) * (2 * pi); % RO [-pi,pi]
    p2 = p2 / encoded_fov(2) * (2 * pi); % PE [-pi,pi]

    %% Calculate a support mask
    support_mask = zeros(Nkx, Nky, Nkz, 'single');
    support_mask((abs(p1) < pi) & (abs(p2) < pi) & (circle_mask > 0)) = 1;

    %% Select spatial positions within a support mask
    p1 = p1(support_mask > 0);
    p2 = p2(support_mask > 0);

    %% Set parameters for Type-1 and Type-2 NUFFTs
    eps = 1e-6;
    iflag = -1;

    %% Calculate "fake" coil sensitivity maps (Nkx x Nky x Nkz)
    sens = complex(ones(Nkx, Nky, Nkz, 'single'));

    %% Type-1 NUFFT based Cartesian MaxGIRF operators
    Ah = @(x) cartesian_maxgirf_2d_adjoint(x, sens, cal_mask, p1, p2, iflag, eps, support_mask, U_cal, V, L, initial_phase);
    AhA = @(x) cartesian_maxgirf_2d_normal(x, sens, cal_mask, p1, p2, iflag, eps, support_mask, U_cal, V, L, initial_phase, lambda);

    %% Perform Cartesian MaxGIRF reconstruction
    cimg = complex(zeros(Nkx, Nky, Nkz, Nc, 'single'));
    for c = 1:Nc
        clear cartesian_maxgirf_2d_normal;
        tstart = tic; fprintf('%s:(c=%2d/%2d) Performing Cartesian MaxGIRF reconstruction:\n', datetime, c, Nc);
        %[img, flag, relres, iter, resvec] = pcg(@(x) AhA(x), Ah(ksp_cal(:,:,:,c)), tol, maxiter);
        %cimg(:,:,:,c) = reshape(img, [Nkx Nky Nkz]);
        cimg(:,:,:,c) = reshape(Ah(ksp_cal(:,:,:,c)), [Nkx Nky Nkz]);
        fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));
    end

    %% Write a .cfl file
    %----------------------------------------------------------------------
    % cimg (Nkx x Nky x Nkz x Nc)
    %----------------------------------------------------------------------
    cimg_filename = sprintf('cimg_cal_slc%d_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_topup%d_i%d_l%4.2f', slice_number, slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, topup_flag, maxiter, lambda);
    cfl_file = fullfile(output_path, cimg_filename);
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, cimg);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Calculate gridded k-space
    %----------------------------------------------------------------------
    % Siemens: k-space <=> image space
    % BART:                image space <=> k-space
    %----------------------------------------------------------------------
    tstart = tic; fprintf('%s: Applying forward FFT to move from image space to k-space... ', datetime);
    kgrid = cimg;
    for dim = 1:3
        kgrid = 1 / sqrt(size(kgrid,dim)) * fftshift(fft(ifftshift(kgrid, dim), [], dim), dim);
    end
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Write a .cfl file
    %----------------------------------------------------------------------
    % kgrid (Nkx x Nky x Nkz x Nc)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('kgrid_slc%d_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_topup%d', slice_number, slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, topup_flag));
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, kgrid);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Estimate coil sensitivity using ESPIRiT calibration
    %----------------------------------------------------------------------
    % ESPIRiT calibration
    %----------------------------------------------------------------------
    % Usage: ecalib [-t f] [-c f] [-k d:d:d] [-r d:d:d] [-m d] [-S] [-W] [-I] [-1]
    %               [-P] [-v f] [-a] [-d d] <kspace> <sensitivities> [<ev-maps>]
    %
    % Estimate coil sensitivities using ESPIRiT calibration.
    % Optionally outputs the eigenvalue maps.
    %
    % -t threshold     This determines the size of the null-space
    % -c crop_value    Crop the sensitivities if the eigenvalue is smaller than {crop_value}
    % -k ksize         kernel size
    % -r cal_size      Limits the size of the calibration region
    % -m maps          Number of maps to compute
    % -S               create maps with smooth transitions (Soft-SENSE)
    % -W               soft-weighting of the singular vectors
    % -I               intensity correction
    % -1               perform only first part of the calibration
    % -P               Do not rotate the phase with respect to the first principal component
    % -v variance      Variance of noise in data
    % -a               Automatically pick thresholds
    % -d level         Debug level
    % -h               help
    %----------------------------------------------------------------------
    kgrid_file   = strcat(bart_output_path, sprintf('kgrid_slc%d_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_topup%d', slice_number, slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, topup_flag));
    sens_file    = strcat(bart_output_path, sprintf('sens_slc%d_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_topup%d', slice_number, slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, topup_flag));
    ev_maps_file = strcat(bart_output_path, sprintf('ev_maps_slc%d_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_topup%d', slice_number, slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, topup_flag));
    command = sprintf('%s ecalib -t 0.001 -c 0 -k6:6:6 -r%d:%d:%d -m 1 -d5 %s %s %s', bart_command, cal_size(1), cal_size(2), cal_size(3), kgrid_file, sens_file, ev_maps_file);
    tstart = tic; fprintf('%s:[BART] Estimating coil sensitivities using ESPIRiT calibration:\n%s\n', datetime, command);
    [status_ecalib,result_ecalib] = system(command);
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    %% Read a .cfl file
    %----------------------------------------------------------------------
    % sens (Nkx x Nky x Nkz x Nc)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('sens_slc%d_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_topup%d', slice_number, slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, topup_flag));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    sens = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Read a .cfl file
    %----------------------------------------------------------------------
    % ev_maps (Nkx x Nky x Nkz)
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('ev_maps_slc%d_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_topup%d', slice_number, slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, topup_flag));
    tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
    ev_maps = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %% Display the magnitude of coil images
    nr_rows = 2;
    nr_cols = 7;
    cimg_montage = complex(zeros(Nkx * nr_rows, Nky * nr_cols, 'single'));
    sens_montage = complex(zeros(Nkx * nr_rows, Nky * nr_cols, 'single'));

    count = 0;
    for idx1 = 1:nr_rows
        for idx2 = 1:nr_cols
            idx1_range = (1:Nkx).' + (idx1 - 1) * Nkx;
            idx2_range = (1:Nky).' + (idx2 - 1) * Nky;
            count = count + 1;
            cimg_montage(idx1_range,idx2_range) = cimg(:,:,1,count);
            sens_montage(idx1_range,idx2_range) = sens(:,:,1,count);
            if count >= Nc
                break;
            end
        end
    end

    title_text1 = sprintf('Cartesian MaxGIRF, SLC = %d, %s slice', slice_number, slice_type);
    title_text2 = sprintf('Gridding/PHC/CFC/SFC/GNC/TOPUP = %d/%d/%d/%d/%d/%d', gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, topup_flag);
    fig_filename = sprintf('cimg_cal_slc%d_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_i%d_l%4.2f', slice_number, slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, maxiter, lambda);

    figure('Color', 'k', 'Position', [1 5 1239 973]);
    imagesc(abs(cimg_montage));
    axis image off;
    colormap(gray(256));
    caxis([0 5]);
    title({'Magnitude of coil images (calibration data)', ...
        title_text1, title_text2}, 'Color', 'w', 'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 16);
    export_fig(fullfile(output_path, sprintf('%s_mag', fig_filename)), '-r300', '-tif', '-c[360, 140, 700, 440]'); % [top,right,bottom,left]
    close gcf;

    %% Display the phase of coil images
    figure('Color', 'k', 'Position', [1 5 1239 973]);
    imagesc(angle(cimg_montage) * 180 / pi);
    axis image off;
    caxis([-180 180]);
    colormap(hsv(256));
    hc = colorbar;
    set(hc, 'Color', 'w', 'FontSize', 14, 'Position', [0.9152 0.2857 0.0121 0.4470], 'TickLabelInterpreter', 'latex');
    title(hc, '[deg]', 'Color', 'w', 'Interpreter', 'latex');
    title({'Phase of coil images (calibration data)', title_text1, title_text2}, 'Color', 'w', 'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 16);
    export_fig(fullfile(output_path, sprintf('%s_phase', fig_filename)), '-r300', '-tif', '-c[360, 140, 700, 440]'); % [top,right,bottom,left]
    close gcf;

    %% Display the magnitude of CSMs
    fig_filename = sprintf('sens_cal_slc%d_%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_i%d_l%4.2f', slice_number, slice_type, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, maxiter, lambda);

    figure('Color', 'k', 'Position', [1 5 1239 973]);
    imagesc(abs(sens_montage));
    axis image off;
    colormap(gray(256));
    caxis([0 1]);
    title({'Magnitude of coil sensitivity maps', title_text1, title_text2}, 'Color', 'w', 'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 16);
    export_fig(fullfile(output_path, sprintf('%s_mag', fig_filename)), '-r300', '-tif', '-c[360, 140, 700, 440]'); % [top,right,bottom,left]
    close gcf;

    %% Display the phase of CSMs
    figure('Color', 'k', 'Position', [1 5 1239 973]);
    imagesc(angle(sens_montage) * 180 / pi);
    axis image off;
    caxis([-180 180]);
    colormap(hsv(256));
    hc = colorbar;
    set(hc, 'Color', 'w', 'FontSize', 14, 'Position', [0.9152 0.2857 0.0121 0.4470], 'TickLabelInterpreter', 'latex');
    title(hc, '[deg]', 'Color', 'w', 'Interpreter', 'latex');
    title({'Phase of coil sensitivity maps', title_text1, title_text2}, 'Color', 'w', 'Interpreter', 'latex', 'FontWeight', 'normal', 'FontSize', 16);
    export_fig(fullfile(output_path, sprintf('%s_phase', fig_filename)), '-r300', '-tif', '-c[360, 140, 700, 440]');
    close gcf;
end
