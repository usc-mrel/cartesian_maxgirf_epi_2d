% demo_fsl_topup.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 04/30/2025, Last modified: 04/30/2025

%% Clean slate
close all; clear all; clc;

%% Start a stopwatch timer
start_time = tic;

%% Set source directories
if ispc
    package_path = 'D:\cartesian_maxgirf_epi_2d';
else
    package_path = '/server/sdata/nlee/cartesian_maxgirf_epi_2d';
end

%% Add source directories to search path
addpath(genpath(package_path));

%% Define the full path of FSL
fsl_path = '/home/image/fsl';

%% Define the full path of dcm2niix
dcm2niix_path = 'D:\dcm2niix_win';

%% Define the full path of data directory
%--------------------------------------------------------------------------
% Figure 2
%--------------------------------------------------------------------------
data1_path = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\dicom\S3_ep2d_se_bw1002_cor_RL';
data2_path = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\dicom\S5_ep2d_se_bw1002_cor_LR';
output_path = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\topup_RL_LR';

%--------------------------------------------------------------------------
% Figure 3
%--------------------------------------------------------------------------
data1_path = 'D:\cartesian_maxgirf_epi_2d\data\acr_phantom_20250209\dicom\S24_ep2d_se_bw976_tra_AP_avg16_fov256_256_pf_R2_S6_DIS2D';
data2_path = 'D:\cartesian_maxgirf_epi_2d\data\acr_phantom_20250209\dicom\S34_ep2d_se_bw976_tra_PA_avg16_fov256_256_pf_R2_S16_DIS2D';
output_path = 'D:\cartesian_maxgirf_epi_2d\data\acr_phantom_20250209\topup_AP_PA';

%% Make an output directory
mkdir(output_path);

%% Set up FSL commands
%--------------------------------------------------------------------------
% Define an FSL command
%--------------------------------------------------------------------------
if ispc
    command_prefix = 'wsl';
else
    command_prefix = '';
end
fsl_command = sprintf('%s . %s/etc/fslconf/fsl.sh; %s/bin', command_prefix, fsl_path, fsl_path);

%--------------------------------------------------------------------------
% Translate from a Windows path to a WSL path
%--------------------------------------------------------------------------
if ispc
    fsl_output_path = strrep(output_path, '\', '/');
    fsl_output_path = sprintf('/mnt/%s/%s/', lower(fsl_output_path(1)), fsl_output_path(4:end));
else
    fsl_output_path = sprintf('%s/', output_path);
end

%% Convert DICOM images to a NIfTI file (magnitude)
command = sprintf('%s%sdcm2niix -f img1 -o %s -z y %s', dcm2niix_path, filesep, output_path, data1_path);
tstart = tic; fprintf('%s:[dcm2niix] Converting DICOM images to a NIfTI file:\n%s\n', datetime, command);
[status_dcm2niix1,result_dcm2niix1] = system(command);
fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

%% Convert DICOM images to a NIfTI file (phase)
command = sprintf('%s%sdcm2niix -f img2 -o %s -z y %s', dcm2niix_path, filesep, output_path, data2_path);
tstart = tic; fprintf('%s:[dcm2niix] Converting DICOM images to a NIfTI file:\n%s\n', datetime, command);
[status_dcm2niix2,result_dcm2niix2] = system(command);
fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

%% Concatenate images in time
output_file = strcat(fsl_output_path, 'img');
img1_file   = strcat(fsl_output_path, 'img1');
img2_file   = strcat(fsl_output_path, 'img2');
command = sprintf('%s/fslmerge -t %s %s %s', fsl_command, output_file, img1_file, img2_file);
tstart = tic; fprintf('%s:[FSL] Concatenating images in time using FSL fslmerge:\n%s\n', datetime, command);
[status_fslmerge,result_fslmerge] = system(command);
fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

%% Create "acq_param.txt"
%--------------------------------------------------------------------------
% Open a. txt file
%--------------------------------------------------------------------------
txt_file = sprintf('%s%sacq_param.txt', output_path, filesep);
fid_txt = fopen(txt_file, 'w');

for idx = 1:2
    %----------------------------------------------------------------------
    % Define the full path of a .json file
    %----------------------------------------------------------------------
    json_file = sprintf('%s%simg%d.json', output_path, filesep, idx);

    %----------------------------------------------------------------------
    % Read a .json file
    %----------------------------------------------------------------------
    tstart = tic; fprintf('%s: Reading a .json file: %s... ', datetime, json_file);
    fid_json = fopen(json_file);
    json_txt = fread(fid_json, [1 inf], 'char=>char');
    fclose(fid_json);
    json = jsondecode(json_txt);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Write a row in a .txt file
    %----------------------------------------------------------------------
    if length(json.PhaseEncodingDirection) == 1
        sign = 1;
    else
        sign = -1;
    end

    if strfind(json.PhaseEncodingDirection, 'i')
        fprintf(fid_txt, '%d %d %d %8.6f\n', sign, 0, 0, json.TotalReadoutTime);
    else
        fprintf(fid_txt, '%d %d %d %8.6f\n', 0, sign, 0, json.TotalReadoutTime);
    end
end
fclose(fid_txt);

%% Calculate a fieldmap using topup
%--------------------------------------------------------------------------
% Part of FSL (ID: "")
% topup
%
% Usage:
% topup --imain=<some 4D image> --datain=<text file> --config=<text file with parameters> --out=my_topup_results
%
%
% Compulsory arguments (You MUST set one or more of):
%         --imain         name of 4D file with images
%
% Optional arguments (You may optionally specify one or more of):
%         --datain        name of text file with PE directions/times
%         --acqp          alternative way to specify text file with PE directions/times
%         --out           base-name of output files (spline coefficients (Hz) and movement parameters)
%         --fout          name of image file with field (Hz)
%         --iout          name of 4D image file with unwarped images
%         --featout       base-name of output for export to FEAT
%         --logout        Name of log-file
%         --warpres       (approximate) resolution (in mm) of warp basis for the different sub-sampling levels, default 10
%         --subsamp       sub-sampling scheme, default 1
%         --fwhm          FWHM (in mm) of gaussian smoothing kernel, default 8
%         --config        Name of config file specifying command line arguments
%         --miter         Max # of non-linear iterations, default 5
%         --lambda        Weight of regularisation, default depending on --ssqlambda and --regmod switches. See user documetation.
%         --ssqlambda     If set (=1), lambda is weighted by current ssq, default 1
%         --regmod        Model for regularisation of warp-field [membrane_energy bending_energy], default bending_energy
%         --estmov        Estimate movements if set, default 1 (true)
%         --minmet        Minimisation method 0=Levenberg-Marquardt, 1=Scaled Conjugate Gradient, default 0 (LM)
%         --splineorder   Order of spline, 2->Qadratic spline, 3->Cubic spline. Default=3
%         --numprec       Precision for representing Hessian, double or float. Default double
%         --interp        Image interpolation model, linear or spline. Default spline
%         --scale         If set (=1), the images are individually scaled to a common mean, default 0 (false)
%         --regrid        If set (=1), the calculations are done in a different grid, default 1 (true)
%         --nthr          Number of threads to use (cannot be greater than numbers of hardware cores), default 1
%         -h,--help       display help info
%         -v,--verbose    Print diagonostic information while running
%         -h,--help       display help info
%--------------------------------------------------------------------------
imain_file  = strcat(fsl_output_path, 'img');
datain_file = strcat(fsl_output_path, 'acq_param.txt');
%\\wsl.localhost\Ubuntu-20.04\home\image\fsl\pkgs\fsl-topup-2203.5-hb6de94e_0\etc\flirtsch
config_file = strcat(fsl_path, '/pkgs/fsl-topup-2203.5-hb6de94e_0/etc/flirtsch/b02b0.cnf');
out_file    = strcat(fsl_output_path, 'topup');
fout_file   = strcat(fsl_output_path, 'topup_fieldmap');
iout_file   = strcat(fsl_output_path, 'topup_img_cor');
command = sprintf('%s/topup --imain=%s --datain=%s --config=%s --out=%s --fout=%s --iout=%s --verbose', fsl_command, imain_file, datain_file, config_file, out_file, fout_file, iout_file);
tstart = tic; fprintf('%s:[FSL] Calculating a fieldmap using FSL topup:\n%s\n', datetime, command);
[status_topup,result_topup] = system(command);
fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

%% Save the output of the topup command as a .txt file
[filepath,output_filename,ext] = fileparts(output_path);
if 1
    txt_file = fullfile(output_path, 'topup_debug.txt');
    [fid,message] = fopen(txt_file, 'w');
    fwrite(fid, result_topup, 'char');
    fclose(fid);
end

%% Calculate a distortion-corrected image using applytopup
%--------------------------------------------------------------------------
% Part of FSL (ID: "")
% applytopup (Version 1.0)
% Copyright(c) 2009, University of Oxford (Jesper Andersson)
%
% Usage:
% applytopup -in=topdn,botup --topup=mytu --inindex=1,2 --out=hifi
% applytopup -in=topdn --topup=mytu --inindex=1 --method=jac --interp=spline --out=hifi
%
%
% Compulsory arguments (You MUST set one or more of):
%         -i,--imain      comma separated list of names of input image (to be corrected)
%         -a,--datain     name of text file with PE directions/times
%         -x,--inindex    comma separated list of indicies into --datain of the input image (to be corrected)
%         -t,--topup      name of field/movements (from topup)
%         -o,--out        basename for output (warped) image
%
% Optional arguments (You may optionally specify one or more of):
%         -m,--method     Use jacobian modulation (jac) or least-squares resampling (lsr), default=lsr.
%         -n,--interp     interpolation method {trilinear,spline}, default=spline
%         -d,--datatype   Force output data type [preserve char short int float double].
%         -v,--verbose    switch on diagnostic messages
%         -h,--help       display this message
%--------------------------------------------------------------------------
for idx = 1:2
    imain_file  = strcat(fsl_output_path, sprintf('img%d', idx));
    datain_file = strcat(fsl_output_path, 'acq_param.txt');
    topup_file  = strcat(fsl_output_path, 'topup');
    out_file    = strcat(fsl_output_path, sprintf('applytopup_img%d_cor', idx));
    command = sprintf('%s/applytopup --imain=%s --inindex=%d --datain=%s --topup=%s --method=jac --out=%s --verbose', fsl_command, imain_file, idx, datain_file, topup_file, out_file);
    tstart = tic; fprintf('%s:[FSL] Calculating a distortion-corrected image using FSL applytopup:\n%s\n', datetime, command);
    [status_applytopup,result_applytopup] = system(command);
    fprintf('%s: done! (%6.4f/%6.4f sec)\n', datetime, toc(tstart), toc(start_time));

    %% Save the output of the topup command as a .txt file
    txt_file = fullfile(output_path, sprintf('applytopup_img%d_debug.txt', idx));
    [fid,message] = fopen(txt_file, 'w');
    fwrite(fid, result_applytopup, 'char');
    fclose(fid);
end

%% Define the full path of a filename
nii_img_file                 = fullfile(output_path, 'img.nii.gz');
nii_topup_img_cor_file       = fullfile(output_path, 'topup_img_cor.nii.gz');
nii_topup_fieldmap_file      = fullfile(output_path, 'topup_fieldmap.nii.gz');
nii_applytopup_img1_cor_file = fullfile(output_path, 'applytopup_img1_cor.nii.gz');
nii_applytopup_img2_cor_file = fullfile(output_path, 'applytopup_img2_cor.nii.gz');
acq_param_file               = fullfile(output_path, 'acq_param.txt');

%% Get the header information
nii_info = niftiinfo(nii_img_file);

%% Read a .nii file
%--------------------------------------------------------------------------
% img
%--------------------------------------------------------------------------
img = single(niftiread(nii_img_file));

%--------------------------------------------------------------------------
% topup_img_cor
%--------------------------------------------------------------------------
topup_img_cor = single(niftiread(nii_topup_img_cor_file));

%--------------------------------------------------------------------------
% topup_fieldmap
%--------------------------------------------------------------------------
topup_fieldmap = niftiread(nii_topup_fieldmap_file); % single

%--------------------------------------------------------------------------
% applytopup_img1_cor
%--------------------------------------------------------------------------
applytopup_img1_cor = single(niftiread(nii_applytopup_img1_cor_file));

%--------------------------------------------------------------------------
% applytopup_img2_cor
%--------------------------------------------------------------------------
applytopup_img2_cor = single(niftiread(nii_applytopup_img2_cor_file));

%% Read a .txt file
fid = fopen(acq_param_file, 'r');
acq_param = fscanf(fid, '%f %f %f %f', [4 inf]).';
fclose(fid);

%% Set the total readout time [sec]
total_readout_time = acq_param(1,4); % [sec]

%% Calculate the recon fov [m]
recon_fov = nii_info.ImageSize(1:3) .* nii_info.PixelDimensions(1:3) * 1e-3; % [m]

%% Calculate the effective echo spacing [sec]
effective_echo_spacing = total_readout_time / nii_info.ImageSize(1);

%% Calculate the bandwidth in the phase encoding direction [Hz/mm]
bandwidth_pe = 1 / (recon_fov(1) * 1e3 * effective_echo_spacing);

%% Calculate a mask
mask = (img(:,:,:,1) > max(img(:)) * 0.1);
mask = bwareaopen(mask, 60); % Keep only blobs with an area of 60 pixels or more.
se = strel('disk',5);
mask_dilated = imdilate(mask,se);

%% Apply a mask
topup_fieldmap_masked = topup_fieldmap .* mask_dilated;

%% Set imaging parameters
[Nx,Ny,nr_slices] = size(topup_fieldmap);

%% Calculate a reorient function (from NIfTI to Siemens RCS)
if strcmp(json.ImageOrientationText, 'Cor')
    reorient = @(x) transpose(flip(x,2)); % tested on Cor RL, Cor LR
elseif strcmp(json.ImageOrientationText, 'Tra')
    reorient = @(x) flip(x,2); % axial
end

%% Write as a .cfl file
for slice_number = 1:nr_slices
    %------------------------------------------------------------------
    % Calculate the actual slice number for Siemens interleaved multislice imaging
    %------------------------------------------------------------------
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

    fieldmap = reorient(topup_fieldmap(:,:,actual_slice_number));

    %----------------------------------------------------------------------
    % Write as a .cfl file
    %----------------------------------------------------------------------
    slice_type = 'flat';
    cfl_file = fullfile(output_path, sprintf('fieldmap_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, fieldmap);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Write as a .cfl file
    %----------------------------------------------------------------------
    cfl_file = fullfile(output_path, sprintf('displacement_slc%d_%s', slice_number, slice_type));
    tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
    writecfl(cfl_file, fieldmap / bandwidth_pe * 1e-3); % [Hz] / [Hz/mm] * [m/1e3mm] => [m]
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Display images
nr_rows = 4;
nr_cols = 5;

img1_montage                = complex(zeros(Nx * nr_rows, Ny * nr_cols, 'single'));
img2_montage                = complex(zeros(Nx * nr_rows, Ny * nr_cols, 'single'));
topup_fieldmap_montage      = complex(zeros(Nx * nr_rows, Ny * nr_cols, 'single'));
topup_img1_cor_montage      = complex(zeros(Nx * nr_rows, Ny * nr_cols, 'single'));
topup_img2_cor_montage      = complex(zeros(Nx * nr_rows, Ny * nr_cols, 'single'));
applytopup_img1_cor_montage = complex(zeros(Nx * nr_rows, Ny * nr_cols, 'single'));
applytopup_img2_cor_montage = complex(zeros(Nx * nr_rows, Ny * nr_cols, 'single'));

count = 0;
for idx1 = 1:nr_rows
    for idx2 = 1:nr_cols
        count = count + 1;
        idx1_range = (1:Nx).' + (idx1 - 1) * Nx;
        idx2_range = (1:Ny).' + (idx2 - 1) * Ny;
        img1_montage(idx1_range,idx2_range) = reorient(img(:,:,count,1));
        img2_montage(idx1_range,idx2_range) = reorient(img(:,:,count,2));
        topup_fieldmap_montage(idx1_range,idx2_range) = reorient(topup_fieldmap(:,:,count));
        topup_img1_cor_montage(idx1_range,idx2_range) = reorient(topup_img_cor(:,:,count,1));
        topup_img2_cor_montage(idx1_range,idx2_range) = reorient(topup_img_cor(:,:,count,2));
        applytopup_img1_cor_montage(idx1_range,idx2_range) = reorient(applytopup_img1_cor(:,:,count));
        applytopup_img2_cor_montage(idx1_range,idx2_range) = reorient(applytopup_img2_cor(:,:,count));
        if count >= nr_slices
            break;
        end
    end
    if count >= nr_slices
        break;
    end
end

%% Display a fieldmap [Hz]
figure('Color', 'w', 'Position', [346 11 1296 939]);
imagesc(topup_fieldmap_montage);
axis image off;
caxis([-20 20]);
colormap(colorcet('D1'));
hc = colorbar;
climits = hc.Limits;
set(hc, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
title(hc, '[Hz]', 'FontSize', 14, 'Interpreter', 'latex');
title('Fieldmap estimated by topup', 'FontSize', 14, 'Interpreter', 'latex');
export_fig(fullfile(output_path, 'fieldmap'), '-r300', '-tif', '-c[100, 100, 300, 320]'); % [top,right,bottom,left]
close gcf;

%% Display a displacement field [mm]
figure('Color', 'w', 'Position', [346 11 1296 939]);
imagesc(topup_fieldmap_montage / bandwidth_pe); % [Hz] / [Hz/mm] => [mm]
axis image off;
caxis([-4 4]);
colormap(flip(brewermap([],"RdBu"),1));
hc = colorbar;
set(hc, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
title(hc, '[mm]', 'FontSize', 14, 'Interpreter', 'latex');
title('Displacement estimated by topup', 'FontSize', 14, 'FontWeight', 'normal', 'Interpreter', 'latex');
export_fig(fullfile(output_path, 'displacement'), '-r300', '-tif', '-c[100, 100, 300, 320]');
close gcf;

%% Display uncorrected images (img 1)
figure('Color', 'w', 'Position', [346 11 1296 939]);
imagesc(img1_montage);
axis image off;
climits = get(gca, 'clim');
colormap(gray(256));
hc = colorbar;
set(hc, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
title('Uncorrected images (img 1)', 'FontSize', 14, 'FontWeight', 'normal', 'Interpreter', 'latex');
export_fig(fullfile(output_path, 'img1'), '-r300', '-tif', '-c[100, 100, 300, 320]');
close gcf;

%% Display uncorrected images (img 2)
figure('Color', 'w', 'Position', [346 11 1296 939]);
imagesc(img2_montage);
axis image off;
caxis(climits);
colormap(gray(256));
hc = colorbar;
set(hc, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
title('Uncorrected images (img 2)', 'FontSize', 14, 'FontWeight', 'normal', 'Interpreter', 'latex');
export_fig(fullfile(output_path, 'img2'), '-r300', '-tif', '-c[100, 100, 300, 320]');
close gcf;

%% Display corrected images by topup (img 1)
figure('Color', 'w', 'Position', [346 11 1296 939]);
imagesc(topup_img1_cor_montage);
axis image off;
caxis(climits);
colormap(gray(256));
hc = colorbar;
set(hc, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
title('Corrected images by topup (img 1)', 'FontSize', 14, 'FontWeight', 'normal', 'Interpreter', 'latex');
export_fig(fullfile(output_path, 'topup_img1_cor'), '-r300', '-tif', '-c[100, 100, 300, 320]');
close gcf;

%% Display corrected images by topup (img 2)
figure('Color', 'w', 'Position', [346 11 1296 939]);
imagesc(topup_img2_cor_montage);
axis image off;
caxis(climits);
colormap(gray(256));
hc = colorbar;
set(hc, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
title('Corrected images by topup (img 2)', 'FontSize', 14, 'FontWeight', 'normal', 'Interpreter', 'latex');
export_fig(fullfile(output_path, 'topup_img2_cor'), '-r300', '-tif', '-c[100, 100, 300, 320]');
close gcf;

%% Display corrected images by applytopup (img 1)
figure('Color', 'w', 'Position', [346 11 1296 939]);
imagesc(applytopup_img1_cor_montage);
axis image off;
caxis(climits);
colormap(gray(256));
hc = colorbar;
set(hc, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
title('Corrected images by applytopup (img 1)', 'FontSize', 14, 'FontWeight', 'normal', 'Interpreter', 'latex');
export_fig(fullfile(output_path, 'applytopup_img1_cor'), '-r300', '-tif', '-c[100, 100, 300, 320]');
close gcf;

%% Display corrected images by topup (img 2)
figure('Color', 'w', 'Position', [346 11 1296 939]);
imagesc(applytopup_img2_cor_montage);
axis image off;
caxis(climits);
colormap(gray(256));
hc = colorbar;
set(hc, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
title('Corrected images by applytopup (img 2)', 'FontSize', 14, 'FontWeight', 'normal', 'Interpreter', 'latex');
export_fig(fullfile(output_path, 'applytopup_img2_cor'), '-r300', '-tif', '-c[100, 100, 300, 320]');
close gcf;
