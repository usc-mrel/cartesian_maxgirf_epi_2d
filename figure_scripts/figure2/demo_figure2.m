% demo_figure2.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 08/20/2024, Last modified: 08/20/2024

%% Clean slate
close all; clear all; clc;

%% Start a stopwatch timer
start_time = tic;

%% Define the full path of an output path
% Phase encoding direction: R >> L
output_path1 = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\meas_MID00057_FID21486_ep2d_se_bw1002_cor_RL_gridding0_phc0_cfc0_sfc0_gnc0_topup0';
output_path2 = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\meas_MID00057_FID21486_ep2d_se_bw1002_cor_RL_gridding1_phc0_cfc0_sfc0_gnc0_topup0';
output_path3 = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\meas_MID00057_FID21486_ep2d_se_bw1002_cor_RL_gridding1_phc1_cfc0_sfc0_gnc0_topup0';
output_path4 = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\meas_MID00057_FID21486_ep2d_se_bw1002_cor_RL_gridding1_phc1_cfc1_sfc0_gnc0_topup0';
output_path5 = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\meas_MID00057_FID21486_ep2d_se_bw1002_cor_RL_gridding1_phc1_cfc1_sfc0_gnc1_topup0';
output_path6 = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\meas_MID00057_FID21486_ep2d_se_bw1002_cor_RL_gridding1_phc1_cfc1_sfc0_gnc1_topup1';
output_path7 = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\meas_MID00074_FID21503_tse_cor_FH_gnc1';

%% Set the full path of a directory
dicom_path = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\dicom\S3_ep2d_se_bw1002_cor_RL';

%% Define variables
Nkx = 256;
Nky = 128;
nr_slices = 20;
Lmax = 20;

N1 = 128;
N2 = 128;
N3 = nr_slices;

%% Read a .cfl file
img_maxgirf1 = complex(zeros(Nkx, Nky, nr_slices, 'single')); % Grid/PHC/CFC/SFC/GNC/TOPUP = 0/0/0/0/0/0
img_maxgirf2 = complex(zeros(Nkx, Nky, nr_slices, 'single')); % Grid/PHC/CFC/SFC/GNC/TOPUP = 1/0/0/0/0/0
img_maxgirf3 = complex(zeros(Nkx, Nky, nr_slices, 'single')); % Grid/PHC/CFC/SFC/GNC/TOPUP = 1/1/0/0/0/0
img_maxgirf4 = complex(zeros(Nkx, Nky, nr_slices, 'single')); % Grid/PHC/CFC/SFC/GNC/TOPUP = 1/1/1/0/0/0
img_maxgirf5 = complex(zeros(Nkx, Nky, nr_slices, 'single')); % Grid/PHC/CFC/SFC/GNC/TOPUP = 1/1/1/0/1/0
img_maxgirf6 = complex(zeros(Nkx, Nky, nr_slices, 'single')); % Grid/PHC/CFC/SFC/GNC/TOPUP = 1/1/1/0/1/1

img_maxgirf_tse = complex(zeros(Nkx, 136, nr_slices, 'single'));

x = zeros(Nkx, Nky, nr_slices, 'single');
y = zeros(Nkx, Nky, nr_slices, 'single');
z = zeros(Nkx, Nky, nr_slices, 'single');

dx = zeros(Nkx, Nky, nr_slices, 'single');
dy = zeros(Nkx, Nky, nr_slices, 'single');
dz = zeros(Nkx, Nky, nr_slices, 'single');

S = zeros(Lmax, nr_slices, 'single');

for slice_number = 1:nr_slices
    %----------------------------------------------------------------------
    % Calculate the actual slice number for Siemens interleaved multislice imaging
    % slice_number: acquisition slice number
    % actual_slice_number is used to retrieve slice information from TWIX format
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
    % Grid/PHC/CFC/SFC/GNC/TOPUP = 0/0/0/0/0/0
    %----------------------------------------------------------------------
    img_filename1 = sprintf('img_maxgirf_slc%d_rep1_flat_gridding0_phc0_cfc0_sfc0_gnc0_topup0_i6_l0.00', slice_number);
    cfl_file = fullfile(output_path1, img_filename1);
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices);
    img_maxgirf1(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Grid/PHC/CFC/SFC/GNC/TOPUP = 1/0/0/0/0/0
    %----------------------------------------------------------------------
    img_filename2 = sprintf('img_maxgirf_slc%d_rep1_flat_gridding1_phc0_cfc0_sfc0_gnc0_topup0_i6_l0.00', slice_number);
    cfl_file = fullfile(output_path2, img_filename2);
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices);
    img_maxgirf2(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Grid/PHC/CFC/SFC/GNC/TOPUP = 1/1/0/0/0/0
    %----------------------------------------------------------------------
    img_filename3 = sprintf('img_maxgirf_slc%d_rep1_flat_gridding1_phc1_cfc0_sfc0_gnc0_topup0_i6_l0.00', slice_number);
    cfl_file = fullfile(output_path3, img_filename3);
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices);
    img_maxgirf3(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Grid/PHC/CFC/SFC/GNC/TOPUP = 1/1/1/0/0/0
    %----------------------------------------------------------------------
    img_filename4 = sprintf('img_maxgirf_slc%d_rep1_flat_gridding1_phc1_cfc1_sfc0_gnc0_topup0_i6_l0.00', slice_number);
    cfl_file = fullfile(output_path4, img_filename4);
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices);
    img_maxgirf4(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Grid/PHC/CFC/SFC/GNC/TOPUP = 1/1/1/0/1/0
    %----------------------------------------------------------------------
    img_filename5 = sprintf('img_maxgirf_slc%d_rep1_flat_gridding1_phc1_cfc1_sfc0_gnc1_topup0_i6_l0.00', slice_number);
    cfl_file = fullfile(output_path5, img_filename5);
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices);
    img_maxgirf5(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % Grid/PHC/CFC/SFC/GNC/TOPUP = 1/1/1/0/1/1
    %----------------------------------------------------------------------
    img_filename6 = sprintf('img_maxgirf_slc%d_rep1_flat_gridding1_phc1_cfc1_sfc0_gnc1_topup1_i6_l0.00', slice_number);
    cfl_file = fullfile(output_path6, img_filename6);
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices);
    img_maxgirf6(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % TSE: Type-1 NUFFFT based CG-SENSE
    %----------------------------------------------------------------------
    img_filename7 = sprintf('img_type1_slc%d_flat_gnc1_i6_l0.00', slice_number);
    cfl_file = fullfile(output_path7, img_filename7);
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices);
    img_maxgirf_tse(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % x (Nkx x Nky)
    %----------------------------------------------------------------------
    img_filename8 = sprintf('x_slc%d_flat', slice_number);
    cfl_file = fullfile(output_path1, img_filename8);
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices);
    x(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % y (Nkx x Nky)
    %----------------------------------------------------------------------
    img_filename9 = sprintf('y_slc%d_flat', slice_number);
    cfl_file = fullfile(output_path1, img_filename9);
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices);
    y(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % z (Nkx x Nky)
    %----------------------------------------------------------------------
    img_filename10 = sprintf('z_slc%d_flat', slice_number);
    cfl_file = fullfile(output_path1, img_filename10);
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices);
    z(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % dx (Nkx x Nky)
    %----------------------------------------------------------------------
    img_filename11 = sprintf('dx_slc%d_flat', slice_number);
    cfl_file = fullfile(output_path1, img_filename11);
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices);
    dx(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % dy (Nkx x Nky)
    %----------------------------------------------------------------------
    img_filename12 = sprintf('dy_slc%d_flat', slice_number);
    cfl_file = fullfile(output_path1, img_filename12);
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices);
    dy(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % dz (Nkx x Nky)
    %----------------------------------------------------------------------
    img_filename13 = sprintf('dz_slc%d_flat', slice_number);
    cfl_file = fullfile(output_path1, img_filename13);
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices);
    dz(:,:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

    %----------------------------------------------------------------------
    % S (Lmax x 1)
    %----------------------------------------------------------------------
    img_filename14 = sprintf('S_img_slc%d_flat_cfc1_sfc0', slice_number);
    cfl_file = fullfile(output_path4, img_filename14);
    tstart = tic; fprintf('%s:(SLC=%2d/%2d) Reading a .cfl file: %s... ', datetime, slice_number, nr_slices);
    S(:,actual_slice_number) = readcfl(cfl_file);
    fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
end

%% Remove readout oversampling
idx1_range = (-floor(N1/2):ceil(N1/2)-1).' + floor(Nkx/2) + 1;
idx2_range = (-floor(N2/2):ceil(N2/2)-1).' + floor(136/2) + 1;

img_maxgirf1 = img_maxgirf1(idx1_range,:,:);
img_maxgirf2 = img_maxgirf2(idx1_range,:,:);
img_maxgirf3 = img_maxgirf3(idx1_range,:,:);
img_maxgirf4 = img_maxgirf4(idx1_range,:,:);
img_maxgirf5 = img_maxgirf5(idx1_range,:,:);
img_maxgirf6 = img_maxgirf6(idx1_range,:,:);

img_maxgirf_tse = img_maxgirf_tse(idx1_range,idx2_range,:);

x = x(idx1_range,:,:);
y = y(idx1_range,:,:);
z = z(idx1_range,:,:);

dx = dx(idx1_range,:,:);
dy = dy(idx1_range,:,:);
dz = dz(idx1_range,:,:);

%% Get directory information
dir_info = dir(fullfile(dicom_path, '*IMA'));
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

%--------------------------------------------------------------------------
% Number of slices
%--------------------------------------------------------------------------
N3 = nr_files;

%% Read a .dicom file
img_dicom_gnl1 = zeros(N1, N2, N3, 'double');

for idx = 1:nr_files
    %----------------------------------------------------------------------
    % Get a DICOM header
    %----------------------------------------------------------------------
    dicom_info = dicominfo(fullfile(dir_info(idx).folder, dir_info(idx).name));

    %----------------------------------------------------------------------
    % Read a DICOM image
    %----------------------------------------------------------------------
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
    img_dicom_gnl1(:,:,idx) = RescaleSlope * dicomread(dicom_info).' + RescaleIntercept; % transpose it!
end

%% Display reconstructed images
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

color_order_rgb = hex2rgb(color_order);

FontSize = 12;

%--------------------------------------------------------------------------
% Calculate a circle
%--------------------------------------------------------------------------
circle_radius = 49;
t = (0:0.01:1).';
x_circle = circle_radius * cos(2 * pi * t) + floor(N1/2) + 1;
y_circle = circle_radius * sin(2 * pi * t) + floor(N1/2) + 1;

slice_number = 9;

for slice_number = 1:nr_slices

figure('Color', 'w', 'Position', [0 2 1138 990]);

%--------------------------------------------------------------------------
% Magnitude: Grid./PHC/CFC/SFC/GNC/TOPUP = 0/0/0/0/0/0
%--------------------------------------------------------------------------
ax1 = subplot(4,4,1);
hold on;
imagesc(abs(img_maxgirf1(:,:,slice_number)));
plot(x_circle + 2, y_circle + 1, '--', 'Color', 'r', 'LineWidth', 1);
axis image ij off;
colormap(ax1, gray(256));
text(ax1, N1 / 2, -12, 'Cartesian MaxGIRF', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax1, N1 / 2,   0, {'no correction'}, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax1, 0, N2 / 2, 'Magnitude', 'Rotation', 90, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax1, 2, 0, '(A)', 'FontSize', 14, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
text(ax1, N1 / 2, N2, sprintf('PE direction (R >> L)'), 'FontSize', FontSize, 'Rotation', 0, 'Interpreter', 'tex', 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(ax1, N1 / 2, 0, sprintf('slice at y = %4.2fmm', y(1,1,slice_number) * 1e3), 'FontSize', FontSize, 'Rotation', 0, 'Interpreter', 'tex', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

%--------------------------------------------------------------------------
% Magnitude: Grid./PHC/CFC/SFC/GNC/TOPUP = 1/0/0/0/0/0
%--------------------------------------------------------------------------
ax2 = subplot(4,4,2);
hold on;
imagesc(abs(img_maxgirf2(:,:,slice_number)));
plot(x_circle + 2, y_circle + 1, '--', 'Color', 'r', 'LineWidth', 1);
axis image ij off;
colormap(ax2, gray(256));
text(ax2, N1 / 2, -12, 'Cartesian MaxGIRF', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax2, N1 / 2,   0, {'Gridding'}, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax2, 2, 0, '(B)', 'FontSize', 14, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Title
%--------------------------------------------------------------------------
text(ax2, N1, -66, '2D SE-EPI: coronal, 2.0 x 2.0 mm^2, R = 1, no partial Fourier, ETL = 128, 1 NSA', 'Color', blue, 'Interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

text(ax2, N1, -27, {sprintf('Gridding/PHC: gridding for ramp sampling/odd-even echo phase correction'), ...
                     sprintf('{\\color[rgb]{%f %f %f}CFC}/{\\color[rgb]{%f %f %f}GNC} vs {\\color[rgb]{%f %f %f}CFC}/{\\color[rgb]{%f %f %f}GNC(2D)}: concomitant field correction/gradient nonlinearity correction {\\color[rgb]{%f %f %f}during recon} vs {\\color[rgb]{%f %f %f}after recon}', ...
                     orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), green_siemens(1,1), green_siemens(1,2), green_siemens(1,3), green_siemens(1,1), green_siemens(1,2), green_siemens(1,3), ...
                     orange_siemens(1,1), orange_siemens(1,2), orange_siemens(1,3), green_siemens(1,1), green_siemens(1,2), green_siemens(1,3)), ...
                     'TOPUP: off-resonance correction using a displacement map from FSL TOPUP'}, 'Rotation', 0, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Create line (top)
annotation(gcf, 'line', [0.1830 0.9429]-0.16, [0.8510 0.8510], 'LineWidth', 2);

% Create line (bottom)
annotation(gcf, 'line', [0.1830 0.9429]-0.16, [0.7846 0.7846], 'LineWidth', 2);

%--------------------------------------------------------------------------
% Magnitude: Grid./PHC/CFC/SFC/GNC/TOPUP = 1/1/0/0/0/0
%--------------------------------------------------------------------------
ax3 = subplot(4,4,3);
hold on;
imagesc(abs(img_maxgirf3(:,:,slice_number)));
plot(x_circle + 2, y_circle + 1, '--', 'Color', 'r', 'LineWidth', 1);
axis image ij off;
colormap(ax3, gray(256));
text(ax3, N1 / 2, -12, 'Cartesian MaxGIRF', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax3, N1 / 2,   0, {'Gridding/PHC'}, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax3, 2, 0, '(C)', 'FontSize', 14, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Magnitude: Grid./PHC/CFC/SFC/GNC/TOPUP = 1/1/1/0/0/0
%--------------------------------------------------------------------------
ax4 = subplot(4,4,4);
hold on;
imagesc(abs(img_maxgirf4(:,:,slice_number)));
plot(x_circle + 2, y_circle + 1, '--', 'Color', 'r', 'LineWidth', 1);
axis image ij off;
colormap(ax4, gray(256));
text(ax4, N1 / 2, -12, 'Cartesian MaxGIRF', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax4, N1 / 2,   0, {'Gridding/PHC/{\color[rgb]{0.9255 0.3882 0.0}CFC}'}, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax4, 2, 0, '(D)', 'FontSize', 14, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Phase: Grid./PHC/CFC/SFC/GNC/TOPUP = 0/0/0/0/0/0
%--------------------------------------------------------------------------
ax5 = subplot(4,4,5);
imagesc(angle(img_maxgirf1(:,:,slice_number)) * 180 / pi);
axis image off;
colormap(ax5, hsv(256));
text(ax5, 0, N2 / 2, 'Phase', 'Rotation', 90, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax5, 2.8, 2.84, '(E)', 'FontSize', 14, 'FontWeight', 'Bold', 'Color', 'w', 'BackgroundColor', 'k', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Phase: Grid./PHC/CFC/SFC/GNC/TOPUP = 1/0/0/0/0/0
%--------------------------------------------------------------------------
ax6 = subplot(4,4,6);
imagesc(angle(img_maxgirf2(:,:,slice_number)) * 180 / pi);
axis image off;
colormap(ax6, hsv(256));
text(ax6, 2.8, 2.84, '(F)', 'FontSize', 14, 'FontWeight', 'Bold', 'Color', 'w', 'BackgroundColor', 'k', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Phase: Grid./PHC/CFC/SFC/GNC/TOPUP = 1/1/0/0/0/0
%--------------------------------------------------------------------------
ax7 = subplot(4,4,7);
imagesc(angle(img_maxgirf3(:,:,slice_number)) * 180 / pi);
axis image off;
colormap(ax7, hsv(256));
text(ax7, 2.8, 2.84, '(G)', 'FontSize', 14, 'FontWeight', 'Bold', 'Color', 'w', 'BackgroundColor', 'k', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Phase: Grid./PHC/CFC/SFC/GNC/TOPUP = 1/1/1/0/0/0
%--------------------------------------------------------------------------
ax8 = subplot(4,4,8);
imagesc(angle(img_maxgirf4(:,:,slice_number)) * 180 / pi);
axis image off;
colormap(ax8, hsv(256));
text(ax8, 2.8, 2.84, '(H)', 'FontSize', 14, 'FontWeight', 'Bold', 'Color', 'w', 'BackgroundColor', 'k', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Magnitude: Grid./PHC/CFC/SFC/GNC/TOPUP = 1/1/1/1/0/0
%--------------------------------------------------------------------------
ax9 = subplot(4,4,9);
hold on;
imagesc(abs(img_maxgirf5(:,:,slice_number)));
plot(x_circle + 2, y_circle + 1, '--', 'Color', 'r', 'LineWidth', 1);
axis image ij off;
colormap(ax9, gray(256));
text(ax9, N1 / 2, -12, 'Cartesian MaxGIRF', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax9, N1 / 2,   0, {'Grid./PHC/{\color[rgb]{0.9255 0.3882 0.0}CFC}/{\color[rgb]{0.9255 0.3882 0.0}GNC}'}, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax9, 0, N2 / 2, 'Magnitude', 'Rotation', 90, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax9, 2, 0, '(I)', 'FontSize', 14, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% EPI DICOM GNC=1
%--------------------------------------------------------------------------
ax10 = subplot(4,4,10);
hold on;
imagesc(img_dicom_gnl1(:,:,slice_number).');
plot(x_circle + 2, y_circle + 1, '--', 'Color', 'r', 'LineWidth', 1);
axis image ij off;
colormap(ax10, gray(256));
text(ax10, N1 / 2, -12, 'Traditional reconstruction', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax10, N1 / 2,   0, {'Grid./PHC/{\color[rgb]{0.0118 0.6000 0.6000}CFC}/{\color[rgb]{0.0118 0.6000 0.6000}GNC}'}, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax10, 2, 0, '(J)', 'FontSize', 14, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% Magnitude: Grid./PHC/CFC/SFC/GNC/TOPUP = 1/1/1/0/1/1
%--------------------------------------------------------------------------
ax11 = subplot(4,4,11);
hold on;
imagesc(abs(img_maxgirf6(:,:,slice_number)));
plot(x_circle + 2, y_circle + 1, '--', 'Color', 'r', 'LineWidth', 1);
axis image ij off;
colormap(ax11, gray(256));
text(ax11, N1 / 2, -12, 'Cartesian MaxGIRF', 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax11, N1 / 2,   0, {'Grid./PHC/{\color[rgb]{0.9255 0.3882 0.0}CFC}/{\color[rgb]{0.9255 0.3882 0.0}GNC}/{\color[rgb]{0.9255 0.3882 0.0}TOPUP}'}, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax11, 2, 0, '(K)', 'FontSize', 14, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

%--------------------------------------------------------------------------
% MaxGIRF TSE
%--------------------------------------------------------------------------
ax12 = subplot(4,4,12);
hold on;
imagesc(abs(flip(img_maxgirf_tse(:,:,slice_number).',1)));
plot(x_circle + 2, y_circle + 1, '--', 'Color', 'r', 'LineWidth', 1);
axis image ij off;
colormap(ax12, gray(256));
text(ax12, N1 / 2, 0, {'TSE reference w/ GNC'}, 'Color', 'k', 'Interpreter', 'tex', 'FontSize', FontSize, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(ax12, 2, 0, '(L)', 'FontSize', 14, 'FontWeight', 'Bold', 'Color', 'w', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

text(ax12, N1 / 2, 0     , 'H', 'FontSize', FontSize, 'Interpreter', 'tex', 'FontWeight', 'bold', 'Color', color_order_rgb(4,:), 'VerticalAlignment', 'top'   , 'HorizontalAlignment', 'center');
text(ax12, N1 / 2, N2    , 'F', 'FontSize', FontSize, 'Interpreter', 'tex', 'FontWeight', 'bold', 'Color', color_order_rgb(4,:), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
text(ax12, 1     , N2 / 2, 'R', 'FontSize', FontSize, 'Interpreter', 'tex', 'FontWeight', 'bold', 'Color', color_order_rgb(4,:), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
text(ax12, N1 - 1, N2 / 2, 'L', 'FontSize', FontSize, 'Interpreter', 'tex', 'FontWeight', 'bold', 'Color', color_order_rgb(4,:), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');

%--------------------------------------------------------------------------
% Singular values
%--------------------------------------------------------------------------
cutoff_percentage = 99.9;
S_scaled = S(:,slice_number) / S(1,slice_number);
cumulative_sum = cumsum(S_scaled) / sum(S_scaled);
indices = find(cumulative_sum > cutoff_percentage * 1e-2);
L = indices(1);

ax13 = axes;
plot(ax13, S_scaled, '.-', 'LineWidth', 1, 'MarkerSize', 12, 'Color', color_order{1});
axis square;
grid on;
grid minor;
set(ax13, 'YScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', FontSize);
xlabel('$$\ell$$', 'Interpreter', 'latex', 'FontSize', FontSize);
title(ax13, 'Singular values $$\sigma_{\ell}$$ for H', 'Interpreter', 'latex', 'FontSize', FontSize);
text(ax13, 0, 1, '(M)', 'FontSize', FontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

ax14 = axes;
hold on;
plot(ax14, cumulative_sum, '.-', 'LineWidth', 1, 'MarkerSize', 12, 'Color', color_order{1});
plot(ax14, [L L], [0 1.1], '--', 'Color', color_order{2}, 'LineWidth', 1);
axis square;
grid on;
grid minor;
set(ax14, 'Box', 'on', 'TickLabelInterpreter', 'latex', 'FontSize', FontSize);
ylim(ax14, [0 1.1]);
legend(ax14, '', sprintf('%5.2f$$\\%%$$ cutoff, L = %d', cutoff_percentage, L), 'Location', 'Southwest', 'Interpreter', 'latex', 'FontSize', 11);
xlabel(ax14, '$$\ell$$', 'Interpreter', 'latex', 'FontSize', FontSize);
title(ax14, {'Cumulative sum', '$$(\Sigma_{k=1}^{\ell}\sigma_{k}) / (\Sigma_{k=1}^{L_{\mathrm{max}}}\sigma_{k})$$'}, 'Interpreter', 'latex', 'FontSize', FontSize);
text(ax14, 0, 1.1, '(N)', 'FontSize', FontSize, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

set(ax13, 'Position', [0.8199 0.4384 0.1714 0.1970]);
set(ax14, 'Position', [0.8199 0.1394 0.1714 0.1970]);

set(ax1 , 'Position', [0.1502-0.16 0.7673-0.0454-0.2 0.2447 0.2181]);
set(ax2 , 'Position', [0.3432-0.16 0.7673-0.0454-0.2 0.2447 0.2181]);
set(ax3 , 'Position', [0.5365-0.16 0.7673-0.0454-0.2 0.2447 0.2181]);
set(ax4 , 'Position', [0.7299-0.16 0.7673-0.0454-0.2 0.2447 0.2181]);

set(ax5 , 'Position', [0.1502-0.16 0.5482-0.048-0.2 0.2447 0.2181]);
set(ax6 , 'Position', [0.3432-0.16 0.5482-0.048-0.2 0.2447 0.2181]);
set(ax7 , 'Position', [0.5365-0.16 0.5482-0.048-0.2 0.2447 0.2181]);
set(ax8 , 'Position', [0.7299-0.16 0.5482-0.048-0.2 0.2447 0.2181]);

set(ax9 , 'Position', [0.1502-0.16 0.2687-0.028-0.2 0.2447 0.2181]);
set(ax10, 'Position', [0.3432-0.16 0.2687-0.028-0.2 0.2447 0.2181]);
set(ax11, 'Position', [0.5365-0.16 0.2687-0.028-0.2 0.2447 0.2181]);
set(ax12, 'Position', [0.7299-0.16 0.2687-0.028-0.2 0.2447 0.2181]);

%export_fig(sprintf('figure2_slc%d', slice_number), '-r250', '-tif', '-c[300, 144, 94, 480]'); % [top,right,bottom,left]
export_fig(sprintf('figure2_slc%d', slice_number), '-r250', '-tif', '-c[300, 0, 94, 0]'); % [top,right,bottom,left]
close gcf;
end
