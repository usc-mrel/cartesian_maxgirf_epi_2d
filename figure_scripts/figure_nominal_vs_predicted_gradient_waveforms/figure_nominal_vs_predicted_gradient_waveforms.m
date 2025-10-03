% figure_nominal_vs_predicted_gradient_waveforms.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 07/20/2025, Last modified: 07/20/2025

%% Clean slate
close all; clear all; clc;

%% Start a stopwatch timer
start_time = tic;

%% Define the full path of an output path
% Phase encoding direction: R >> L
output_path = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828\meas_MID00057_FID21486_ep2d_se_bw1002_cor_RL_gridding0_phc0_cfc0_sfc0_gnc0_topup0';

slice_number = 15;
slice_type = 'flat';

%% Read a .cfl file
%--------------------------------------------------------------------------
% t_grt (grad_samples x 1)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 't_img_grt');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
t_grt = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% g_gcs_grt_nominal (3 x grad_samples)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'g_img_gcs_grt_nominal');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
g_gcs_grt_nominal = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% g_gcs_grt_predicted (3 x grad_samples)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'g_img_gcs_grt_predicted');
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
g_gcs_grt_predicted = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% k_gcs_grt_nominal (3 x grad_samples)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'k_img_gcs_grt_nominal');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
k_gcs_grt_nominal = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% k_grt_nominal (grad_samples x Nl)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'k_img_grt_nominal');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
k_grt_nominal = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% k_grt_predicted (grad_samples x Nl)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'k_img_grt_predicted');
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
k_grt_predicted = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% t_adc (adc_samples x 1)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 't_img_adc');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
t_adc = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% k_gcs_adc_nominal (3 x adc_samples)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'k_img_gcs_adc_nominal');
tstart = tic; fprintf('%s: Reading a .cfl file: %s... ', datetime, cfl_file);
k_gcs_adc_nominal = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% R_gcs2dcs (9 x 1)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, 'R_gcs2dcs');
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
R_gcs2dcs = reshape(readcfl(cfl_file), [3 3]);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% parabolic_shift (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('parabolic_shift_slc%d_%s', slice_number, slice_type));
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
parabolic_shift = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% parabolic_shift_nominal (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('parabolic_shift_nominal_slc%d_%s', slice_number, slice_type));
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
parabolic_shift_nominal = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%--------------------------------------------------------------------------
% parabolic_shift_predicted (Nkx x Nky x Nkz)
%--------------------------------------------------------------------------
cfl_file = fullfile(output_path, sprintf('parabolic_shift_predicted_slc%d_%s', slice_number, slice_type));
tstart = tic; fprintf('%s: Writing a .cfl file: %s... ', datetime, cfl_file);
parabolic_shift_predicted = readcfl(cfl_file);
fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));

%% Calculate nominal gradient waveforms in the DCS (GRT) [mT/m] [x,y,z]
g_dcs_grt_nominal = pagemtimes(R_gcs2dcs, g_gcs_grt_nominal);

%% Calculate nominal k-space trajectories in the DCS (GRT) [rad/m] [x,y,z]
k_dcs_grt_nominal = pagemtimes(R_gcs2dcs, k_gcs_grt_nominal);

%% Calculate nominal k-space trajectories in the DCS (ADC) [rad/m] [x,y,z]
k_dcs_adc_nominal = pagemtimes(R_gcs2dcs, k_gcs_adc_nominal);

%% Calculate GIRF-predicted gradient waveforms in the DCS (GRT) [mT/m] [x,y,z]
g_dcs_grt_predicted = pagemtimes(R_gcs2dcs, g_gcs_grt_predicted);

%% Remove readout oversampling
[Nkx,Nky] = size(parabolic_shift_nominal);
Nx = Nkx / 2;
Ny = Nky;

idx1_range = (-floor(Nx/2):ceil(Nx/2)-1).' + floor(Nkx/2) + 1;

parabolic_shift           = parabolic_shift(idx1_range,:);
parabolic_shift_nominal   = parabolic_shift_nominal(idx1_range,:);
parabolic_shift_predicted = parabolic_shift_predicted(idx1_range,:);

%% Dislay nominal vs predicted gradient waveforms
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

figure('Color', 'w', 'Position', [1 2 1086 990]);
ax1 = subplot(4,2,1);
hold on;
plot(t_grt * 1e3, g_dcs_grt_nominal(1,:), 'LineWidth', 1, 'Color', color_order{1});
plot(t_grt * 1e3, g_dcs_grt_predicted(1,:), 'LineWidth', 1, 'Color', color_order{2});
grid on;
grid minor;
set(gca, 'Box', 'On');
xlabel('Time [msec]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$G_{x}$ [mT/m]', 'Interpreter', 'latex', 'FontSize', 12);
title(ax1, '\textbf{(A)} Phase-encoding gradient waveform (physical x-axis)', 'Interpreter', 'latex', 'FontSize', 12);
legend(ax1, 'nominal', 'predicted', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'southeast');
ylim([-15 15]);

%--------------------------------------------------------------------------
% Title
%--------------------------------------------------------------------------
text(ax1, 260, 25, 'Comparision of nominal and GIRF-predicted gradient waveforms', 'Color', blue, 'Interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

text(ax1, 260, 20, '2D SE-EPI: coronal, 2.0 x 2.0 mm^2, R = 1, no partial Fourier, ETL = 128, 1 NSA', 'Color', blue, 'Interpreter', 'tex', 'FontSize', 14, 'FontWeight', 'normal', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

ax2 = subplot(4,2,2);
hold on;
plot(t_grt * 1e3, g_dcs_grt_nominal(1,:), 'LineWidth', 1, 'Color', color_order{1});
plot(t_grt * 1e3, g_dcs_grt_predicted(1,:), 'LineWidth', 1, 'Color', color_order{2});
plot(t_grt * 1e3, g_dcs_grt_nominal(1,:) - g_dcs_grt_predicted(1,:), 'LineWidth', 1, 'Color', color_order{3});
grid on;
grid minor;
set(gca, 'Box', 'On');
xlabel('Time [msec]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$G_{x}$ [mT/m]', 'Interpreter', 'latex', 'FontSize', 12);
title(ax2, '\textbf{(B)} Zoomed-in of (A)', 'Interpreter', 'latex', 'FontSize', 12);
legend(ax2, 'nominal', 'predicted', 'nominal - predicted', 'Interpreter', 'latex', 'FontSize', 10, 'Location', 'best');

xlim([81.6e-3 81.8e-3] * 1e3);
xlim([79.5e-3 84.5e-3] * 1e3);
%ylim([-1 3]);

ax3 = subplot(4,2,3);
hold on;
plot(t_grt * 1e3, g_dcs_grt_nominal(3,:), 'LineWidth', 1, 'Color', color_order{1});
plot(t_grt * 1e3, g_dcs_grt_predicted(3,:), 'LineWidth', 1, 'Color', color_order{2});
grid on;
grid minor;
set(gca, 'Box', 'On');
ylim([-15 15]);
xlabel('Time [msec]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$G_{z}$ [mT/m]', 'Interpreter', 'latex', 'FontSize', 12);
title(ax3, '\textbf{(C)} Readout gradient waveform (physical z-axis)', 'Interpreter', 'latex', 'FontSize', 12);

ax4 = subplot(4,2,4);
hold on;
plot(t_grt * 1e3, g_dcs_grt_nominal(3,:), 'LineWidth', 1, 'Color', color_order{1});
plot(t_grt * 1e3, g_dcs_grt_predicted(3,:), 'LineWidth', 1, 'Color', color_order{2});
plot(t_grt * 1e3, g_dcs_grt_nominal(3,:) - g_dcs_grt_predicted(3,:), 'LineWidth', 1, 'Color', color_order{3});
grid on;
grid minor;
set(gca, 'Box', 'On');
xlim([79.5e-3 84.5e-3] * 1e3);
%xlim([81.6e-3 81.8e-3] * 1e3);
ylim([-15 15]);
xlabel('Time [msec]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$G_{z}$ [mT/m]', 'Interpreter', 'latex', 'FontSize', 12);
title(ax4, '\textbf{(D)} Zoomed-in of (C)', 'Interpreter', 'latex', 'FontSize', 12);

ax5 = subplot(4,2,5);
hold on;
plot(t_grt * 1e3, k_grt_nominal(:,4), 'LineWidth', 1, 'Color', color_order{1});
plot(t_grt * 1e3, k_grt_predicted(:,4), 'LineWidth', 1, 'Color', color_order{2});
grid on;
grid minor;
set(gca, 'Box', 'On');
xlabel('Time [sec]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$k_{4}$ [rad/m$^2$]', 'Interpreter', 'latex', 'FontSize', 12);
text(65, 1000, '$$k_{4}(t) = \frac{1}{8 B_0}\int_{0}^{t} G_{z}^2(\tau) d\tau$$', 'Interpreter', 'latex', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'FontSize', 12);
title(ax5, '\textbf{(E)} 4th concomitant-field phase coefficient', 'Interpreter', 'latex', 'FontSize', 12);

ax6 = subplot(4,2,6);
hold on;
%plot(t_grt * 1e3, k_grt_nominal(:,4), 'LineWidth', 1, 'Color', color_order{1});
%plot(t_grt * 1e3, k_grt_predicted(:,4), 'LineWidth', 1, 'Color', color_order{2});
plot(t_grt * 1e3, k_grt_nominal(:,4) - k_grt_predicted(:,4), 'LineWidth', 1, 'Color', color_order{3});
grid on;
grid minor;
set(gca, 'Box', 'On');
xlim([79.5e-3 84.5e-3] * 1e3);
xlabel('Time [msec]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$k_{4}$ [rad/m$^2$]', 'Interpreter', 'latex', 'FontSize', 12);
title(ax6, '\textbf{(F)} Zoomed-in of (E)', 'Interpreter', 'latex', 'FontSize', 12);

ax7 = subplot(4,2,7);
hold on;
plot(t_grt * 1e3, k_grt_nominal(:,6), 'LineWidth', 1, 'Color', color_order{1});
plot(t_grt * 1e3, k_grt_predicted(:,6), 'LineWidth', 1, 'Color', color_order{2});
grid on;
grid minor;
set(gca, 'Box', 'On');
xlabel('Time [sec]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$k_{6}$ [rad/m$^2$]', 'Interpreter', 'latex', 'FontSize', 12);
text(65, 14.5, '$$k_{6}(t) = \frac{1}{2 B_0}\int_{0}^{t} (G_{x}^2(\tau) + G_{y}^2(\tau)) d\tau$$', 'Interpreter', 'latex', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'FontSize', 12);
title(ax7, '\textbf{(G)} 6th concomitant-field phase coefficient', 'Interpreter', 'latex', 'FontSize', 12);

ax8 = subplot(4,2,8);
hold on;
plot(t_grt * 1e3, k_grt_nominal(:,6) - k_grt_predicted(:,6), 'LineWidth', 1, 'Color', color_order{3});
grid on;
grid minor;
set(gca, 'Box', 'On');
xlim([79.5e-3 84.5e-3] * 1e3);
xlabel('Time [msec]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$k_{6}$ [rad/m$^2$]', 'Interpreter', 'latex', 'FontSize', 12);
title(ax8, '\textbf{(H)} Zoomed-in of (G)', 'Interpreter', 'latex', 'FontSize', 12);

ax9 = axes;
imagesc(ax9, parabolic_shift_nominal * 1e3);
axis image;
axis off;
colormap(ax9, turbo(256));
hc9 = colorbar;
title(ax9, {'\textbf{(I)} Parabolic shift', '(nominal)'}, 'Interpreter', 'latex', 'FontSize', 12);
set(hc9, 'Position', [0.9535 0.6374 0.0105 0.1909], 'TickLabelInterpreter', 'latex', 'FontSize', 10);
title(hc9, '[mm]', 'Interpreter', 'latex', 'FontSize', 10);

ax10 = axes;
imagesc(ax10, parabolic_shift_predicted * 1e3);
axis image;
axis off;
colormap(ax10, turbo(256));
hc10 = colorbar;
title(ax10, {'\textbf{(J)} Parabolic shift', '(predicted)'}, 'Interpreter', 'latex', 'FontSize', 12);
set(hc10, 'Position', [0.9535 0.3960 0.0105 0.1909], 'TickLabelInterpreter', 'latex', 'FontSize', 10);
title(hc10, '[mm]', 'Interpreter', 'latex', 'FontSize', 10);

ax11 = axes;
imagesc(ax11, (parabolic_shift_nominal - parabolic_shift_predicted) * 1e3);
axis image;
axis off;
colormap(ax11, turbo(256));
hc11 = colorbar;
title(ax11, {'\textbf{(K)} Difference', '(nominal - predicted)'}, 'Interpreter', 'latex', 'FontSize', 12);
set(hc11, 'Position', [0.9535 0.1525 0.0105 0.1909], 'TickLabelInterpreter', 'latex', 'FontSize', 10);
title(hc11, '[mm]', 'Interpreter', 'latex', 'FontSize', 10);

set(ax1, 'Position', [0.0830 0.7690 0.3580 0.1545]);
set(ax2, 'Position', [0.5028 0.7690 0.2541 0.1545]);

set(ax3, 'Position', [0.0830 0.5427 0.3580 0.1545]);
set(ax4, 'Position', [0.5028 0.5427 0.2541 0.1545]);

set(ax5, 'Position', [0.0830 0.3164 0.3580 0.1545]);
set(ax6, 'Position', [0.5028 0.3164 0.2541 0.1545]);

set(ax7, 'Position', [0.0830 0.0901 0.3580 0.1545]);
set(ax8, 'Position', [0.5028 0.0901 0.2541 0.1545]);

set(ax9 , 'Position', [0.7753 0.5504 0.1719 0.3668]);
set(ax10, 'Position', [0.7753 0.3070 0.1719 0.3668]);
set(ax11, 'Position', [0.7753 0.0636 0.1719 0.3668]);

export_fig(sprintf('figure_nominal_vs_predicted_gradient_waveforms'), '-r250', '-tif', '-c[0, 0, 130, 90]'); % [top,right,bottom,left]
