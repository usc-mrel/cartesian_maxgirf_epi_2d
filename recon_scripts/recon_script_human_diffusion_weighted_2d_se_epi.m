% recon_script_human_diffusion_weighted_2d_se_epi.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 07/01/2025, Last modified: 07/1/2025

%% Clean slate
close all; clear all; clc;

%% Set source directories
if ispc
    package_path   = 'D:\cartesian_maxgirf_epi_2d';
    ismrmrd_path   = 'D:\ismrmrd';
    grad_file_path = 'D:\cartesian_maxgirf_epi_2d\GradientCoils';
else
    package_path = '/server/sdata/nlee/cartesian_maxgirf_epi_2d';
    ismrmrd_path = '/server/sdata/nlee/ismrmrd';
    grad_file_path = '/server/sdata/nlee/cartesian_maxgirf_epi_2d/GradientCoils';
end

%% Add source directories to search path
addpath(genpath(package_path));
addpath(genpath(ismrmrd_path));

%% Define a list of .json files
if ispc
    json_files{1} = 'D:\cartesian_maxgirf_epi_2d\data\GNL_MPRAGE_EPI\json_files\meas_MID00152_FID61451_ep2d_diff_4Trace_p2_fatsat_gridding1_phc1_cfc0_sfc0_gnc0_topup0.json';
    json_files{2} = 'D:\cartesian_maxgirf_epi_2d\data\GNL_MPRAGE_EPI\json_files\meas_MID00152_FID61451_ep2d_diff_4Trace_p2_fatsat_gridding1_phc1_cfc1_sfc0_gnc0_topup0.json';
    json_files{3} = 'D:\cartesian_maxgirf_epi_2d\data\GNL_MPRAGE_EPI\json_files\meas_MID00152_FID61451_ep2d_diff_4Trace_p2_fatsat_gridding1_phc1_cfc1_sfc0_gnc1_topup0.json';
else
    json_files{1} = '/server/sdata/nlee/cartesian_maxgirf_epi_2d/data/GNL_MPRAGE_EPI/json_files/meas_MID00152_FID61451_ep2d_diff_4Trace_p2_fatsat_gridding1_phc1_cfc1_sfc0_gnc0_topup0_server.json';
    json_files{2} = '/server/sdata/nlee/cartesian_maxgirf_epi_2d/data/GNL_MPRAGE_EPI/json_files/meas_MID00152_FID61451_ep2d_diff_4Trace_p2_fatsat_gridding1_phc1_cfc1_sfc0_gnc1_topup0_server.json';
end

%% Calculate the number of json files
nr_json_files = length(json_files);

%% Process per json file
for json_number = 1%1:nr_json_files

    %% Define the name of a .json file
    json_file = json_files{json_number};

    %% Calculate voxel coordinates
    demo_cartesian_maxgirf_2d_calculate_voxel_coordinates;

    %% Prepare "calibration" k-space data
    demo_cartesian_maxgirf_2d_prepare_ksp_calibration;

    %% Prepare "imaging" k-space data
    demo_cartesian_maxgirf_2d_prepare_ksp_imaging;

    %% Estimate CSMs
    demo_cartesian_maxgirf_2d_estimate_csm;

    %% Cartesian MaxGIRF reconstruction
    demo_cartesian_maxgirf_2d_recon;

    %% Cartesian MaxGIRF reconstruction (low resolution)
    demo_cartesian_maxgirf_2d_recon_lowres;

    %% Phase-constrained, low-rank and subspace model-based MaxGIRF reconstruction
    demo_cartesian_maxgirf_2d_recon_pclr_diffusion;
end
