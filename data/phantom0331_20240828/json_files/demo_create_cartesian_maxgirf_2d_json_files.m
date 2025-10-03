% demo_create_cartesian_maxgirf_2d_json_files.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 06/29/2025, Last modified: 06/29/2025

%% Clean slate
close all; clear all; clc;

%% Define flags
server_flag = 0; %1=server

%% Define input filenames
%--------------------------------------------------------------------------
% data path
%--------------------------------------------------------------------------
if server_flag
    data_path = '/server/sdata/nlee/cartesian_maxgirf_epi_2d/data/phantom0331_20240828';
else
    data_path = 'D:\cartesian_maxgirf_epi_2d\data\phantom0331_20240828';
end

%--------------------------------------------------------------------------
% data filename (EPI)
%--------------------------------------------------------------------------
data_filenames{1} = 'meas_MID00057_FID21486_ep2d_se_bw1002_cor_RL';
data_filenames{2} = 'meas_MID00061_FID21490_ep2d_se_bw1002_cor_LR';

%--------------------------------------------------------------------------
% topup filename (EPI)
%--------------------------------------------------------------------------
topup_filenames{1} = 'topup_RL_LR';
topup_filenames{2} = 'topup_RL_LR';

%--------------------------------------------------------------------------
% bart_path
%--------------------------------------------------------------------------
if server_flag
    bart_path = '/server/home/nlee/bart';
else
    bart_path = '/home/image/bart';
end

%% Set reconstruction options
Lmax              = 20;         % maximum rank of the SVD approximation of a higher-order encoding matrix
cutoff_percentage = 99.9;       % Cutoff percentage to calculate an approximate full rank
lambda            = 0.0;        % regularization parameter for spatial domain
tol               = 1e-5;       % tolerance for PCG
maxiter           = 6;          % max. number of CG iterations
slice_type        = 'flat';     % slice type
cal_size          = [48 24 1];  % size of the calibration region

%--------------------------------------------------------------------------
% coil_selection
%--------------------------------------------------------------------------
coil_selection = [];

%--------------------------------------------------------------------------
% nr_slices
%--------------------------------------------------------------------------
nr_slices = 20;

%--------------------------------------------------------------------------
% nr_repetitions
%--------------------------------------------------------------------------
nr_repetitions = 1;

%--------------------------------------------------------------------------
% main_orientation (SAGITTAL/CORONAL/TRANSVERSAL = 0/1/2)
%--------------------------------------------------------------------------
main_orientation = 1;

%% Start a stopwatch timer
start_time = tic;

%% Calculate a list of flags
% gridding/phc/cfc/sfc/gnc/topup
flag_list = zeros(6, 6, 'double');

flag_list(1,:) = [0, 0, 0, 0, 0, 0]; % no correction
flag_list(2,:) = [1, 0, 0, 0, 0, 0]; % gridding
flag_list(3,:) = [1, 1, 0, 0, 0, 0]; % gridding/phc
flag_list(4,:) = [1, 1, 1, 0, 0, 0]; % gridding/phc/cfc
flag_list(5,:) = [1, 1, 1, 0, 1, 0]; % gridding/phc/cfc/gnc
flag_list(6,:) = [1, 1, 1, 0, 1, 1]; % gridding/phc/cfc/gnc/topup

%% Calculate the number of datasets
nr_datasets = length(data_filenames);

%% Calculate the number of .json files
nr_json_files = size(flag_list,1);

%% Create a .json file
first_indent_size = blanks(2);
second_indent_size = blanks(4);
third_indent_size = blanks(6);

for data_number = 1:nr_datasets
    %----------------------------------------------------------------------
    % Set filenames
    %----------------------------------------------------------------------
    data_filename  = data_filenames{data_number};
    topup_filename = topup_filenames{data_number};

    for json_number = 1:nr_json_files
        %------------------------------------------------------------------
        % Set flags
        %------------------------------------------------------------------
        gridding_flag = flag_list(json_number,1);
        phc_flag      = flag_list(json_number,2);
        cfc_flag      = flag_list(json_number,3);
        sfc_flag      = flag_list(json_number,4);
        gnc_flag      = flag_list(json_number,5);
        topup_flag    = flag_list(json_number,6);

        %------------------------------------------------------------------
        % Set the name of a .json file
        %------------------------------------------------------------------
        json_filename = sprintf('%s_gridding%d_phc%d_cfc%d_sfc%d_gnc%d_topup%d', data_filename, gridding_flag, phc_flag, cfc_flag, sfc_flag, gnc_flag, topup_flag);

        %------------------------------------------------------------------
        % Set the full path of an input filename
        %------------------------------------------------------------------
        siemens_twix_file  = fullfile(data_path, 'raw', sprintf('%s.dat', data_filename));
        ismrmrd_data_file  = fullfile(data_path, 'raw', sprintf('%s.h5', data_filename));
        ismrmrd_noise_file = fullfile(data_path, 'raw', sprintf('noise_%s.h5', data_filename));
        output_path        = fullfile(data_path, json_filename);
        sens_path          = fullfile(data_path, json_filename);
        topup_path         = fullfile(data_path, topup_filename);

        %------------------------------------------------------------------
        % Open a .json file
        %------------------------------------------------------------------
        if server_flag
            json_file = fullfile(pwd, sprintf('%s_server.json', json_filename));
        else
            json_file = fullfile(pwd, sprintf('%s.json', json_filename));
        end
        tstart = tic; fprintf('%s: Creating a .json file: %s... ', datetime, json_file);
        fid = fopen(json_file, 'w');

        %------------------------------------------------------------------
        % first line
        %------------------------------------------------------------------
        fprintf(fid, '{\n');

        %------------------------------------------------------------------
        % Set filenames
        %------------------------------------------------------------------
        fprintf(fid, '%s"siemens_twix_file": "%s",\n' , first_indent_size, strrep(siemens_twix_file, '\', '/'));
        fprintf(fid, '%s"ismrmrd_data_file": "%s",\n' , first_indent_size, strrep(ismrmrd_data_file, '\', '/'));
        fprintf(fid, '%s"ismrmrd_noise_file": "%s",\n', first_indent_size, strrep(ismrmrd_noise_file, '\', '/'));
        fprintf(fid, '%s"output_path": "%s",\n'       , first_indent_size, strrep(output_path, '\', '/'));
        fprintf(fid, '%s"sens_path": "%s",\n'         , first_indent_size, strrep(sens_path, '\', '/'));
        fprintf(fid, '%s"topup_path": "%s",\n'        , first_indent_size, strrep(topup_path, '\', '/'));
        fprintf(fid, '%s"bart_path": "%s",\n'         , first_indent_size, bart_path);

        %------------------------------------------------------------------
        % Recon parameters
        %------------------------------------------------------------------
        fprintf(fid, '%s"recon_parameters": [\n'          , first_indent_size);
        fprintf(fid, '%s{\n'                              , second_indent_size);

        if cfc_flag || sfc_flag
            fprintf(fid, '%s"Lmax": %d,\n'                , third_indent_size, Lmax);
            fprintf(fid, '%s"cutoff_percentage": %e,\n'   , third_indent_size, cutoff_percentage);
        else
            fprintf(fid, '%s"Lmax": %d,\n'                , third_indent_size, 1);
            fprintf(fid, '%s"cutoff_percentage": %e,\n'   , third_indent_size, 0);
        end

        fprintf(fid, '%s"lambda": %g,\n'                  , third_indent_size, lambda);
        fprintf(fid, '%s"tol": %g,\n'                     , third_indent_size, tol);
        fprintf(fid, '%s"maxiter": %d,\n'                 , third_indent_size, maxiter);
        fprintf(fid, '%s"slice_type": "%s",\n'            , third_indent_size, slice_type);
        fprintf(fid, '%s"cal_size": [%d, %d, %d],\n'      , third_indent_size, cal_size(1), cal_size(2), cal_size(3));
        fprintf(fid, '%s"remove_oversampling_flag": %d,\n', third_indent_size, 0);
        fprintf(fid, '%s"gridding_flag": %d,\n'           , third_indent_size, gridding_flag);
        fprintf(fid, '%s"phc_flag": %d,\n'                , third_indent_size, phc_flag);
        fprintf(fid, '%s"cfc_flag": %d,\n'                , third_indent_size, cfc_flag);
        fprintf(fid, '%s"sfc_flag": %d,\n'                , third_indent_size, sfc_flag);
        fprintf(fid, '%s"gnc_flag": %d,\n'                , third_indent_size, gnc_flag);
        fprintf(fid, '%s"topup_flag": %d\n'               , third_indent_size, topup_flag);
        fprintf(fid, '%s}\n'                              , second_indent_size);
        fprintf(fid, '%s],\n'                             , first_indent_size);

        %------------------------------------------------------------------
        % coil_selection
        %------------------------------------------------------------------
        if exist('coil_selection', 'var') && ~isempty(coil_selection)
            fprintf_text = sprintf('%%s"coil_selection": [ %s%s ],\n', repmat('%d, ', 1, length(coil_selection)-1), '%d');
            fprintf(fid, fprintf_text, first_indent_size, coil_selection);
        end

        %------------------------------------------------------------------
        % nr_slices
        %------------------------------------------------------------------
        if exist('nr_slices', 'var')
            fprintf(fid, '%s"nr_slices": %d,\n', first_indent_size, nr_slices);
        end

        %------------------------------------------------------------------
        % nr_repetitions
        %------------------------------------------------------------------
        if exist('nr_repetitions', 'var')
            fprintf(fid, '%s"nr_repetitions": %d,\n', first_indent_size, nr_repetitions);
        end

        %------------------------------------------------------------------
        % main_orientation
        %------------------------------------------------------------------
        if exist('main_orientation', 'var')
            fprintf(fid, '%s"main_orientation": %d\n', first_indent_size, main_orientation);
        end

        %------------------------------------------------------------------
        % last line
        %------------------------------------------------------------------
        fprintf(fid, '}\n');

        %------------------------------------------------------------------
        % Close a .json file
        %------------------------------------------------------------------
        fclose(fid);
        fprintf('done! (%6.4f/%6.4f sec)\n', toc(tstart), toc(start_time));
    end
end
