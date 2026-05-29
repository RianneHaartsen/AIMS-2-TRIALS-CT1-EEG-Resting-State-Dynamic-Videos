%% RCT1_01B: pre-processing - version B) 1/f (one-over-f)

% This script pre-processes the RCT1 data using the miniMADE pipeline with adjustments for the current project.
% The script takes the filtered continuous data from
% RCT1_01_PreprocFilterContinuous (paths in the table miniMADE_filtered_report) 
% and applies the next preprocessing steps, including:
% - segmenting stimuli data into epochs (2-sec duration with 50% overlap for FOOOF later)
% - demeaning of the data (baseline correction across the whole epoch)
% - artefact identification: blinks, drift, flat signal, exceeding
% threshold, jumps
% - epoch rejection/ interpolation
% - re-reference to average to remove noise common in all channels

% The output of the scripts are cleaned time series per channel per epoch,
% and an overview of all preprocessed datasets in the
% miniMADE_preprocessing_report table. 

% Created by dr. Rianne Haartsen, 2024

% Additions/ changes to miniMADE for RCT1:
% - 2-sec duration segements with 50% overlap
% - AR1: blinks and drift based on Luke Mason's approach + 2/3 channels for
% artefactnidentification
% - AR2: jump for 100µV within 4ms
% - Interpolation only for ≥ 20% of 19 channels otherwise epoch is excluded
% 2-sec epochs, 50% overlap, with 20% as threshold for interpolation
% based on the filtered data


% ************************************************************************
% The Maryland Analysis of Developmental EEG (MADE) Pipeline
% Version 1.1
% Developed at the Child Development Lab, University of Maryland, College Park

% Contributors to MADE pipeline:
% Ranjan Debnath (ranjan.ju@gmail.com)
% George A. Buzzell (georgebuzzell@gmail.com)
% Santiago Morales Pamplona (moraless@umd.edu)
% Stephanie C. Leach (sleach12@umd.edu)
% Maureen Elizabeth Bowers (mbowers1@umd.edu)
% Nathan A. Fox (fox@umd.edu)

% ----------------------------------------------------------------------- %
% Versions Log
%
% v1_1: Sept. 23, 2020 - modified by Stephanie C. Leach
%       -Fixed bugs arising from changes in EEGLab and Matlab functions
%           * pop_rmbase() no longer works with [] as the second argument.
%             The baseline removal code has been updated to change [] to 
%             an accepted argument in the baseline removal section
%       -Fixed a bug with the FASTER bad channel section
%           * removing the reference channel before removing the channels
%             marked as bad by FASTER causes an error if the reference
%             channel isn't the last row in the data matrix. We've now 
%             moved the reference channel removal to after bad channel 
%             removal (reference channel identified again)
%       -Added a flat channel check to the ICA prep (cleaning) code
%           * helps prevent ICA decompositions with less ICs than electrodes
%       -Changed ICA prep so that participants who lost more than 20% of
%        electrodes are saved as not having any usable data
%           * prevents crashing during automated IC rejection
%           * also prevents datasets with >20% channels interpolated
%       -Added a rank check before ADJUST or Adjusted-ADJUST
%           * avoids crashing when rank is < total number of channels
%       -Added option to run a version of MADE optimized for low-density EEG
%        systems (<32 channels)... called miniMADE
%           * skips the FASTER and ICA
%           * skipping interpolation is optional, but recommended
%           * modifies artifact rejection steps to include checks for flat
%             channels and voltage jumps
%           * advanced users have the option to replace bad channels with NaNs
%             instead of removing them from the code (NOTE: this will require
%             special code at later analysis steps to remove the NaNs before 
%             doing any more preprocessing steps or analyses)
%           * advanced users also have the option to use the user entered frontal
%             electrodes for additional epoch rejection before NaN replacement
%       -Added BIDS as a 3rd formatting option for saving data
%           * saves raw AND preprocessed data in BIDS format
%           * requires the user to fill out a few addition fields
% ----------------------------------------------------------------------- %

% MADE uses EEGLAB toolbox and some of its plugins. Before running the pipeline, you have to install the following:
% EEGLab:  https://sccn.ucsd.edu/eeglab/downloadtoolbox.php/download.php

% You also need to download the following plugins/extensions from here: https://sccn.ucsd.edu/wiki/EEGLAB_Extensions
% Specifically, download:
% MFFMatlabIO: https://github.com/arnodelorme/mffmatlabio/blob/master/README.txt
% FASTER: https://sourceforge.net/projects/faster/
% ADJUST: https://www.nitrc.org/projects/adjust/ AND Adjusted-ADJUST: included with this pipeline

% After downloading these plugins (as zip files), you need to place it in the eeglab/plugins folder.
% For instance, for FASTER, you uncompress the downloaded extension file (e.g., 'FASTER.zip') and place it in the main EEGLAB "plugins" sub-directory/sub-folder.
% After placing all the required plugins, add the EEGLAB folder to your path by using the following code:

% addpath(genpath(('...')) % Enter the path of the EEGLAB folder in this line

% Please cite the following references for in any manuscripts produced utilizing MADE pipeline:

% EEGLAB: A Delorme & S Makeig (2004) EEGLAB: an open source toolbox for
% analysis of single-trial EEG dynamics. Journal of Neuroscience Methods, 134, 9?21.

% firfilt (filter plugin): developed by Andreas Widmann (https://home.uni-leipzig.de/biocog/content/de/mitarbeiter/widmann/eeglab-plugins/)

% FASTER: Nolan, H., Whelan, R., Reilly, R.B., 2010. FASTER: Fully Automated Statistical
% Thresholding for EEG artifact Rejection. Journal of Neuroscience Methods, 192, 152?162.

% ADJUST: Mognon, A., Jovicich, J., Bruzzone, L., Buiatti, M., 2011. ADJUST: An automatic EEG
% artifact detector based on the joint use of spatial and temporal features. Psychophysiology, 48, 229?240.
%   Our group has modified ADJUST plugin to improve selection of ICA components containing artifacts.
%   If using our modified version, please cite the following reference.
%   Adjusted-ADJUST: Leach, S.C., Morales, S., Bowers, M. E., Buzzell, G. A., Debnath, R., Beall, D., 
%   Fox, N. A., 2020. Adjusting ADJUST: Optimizing the ADJUST Algorithm for Pediatric Data Using Geodesic 
%   Nets. Psychophysiology, 57(8), e13566.

% This pipeline is released under the GNU General Public License version 3.

% ************************************************************************

%% User input: user provide relevant information to be used for data processing
% Preprocessing of EEG data involves using some common parameters for
% every subject. This part of the script initializes the common parameters.


clear % clear matlab workspace
clc % clear matlab command window

% add path to RCT1 scripts
addpath('xxx/Rbaclofen_RCT1')
% add path to MADE pipeline
addpath(genpath('xxx/MADE-EEG-preprocessing-pipeline'))
% for example: addpath(genpath('C:\Users\Berger\Documents\MADE-EEG-preprocessing-pipeline-master'));

% addpath('/Users/riannehaartsen/Documents/MATLAB/eeglab2024.0');% enter the path of the EEGLAB folder in this line
addpath('xxx/eeglab2022.0/'); % use this path for most of the processing
% but use the 2024 version for reading in the layout
% for example: addpath(genpath('C:\Users\Berger\Documents\eeglab13_4_4b'));
eeglab nogui;% open eeglab

% Do you want to use miniMADE (recommended for low density (<32 channels) systems)
run_miniMADE = 1; % 0 = NO (run full MADE pipeline),  = YES (run MADE pipeline with minimal preprocessing steps)
% Note: Running miniMADE will skip the FASTER and ICA steps. Epoch level interpolation can still be performed, but is not recommended
% miniMADE also skips interim saving regardless of user selection

% 1. Enter the path of the folder that has the raw data to be analyzed
rawdata_location = 'xxx';

% 2. Enter the path of the folder where you want to save the processed data
output_location = 'xxx/DATA/01B_PreprocOOF';

% 3. Enter the path of the channel location file
channel_locations = 'xxx/EEG1010_RCT1.lay';

% 4. Do your data need correction for anti-aliasing filter and/or task related time offset?
adjust_time_offset = 0; % 0 = NO (no correction), 1 = YES (correct time offset)
% If your data need correction for time offset, initialize the offset time (in milliseconds)
filter_timeoffset   = 0; % anti-aliasing time offset (in milliseconds). 0 = No time offset
stimulus_timeoffset = 0; % stimulus related time offset (in milliseconds). 0 = No time offset
response_timeoffset = 0; % response related time offset (in milliseconds). 0 = No time offset
stimulus_markers = {'xxx', 'xxx'}; % enter the stimulus makers that need to be adjusted for time offset
respose_markers  = {'xxx', 'xxx'}; % enter the response makers that need to be adjusted for time offset

% 5. Do you want to down sample the data?
down_sample = 0; % 0 = NO (no down sampling), 1 = YES (down sampling)
sampling_rate = 500; % set sampling rate (in Hz), if you want to down sample

% 6. Do you want to delete the outer layer of the channels? (Rationale has been described in MADE manuscript)
%    This function can also be used to down sample electrodes. For example, if EEG was recorded with 128 channels but you would
%    like to analyse only 64 channels, you can assign the list of channnels to be excluded in the 'outerlayer_channel' variable.    
delete_outerlayer = 1; % 0 = NO (do not delete outer layer), 1 = YES (delete outerlayer);
% If you want to delete outer layer, make a list of channels to be deleted
outerlayer_channel = {'O2'}; % list of channels
% recommended list for EGI 128 channel net: {'E17' 'E38' 'E43' 'E44' 'E48' 'E49' 'E113' 'E114' 'E119' 'E120' 'E121' 'E125' 'E126' 'E127' 'E128' 'E56' 'E63' 'E68' 'E73' 'E81' 'E88' 'E94' 'E99' 'E107'}

% 7. Initialize the filters
highpass = 1; % High-pass frequency
lowpass  = 120; % Low-pass frequency. We recommend low-pass filter at/below line noise frequency (see manuscript for detail)

% 8. Are you processing task-related or resting-state EEG data?
task_eeg = 4; %RH:  0 = resting, 1 = task, 2 = resting state with 2 conditions, 4 = arbaclofen task with 4 rs/dynamic videos
taskonset_event_markers = {'10', '11', '17', '18', '1', '3'}; % enter all the event/condition markers; reflecting the onset of the condition - study specific

% 9. Do you want to epoch/segment your data?
epoch_data = 1; % 0 = NO (do not epoch), 1 = YES (epoch data)
rest_epoch_length = 1; % for resting EEG continuous data will be segmented into consecutive epochs of a specified length (here 1 second) by adding dummy events
taskoffset_event_markers = {'19', '19', '19', '19', '2', '4'}; % enter all the event/condition markers; reflecting the offset of the condition - study specific
overlap_epoch = 1;     % 0 = NO (do not create overlapping epoch), 1 = YES (50% overlapping epoch)
dummy_events ={'910','911','917', '918', '901', '903'}; % enter dummy events name

% 10. Do you want to remove/correct baseline?
remove_baseline = 1; % 0 = NO (no baseline correction), 1 = YES (baseline correction)
baseline_window = []; % baseline period in milliseconds (MS), [] = entire epoch

% 11. Do you want to remove artifact laden epoch based on voltage threshold?
voltthres_rejection = 1; % 0 = NO, 1 = YES
volt_threshold = [-120 120]; % lower and upper threshold (in uV)

% 12. Do you want to perform epoch level channel interpolation for artifact laden epoch? (see manuscript for detail)
% Note: interpolation is not recommended for systems with less than 20 channels
interp_epoch = 1; % 0 = NO, 1 = YES.
frontal_channels = {'AF8', 'AF7', 'Fpz'}; % If you set interp_epoch = 1, enter the list of frontal channels to check (see manuscript for detail)
% recommended list for EGI 128 channel net: {'E1', 'E8', 'E14', 'E21', 'E25', 'E32', 'E17'}
% recommended list for EGI 64 channel net: {'E1', 'E5', 'E10', 'E17'}

%13. Do you want to interpolate the bad channels that were removed from data?
% Note: because miniMADE automatically skips FASTER and ICA, this field will not affect miniMADE preprocessing
interp_channels = 1; % 0 = NO (Do not interpolate), 1 = YES (interpolate missing channels)

% 14. Do you want to rereference your data?
rerefer_data = 1; % 0 = NO, 1 = YES
reref=[]; % Enter electrode name/s or number/s to be used for rereferencing
% For channel name/s enter, reref = {'channel_name', 'channel_name'};
% For channel number/s enter, reref = [channel_number, channel_number];
% For average rereference enter, reref = []; default is average rereference

% 15. Do you want to save interim results?
save_interim_result = 0; % 0 = NO (Do not save) 1 = YES (save interim results)

% 16. How do you want to save your data? .set or .mat
output_format = 1; % 1 = .set (EEGLAB data structure), 2 = .mat (Matlab data structure), 3 = BIDS format
% If you chose BIDS format, specify subject number location in the file name and the task name
subject_number_loc = [1 2]; % should enter as [start stop] locations (e.g., par001_eeg.mff would be entered as [4 6])
task_name = 'task_name'; % should enter the eeg task name you want included in the file name

% ---------------- ADVANCED OPTIONS ---------------- %
% 17. Do you want to allow missing channels in epochs?
% Note: matlab matrices do not allow for missing rows (channels for data matrix). As such, channels removed from epochs will be replaced with
%       NaNs, which will need to be removed when averaging epochs (in the case of ERPs) or calculating other metrics (e.g., spectral power)
allow_missing_chans = 0; % this will replace bad channels with NaNs, 0 = NO and 1 = YES
% If allow_missing_chans = 1 (YES), volt_threshold values will be used to determine which channels are bad and will be replaced with NaN
% If allow_missing_chans = 1 (YES), interp_epoch & interp_channels CANNOT also = 1 (YES)
% If allow_missing_chans = 1 (YES), an average rereference cannot be used (reref cannot = [])
blink_check = 0; % this will check for blinks using the frontal_channels (defined in #12) and remove epochs containing them before replacing bad channels with NaNs
% This field will only be considered if allow_missing_chans = 1 (YES)
% WARNING: Make sure frontal_channels contains a list of frontal channels to check... If this variable is not properly defined the code will crash
chan_thresh = 0.8; % acceptable values are 0-1 and represent the percent of channels that must be good to keep an epoch
% for example: chan_thresh = 0.8 would remove epochs where greater than 80% of channels were replaced by NaNs

% ********* no need to edit beyond this point for EGI .mff data **********
% ********* for non-.mff data format edit data import function ***********
% ********* below using relevant data import plugin from EEGLAB **********

%% Read files to analyses
datafile_names=dir(rawdata_location);
datafile_names=datafile_names(~ismember({datafile_names.name},{'.', '..', '.DS_Store'}));
datafile_names=datafile_names(~contains({datafile_names.name},{'.fdt'}));
datafile_names={datafile_names.name};
[filepath,name,ext] = fileparts(char(datafile_names{1}));

%% Check whether EEGLAB and all necessary plugins are in Matlab path.
if exist('eeglab','file')==0
    error(['Please make sure EEGLAB is on your Matlab path. Please see EEGLAB' ...
        'wiki page for download and instalation instructions']);
end

if strcmp(ext, '.mff')==1
    if exist('mff_import', 'file')==0
        error(['Please make sure "mffmatlabio" plugin is in EEGLAB plugin folder and on Matlab path.' ...
            ' Please see EEGLAB wiki page for download and instalation instructions of plugins.' ...
            ' If you are not analysing EGI .mff data, edit the data import function below.']);
    end
else
    warning('Your data are not EGI .mff files. Make sure you edit data import function before using this script');
end

if exist('pop_firws', 'file')==0
    error(['Please make sure  "firfilt" plugin is in EEGLAB plugin folder and on Matlab path.' ...
        ' Please see EEGLAB wiki page for download and instalation instructions of plugins.']);
end

if exist('channel_properties', 'file')==0
    error(['Please make sure "FASTER" plugin is in EEGLAB plugin folder and on Matlab path.' ...
        ' Please see EEGLAB wiki page for download and instalation instructions of plugins.']);
end

if exist('ADJUST', 'file')==0
    error(['Please make sure you download modified "ADJUST" plugin from GitHub (link is in MADE manuscript)' ...
        ' and ADJUST is in EEGLAB plugin folder and on Matlab path.']);
end

%% Check that ADVANCED pipeline selections are compatible with other preprocessing selections
if allow_missing_chans == 1 && (interp_epoch == 1 || interp_channels == 1)
    error(['The allow_missing_chans option (ADVANCED) cannot be turned on if channel interpolation is on...' ...
        ' allow_missing_chans does not allow for channel interpolation. Please make sure interp_epoch and interp_channels are off']);
end

if allow_missing_chans == 1 && rerefer_data == 1 && isempty(reref)
    error(['An average rereference cannot be used if the allow_missing_chans option (ADVANCED) is on...' ...
        ' allow_missing_chans does not allow for average reference. Please ensure only a subset of channels are used for rereferencing']);
end

if allow_missing_chans == 1 && voltthres_rejection == 0
    warning('voltage threshold rejection thresholds will still be used to select bad channels and replace them with NaNs');
end

%% Create output folders to save data
if output_format < 3 % if not BIDS format
    if save_interim_result ==1
        if exist([output_location filesep 'filtered_data'], 'dir') == 0
            mkdir([output_location filesep 'filtered_data'])
        end
        if exist([output_location filesep 'ica_data'], 'dir') == 0
            mkdir([output_location filesep 'ica_data'])
        end
    end
    if exist([output_location filesep 'processed_data'], 'dir') == 0
        mkdir([output_location filesep 'processed_data'])
    end
elseif output_format == 3 % if BIDS format
    % check if derivatives folder already exists and create it if not (where we will save preprocessed data)
    if exist([ output_location filesep 'derivatives']) == 0
        mkdir([ output_location filesep 'derivatives'])
    end
    % check if eegpreprocess folder already exists and create it if not (where we will save preprocessed data)
    if exist([ output_location filesep 'derivatives' filesep 'eegpreprocess']) == 0
        mkdir([ output_location filesep 'derivatives' filesep 'eegpreprocess'])
    end
    % create file path for derivatives folder
    output_location_derivatives = [output_location filesep 'derivatives'];
end

%% Initialize output variables

% epoch rejection
total_epochs_before_artifact_rejection=[];
total_epochs_after_artifact_rejection=[];
% light sensor check
LS_check = [];
% marker presence
N_TaskOnsetMarkers = [];
% AR info
AR1_blinks = [];
Neps_postAR1 = [];
AR2a_thresholds = [];
AR2b_flat = [];
AR2c_jumps = [];
AR2_thr_flat_jump = []; % contains the channels interpolated per trial
Eps_BAD_InvalidInterp = []; % trials with more thatn 20% channels interpolated
Neps_postAR2 = [];
Chan_labels_curr = [];

% channels interpolated
% total_channels_interpolated=[]; % total_channels_interpolated=faster_bad_channels+ica_preparation_bad_channels
curdate=datestr(now,'dd-mm-yyyy'); % set current date here so that table won't crash if date changes

% load the preprocessing filtering table
load xxx/DATA/01_Preproc_Filtered_Continuous/miniMADE_filtered_report.mat
filtered_location = 'xxx/DATA/01_Preproc_Filtered_Continuous';

%% Loop over all data files
for subject = 1:height(miniMADE_filtered_table)
    EEG=[];
    
    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', subject, miniMADE_filtered_table.datafile_names{subject});

        %% Initialize output variables
    
        % epoch rejection
        total_epochs_before_artifact_rejection(subject)=0;
        total_epochs_after_artifact_rejection(subject)=0;
        % light sensor check
        LS_check = miniMADE_filtered_table.LS_check(subject);
        % marker presence
        N_TaskOnsetMarkers = [];
        % AR info
        AR1_badFrontCh = [];
        AR1_blinks = [];
        Neps_postAR1 = 0;
        AR2a_thresholds = [];
        AR2b_flat = [];
        AR2c_jumps = [];
        AR2_thr_flat_jump = []; % contains the channels interpolated per trial
        Eps_BAD_InvalidInterp = []; % trials with more thatn 20% channels interpolated
        Neps_postAR2 = 0;
        Chan_labels_curr = [];
        
        % channels interpolated
        % total_channels_interpolated=[]; % total_channels_interpolated=faster_bad_channels+ica_preparation_bad_channels
        curdate=datestr(now,'dd-mm-yyyy'); % set current date here so that table won't crash if date changes
        savedDataPath = NaN; 

        cd(output_location)
        if exist('miniMADE_preprocessing_report.mat','file')
            load miniMADE_preprocessing_report.mat
        end

        report_table_newrow=table(miniMADE_filtered_table.datafile_names(subject)', {curdate}, {savedDataPath}, total_epochs_before_artifact_rejection(subject)', total_epochs_after_artifact_rejection(subject)',...
            LS_check, {N_TaskOnsetMarkers}, {AR1_badFrontCh}, {AR1_blinks}, Neps_postAR1, {AR2a_thresholds}, {AR2b_flat}, {AR2c_jumps}, {AR2_thr_flat_jump}, {Eps_BAD_InvalidInterp}, Neps_postAR2, {Chan_labels_curr});
        report_table_newrow.Properties.VariableNames={'datafile_names', 'CurDate',  'DataPath', 'total_epochs_before_artifact_rejection', 'total_epochs_after_artifact_rejection',...
            'LS','EEGevents','AR1_badfrontCh','AR1_blinks','Neps_postAR1','AR2a_thres','AR2b_flat','AR2c_jump','AR2_summary','Eps_toomuchinterpolation','Neps_postAR2','Chan_labels'};
        miniMADE_report_table(subject,:) = report_table_newrow;
        save('miniMADE_preprocessing_report.mat','miniMADE_report_table')

    try
        %% STEP 1: Import EGI data file and relevant information
        % EEG = mff_import([rawdata_location filesep datafile_names{subject}]);
        EEG = pop_loadset([miniMADE_filtered_table.DataPath{subject}]);
        EEG = eeg_checkset(EEG);

        %% STEP 12: Segment data into fixed length epochs
        if epoch_data==1
            if task_eeg ==1 % task eeg
                EEG = eeg_checkset(EEG);
                EEG = pop_epoch(EEG, task_event_markers, task_epoch_length, 'epochinfo', 'yes');
            elseif task_eeg==0 % resting eeg
                if overlap_epoch==1
                    EEG=eeg_regepochs(EEG,'recurrence',(rest_epoch_length/2),'limits',[0 rest_epoch_length], 'rmbase', [NaN], 'eventtype', char(dummy_events));
                    EEG = eeg_checkset(EEG);
                else
                    EEG=eeg_regepochs(EEG,'recurrence',rest_epoch_length,'limits',[0 rest_epoch_length], 'rmbase', [NaN], 'eventtype', char(dummy_events));
                    EEG = eeg_checkset(EEG);
                end
    
            elseif task_eeg==4 % resting & dynamic videos eeg - 4 conditions
    
                    % Often resting state EEG is collected in eyes close and eyes open conditions
                    % This script inserts event at specific time interval for resting state EEG
                    % and creates separate events for eyes close and eyes open resting EEG
                    % data. See MADE paper for details.
    
                    % Initialize variables
    
                    % 1. Name of eyes close and eyes open markers
                    rest_event_markers = taskonset_event_markers; %RH: video markers
                    % 2. Name of new eyes close and eyes open markers
                    new_rest_markers = dummy_events; % enter markers for epochs
                    % 3. Does the data have a trial end marker?
                    trial_end_marker = 1; % 0=NO (no trial end marker in data), 1=YES (data have trial end marker)
                    trial_end_marker_name= taskoffset_event_markers; % enter trial end marker name
                    % 4. Length of epochs in seconds
                    rest_epoch_length=2;
                    % 5. Do you want to create overlapping epoch?
                    overlap_epoch = 1; % 0 = NO (do not create overlapping epoch), 1 = YES (50% overlapping epoch)
    
                    % Insert markers
                    if overlap_epoch==1
                        time_samples = (rest_epoch_length/2)*EEG.srate; % convert time window into samples or data points
                    else
                        time_samples = rest_epoch_length*EEG.srate;
                    end
                    
                    EEG.urevent_all = EEG.event; 
                    EEG.urevent = EEG.event;
                    EEG.event = [];

                    tm=1;
                    for ue=1:length(EEG.urevent)
                        for rm=1:length(rest_event_markers)
                            if EEG.urevent(ue).type == str2double(rest_event_markers{rm})
                                if trial_end_marker == 1
                                    for te=ue:length(EEG.urevent)
                                        if EEG.urevent(te).type == str2double(trial_end_marker_name{rm})
                                            trial_end_latency = EEG.urevent(te).latency;
                                            break;
                                        end
                                    end
                                else
                                    if ue < length(EEG.urevent)
                                        trial_end_latency=EEG.urevent(ue+1).latency-rest_epoch_length-time_samples;
                                    else
                                        trial_end_latency=length(EEG.times)-rest_epoch_length-time_samples;
                                    end
                                end
                                event_times = EEG.urevent(ue).latency;
                                while (event_times+time_samples) < trial_end_latency
                                    EEG.event(tm).type = str2double(new_rest_markers{rm});
                                    EEG.event(tm).latency = event_times;
                                    event_times = event_times+time_samples;
                                    tm=tm+1;
                                end
                            end
                        end
                    end
                    
                    % adjust the new EEG.events.value and duration
                    for eps = 1:length(EEG.event)
                            EEG.event(eps).value = 'EEGdummy';
                            EEG.event(eps).duration = 1;
                    end
                    % adjust the latencies to round numbers
                    for eps = 1:length(EEG.event)
                            EEG.event(eps).latency = round(EEG.event(eps).latency);
                    end
                    
                    % restore and adjust the EEG.urevent 
                        Nlasturevent = length(EEG.urevent_all);
                        NewUrevents = EEG.urevent_all;
                        NewEvents = EEG.event;
                        for eps = 1:length(EEG.event)
                            % add event into EEG.urevent
                            NewUrevents(Nlasturevent + eps).type = NewEvents(eps).type;
                            NewUrevents(Nlasturevent + eps).latency = NewEvents(eps).latency;
                            NewUrevents(Nlasturevent + eps).duration = NewEvents(eps).duration;
                            NewUrevents(Nlasturevent + eps).urevent = Nlasturevent + eps;
                            % add in index urevent_orig into event field
                            EEG.event(eps).urevent = Nlasturevent + eps;
                        end
                        EEG.urevent = NewUrevents;      

                    % create epoch
                        EEG = eeg_checkset(EEG);
                        EEG = pop_epoch(EEG, new_rest_markers, [0 rest_epoch_length], 'epochinfo', 'yes');
                        EEG = eeg_checkset(EEG);

                        clear tm ue rm trial_end_marker rest_event_markers new_rest_markers event_times trial_end+latency
                        clear Nlasturevent NewUrevents NewEvents
                end
                    
        end
        
        total_epochs_before_artifact_rejection(subject)=EEG.trials;
        
        % check onset markers present for later
        N_TaskOnsetMarkers = zeros(2,size(taskonset_event_markers,2));
        for ss = 1:size(taskonset_event_markers,2)
            curr_marker_onset = str2double(taskonset_event_markers{1,ss});
            curr_offset = str2double(taskoffset_event_markers{1,ss});
            count = 0;
             for rr = 1:size(EEG.urevent,1)
                if EEG.urevent(rr).type == curr_marker_onset && EEG.urevent(rr+1).type == curr_offset
                    count = count + 1;
                elseif EEG.urevent(rr).type == curr_marker_onset && EEG.urevent(rr+1).type == 5 && EEG.urevent(rr+2).type == curr_offset
                    count = count + 1;
                elseif EEG.urevent(rr).type == curr_marker_onset && EEG.urevent(rr+1).type == 6 && EEG.urevent(rr+2).type == curr_offset
                    count = count + 1;
                end
            end
            N_TaskOnsetMarkers(:,ss) = [curr_marker_onset; count];
        end
        clear ss curr_marker_onset curr_offset count

        % updata miniMADE report table
        report_table_newrow=table(datafile_names(subject)', {curdate}, {savedDataPath}, total_epochs_before_artifact_rejection(subject)', total_epochs_after_artifact_rejection(subject)',...
            LS_check, {N_TaskOnsetMarkers}, {AR1_badFrontCh}, {AR1_blinks}, Neps_postAR1, {AR2a_thresholds}, {AR2b_flat}, {AR2c_jumps}, {AR2_thr_flat_jump}, {Eps_BAD_InvalidInterp}, Neps_postAR2, {Chan_labels_curr});
        report_table_newrow.Properties.VariableNames={'datafile_names', 'CurDate',  'DataPath', 'total_epochs_before_artifact_rejection', 'total_epochs_after_artifact_rejection',...
            'LS','EEGevents','AR1_badfrontCh','AR1_blinks','Neps_postAR1','AR2a_thres','AR2b_flat','AR2c_jump','AR2_summary','Eps_toomuchinterpolation','Neps_postAR2','Chan_labels'};
        miniMADE_report_table(subject,:) = report_table_newrow;
        save('miniMADE_preprocessing_report.mat','miniMADE_report_table')

        %% STEP 13: Remove baseline
        if remove_baseline==1
            if isempty(baseline_window) % set up for entire epoch
                baseline_window = [EEG.times(1) EEG.times(end)]; % set start and stop times for the baseline window
            end
            EEG = eeg_checkset( EEG );
            EEG = pop_rmbase( EEG, baseline_window);
        end
        
        %% STEP 14: Artifact rejection

        for cc = 1:size(EEG.chanlocs,2)
            Chan_labels_curr{1,cc} = EEG.chanlocs(cc).labels;
        end

        all_bad_epochs=0;
        if allow_missing_chans == 0 
            if voltthres_rejection==1 % check voltage threshold rejection
                if interp_epoch==1 % check epoch level channel interpolation
                    % first pass: loop through frontal channels and reject bad epochs

                    %   - removes epochs with (possible) blinks - based on
                    %   filtering and model fitting
                    blinkLen = 0.05;
                    maxsd = 2.5;
                    maxr2 = 0.6;
                    EEG_foreog = EEG;

                    % get number of chans/trials
                    numChans = size(EEG_foreog.data, 1);
                    numTrials = EEG_foreog.trials;

                    % Filter data 3-10Hz: as previously in this pipeline  
                    highpass_eog = 3; lowpass_eog = 10;

                    high_transband_eog = highpass_eog; % high pass transition band
                    low_transband_eog = 10; % low pass transition band
                    hp_eog_fl_order = 3.3 / (high_transband_eog / EEG.srate);
                    lp_eog_fl_order = 3.3 / (low_transband_eog / EEG.srate);
                    
                    % Round filter order to next higher even integer. Filter order is always even integer.
                    if mod(floor(hp_eog_fl_order),2) == 0
                        hp_eog_fl_order=floor(hp_eog_fl_order);
                    elseif mod(floor(hp_eog_fl_order),2) == 1
                        hp_eog_fl_order=floor(hp_eog_fl_order)+1;
                    end
                    
                    if mod(floor(lp_eog_fl_order),2) == 0
                        lp_eog_fl_order=floor(lp_eog_fl_order)+2;
                    elseif mod(floor(lp_eog_fl_order),2) == 1
                        lp_eog_fl_order=floor(lp_eog_fl_order)+1;
                    end
                    
                    % Calculate cutoff frequency
                    high_cutoff_eog = highpass_eog/2;
                    low_cutoff_eog = lowpass_eog + (low_transband_eog/2);

                    
                    % Performing high pass filtering
                    EEG_foreog = eeg_checkset( EEG_foreog );
                    EEG_foreog = pop_firws(EEG_foreog, 'fcutoff', high_cutoff_eog, 'ftype', 'highpass', 'wtype', 'hamming', 'forder', hp_eog_fl_order, 'minphase', 0);
                    EEG_foreog = eeg_checkset( EEG_foreog );
                    
                    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                    
                    % pop_firws() - filter window type hamming ('wtype', 'hamming')
                    % pop_firws() - applying zero-phase (non-causal) filter ('minphase', 0)
                    
                    % Performing low pass filtering
                    EEG_foreog = eeg_checkset( EEG_foreog );
                    EEG_foreog = pop_firws(EEG_foreog, 'fcutoff', low_cutoff_eog, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder', lp_eog_fl_order, 'minphase', 0);
                    EEG_foreog = eeg_checkset( EEG_foreog );

    
                    % compute channel z-scores for each trial
                    lens = repmat(size(EEG_foreog.data,2),[1,size(EEG_foreog.data,3)]);
                    for tt = 1:size(EEG_foreog.data,3)
                        if tt == 1
                            cont = EEG_foreog.data(:,:,1);
                        else
                            cont = [cont EEG_foreog.data(:,:,tt)];
                        end
                    end
                    % calculate z-scores for each channel across samples
                    zcont = zscore(cont, [], 2);
                    zData = EEG_foreog;
                    for tr = 1:zData.trials
                        s1 = 1 + ((tr - 1) * lens(tr));
                        s2 = tr * lens(tr);
                        zData.data(:,:,tr) = zcont(:, s1:s2);
                    end
                    % find samples exceeding maxsd
                    zcrit = zData.data > maxsd;
                    
                    % % compute channel zscores for each trial
                    % zdata = eegZScoreSegs(data_blink);
                    % zcrit = cellfun(@(x) x > maxsd, zdata.trial, 'uniform', false); 
                    % separate matrices
                    blink = false(numChans, numTrials);
                    drift = false(numChans, numTrials);

                    % get indices for frontal channels
                    chans=[]; chansidx=[];chans_labels2=[];
                    chans_labels2=cell(1,EEG.nbchan);
                    for i=1:EEG.nbchan
                        chans_labels2{i}= EEG.chanlocs(i).labels;
                    end
                    [chans,chansidx] = ismember(frontal_channels, chans_labels2);
                    frontal_channels_idx = chansidx(chansidx ~= 0);

                    for tr = 1:numTrials
                        for ch = 1:numChans
                            % check excluded electrodes
                            if any(ismember(frontal_channels_idx,ch))
                                % get time range
                                s1 = 1; s2 = length(zData.times); 
                                % detect blink artefacts                
                                ct = findcontig2(zcrit(ch, s1:s2, tr)', 1);
                                if ~isempty(ct)
                                    len = ct(:, 3) / zData.srate;
                                    blink(ch, tr) = any(len > blinkLen);
                                end             
                                % detect drift
                                if ~blink(ch, tr)
                                    [~, gof] = fit(EEG_foreog.times', double(EEG.data(ch, :, tr)'),...
                                        'poly2');
                                    drift(ch, tr) = gof.rsquare >= .65;
                                    clear gof
                                end
                            end
                        end
                    end
                    
                    % collate for tracking
                    AR1_badFrontCh = [blink(frontal_channels_idx,:); drift(frontal_channels_idx,:)];
                    badepoch=zeros(1, EEG.trials);
                        for ii=1:size(AR1_badFrontCh, 2)
                            bad_blink = sum(AR1_badFrontCh([1:3],ii),1);
                            bad_drift = sum(AR1_badFrontCh([4:6],ii),1);
                            if bad_blink >= 2 
                                badepoch(ii)= 1;
                            elseif bad_drift >= 2 
                                badepoch(ii)= 1;
                            elseif bad_blink >= 2 && bad_drift >= 2
                                badepoch(ii)= 1;
                            end
                        end
                        badepoch=logical(badepoch);
                    AR1_blinks = badepoch;

                    clear blinkLen maxsd maxr2 EEG_foreog numChans numTrials
                    clear highpass_eog lowpass_eog low_transband high_transband low_transband
                    clear hp_eog_fl_order lp_eog_fl_order high_cutoff low_cutoff 
                    clear lens tt cont zcont zData zcrit tr s1 s2 blink drift tr ch s1 s2 lent
    
                    % If all epochs are artifacted, save the dataset and ignore rest of the preprocessing for this subject.
                    if sum(badepoch)==EEG.trials || sum(badepoch)+1==EEG.trials
                        all_bad_epochs=1;
                        warning(['No usable data for datafile', datafile_names{subject}]);
                        Neps_postAR1 = 0;
                    else
                        EEG = pop_rejepoch( EEG, badepoch, 0);
                        EEG = eeg_checkset(EEG);
                        Neps_postAR1 = EEG.trials;
                    end
    
                    % second pass: loop through all channels and interpolate remaining bad channels at the epoch level
                    %   - miniMADE has extra artifact checks at this step
                    if all_bad_epochs==1
                        warning(['No usable data for datafile', datafile_names{subject}]);
                    else
                        % Interpolate artifacted data for all remaing channels
                        badChans = zeros(EEG.nbchan, EEG.trials);
                        % Find artifacted epochs by detecting outlier voltage but don't remove
                        for ch=1:EEG.nbchan
                            EEG = pop_eegthresh(EEG,1, ch, volt_threshold(1), volt_threshold(2), EEG.xmin, EEG.xmax,0,0);
                            EEG = eeg_checkset(EEG);
                            EEG = eeg_rejsuperpose(EEG, 1, 1, 1, 1, 1, 1, 1, 1);
                            badChans(ch,:) = EEG.reject.rejglobal;
                        end
                        AR2a_thresholds = badChans; 

                        AR2b_flat = zeros(EEG.nbchan, EEG.trials);
                        AR2c_jumps = zeros(EEG.nbchan, EEG.trials);
                        tmpData = zeros(EEG.nbchan, EEG.pnts, EEG.trials);
                        if run_miniMADE == 0
                            for e = 1:EEG.trials
                                % Initialize variables EEGe and EEGe_interp;
                                EEGe = []; EEGe_interp = []; badChanNum = [];
                                % Select only this epoch (e)
                                EEGe = pop_selectevent( EEG, 'epoch', e, 'deleteevents', 'off', 'deleteepochs', 'on', 'invertepochs', 'off');
                                badChanNum = find(badChans(:,e)==1); % find which channels are bad for this epoch
                                EEGe_interp = eeg_interp(EEGe,badChanNum); %interpolate the bad channels for this epoch
                                tmpData(:,:,e) = EEGe_interp.data; % store interpolated data into matrix
                            end
                        elseif run_miniMADE == 1
                            for e = 1:EEG.trials
                                EEGe = []; EEGe_interp = []; badChanNum = []; % Initialize variables EEGe and EEGe_interp;
                                %select only this epoch (e)
                                EEGe = pop_selectevent( EEG, 'epoch',e,'deleteevents','off','deleteepochs','on','invertepochs','off');
                                badChanNum = find(badChans(:,e)==1); %find which channels are bad for this epoch

                                % find and add flat chans to the bad chans list
                                flatChanNum = find(range(EEGe.data,2) < 1);
                                badChanNum  = unique([badChanNum; flatChanNum]);

                                % find chans with a large jump/deflection
                                % by taking the difference between 2
                                % samples 4ms apart
                                Timepersample = 1/EEGe.srate; Steps = .004/Timepersample;
                                % create vector with differences
                                Npoints = length(EEGe.data);
                                differences = zeros(size(EEGe.data,1), Npoints - Steps);
                                for ii = 1:(Npoints - Steps)
                                    differences(:,ii) = EEGe.data(:,(ii + Steps)) - EEGe.data(:,ii);
                                end
                                % find channels with jumps exceeding
                                % 100µV
                                [jump_chans, ~] = find(abs(differences) > 100); 
                                badChanNum = unique([badChanNum; unique(jump_chans)]);

                                % interpolate using bad channel list with extra checks
                                if length(badChanNum) < EEGe.nbchan - 1 % script will crash if we try to interpolate with 0 or 1 channels left 
                                    EEGe_interp = eeg_interp(EEGe,badChanNum); %interpolate the bad channels for this epoch
                                    tmpData(:,:,e) = EEGe_interp.data; % store interpolated data into matrix
                                end
                                badChans(badChanNum,e) = 1; % modify

                                % keep track of flat and jump channel information
                                if ~isempty(flatChanNum)
                                    AR2b_flat(flatChanNum,e) = 1;
                                end
                                if ~isempty(jump_chans)
                                    AR2c_jumps(jump_chans,e) = 1;
                                end
                                AR2_thr_flat_jump = badChans;
                            end
                        end
                        EEG.data = tmpData; % now that all of the epochs have been interpolated, write the data back to the main file
    
                        % If more than 20% of channels in an epoch were interpolated, reject that epoch
                        badepoch=zeros(1, EEG.trials);
                        for ei=1:EEG.trials
                            NumbadChan = badChans(:,ei); % find how many channels are bad in an epoch
                            if sum(NumbadChan) > round((20/100)* size(chans_labels2,2))% check if more than 20% of 19 channels are bad
                                badepoch (ei)= sum(NumbadChan);
                            end
                        end
                        badepoch=logical(badepoch);
                        Eps_BAD_InvalidInterp = badepoch;
                    end
                    % If all epochs are artifacted, save the dataset and ignore rest of the preprocessing for this subject.
                    if sum(badepoch)==EEG.trials || sum(badepoch)+1==EEG.trials
                        all_bad_epochs=1;
                        warning(['No usable data for datafile', datafile_names{subject}]);
                        Neps_posAR2 = 0;
                    else
                        EEG = pop_rejepoch(EEG, badepoch, 0);
                        EEG = eeg_checkset(EEG);
                        Neps_postAR2 = EEG.trials;
                    end

                else % if no epoch level channel interpolation
                    if run_miniMADE == 0
                        EEG = pop_eegthresh(EEG, 1, (1:EEG.nbchan), volt_threshold(1), volt_threshold(2), EEG.xmin, EEG.xmax, 0, 0);
                        EEG = eeg_checkset(EEG);
                        EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
                    elseif run_miniMADE == 1
                        badChans = zeros(EEG.nbchan, EEG.trials);
                        % Find artifacted epochs by detecting outlier voltage but don't remove
                        for ch=1:EEG.nbchan
                            EEG = pop_eegthresh(EEG,1, ch, volt_threshold(1), volt_threshold(2), EEG.xmin, EEG.xmax,0,0);
                            EEG = eeg_checkset(EEG);
                            EEG = eeg_rejsuperpose(EEG, 1, 1, 1, 1, 1, 1, 1, 1);
                            badChans(ch,:) = EEG.reject.rejglobal;
                        end
                        tmpData = zeros(EEG.nbchan, EEG.pnts, EEG.trials);
                        for e = 1:EEG.trials
                            EEGe = []; badChanNum = []; % Initialize variables EEGe and badChanNum;
                            %select only this epoch (e)
                            EEGe = pop_selectevent( EEG, 'epoch',e,'deleteevents','off','deleteepochs','on','invertepochs','off');
                            badChanNum = find(badChans(:,e)==1); %find which channels are bad for this epoch
                            % find and add flat chans to the bad chans list
                            flatChanNum = find(range(EEGe.data,2) < 1);
                            badChanNum  = unique([badChanNum; flatChanNum]);
                            % find chans with a large jump/deflection (also bad chans) by taking 1st derivative
                            [jump_chans, ~] = find( abs(diff(EEGe.data,1,2) ./ repmat(diff(1:EEGe.pnts),EEGe.nbchan,1)) > 50);
                            badChanNum = unique([badChanNum; unique(jump_chans)]);
                            % add any new bad channels back to the bad chan list
                            badChans(badChanNum,e) = 1; % modify
                            % keep track of flat and jump channel information
                            %flat_mat(e) = ~isempty(flatChanNum);
                            %jump_mat(e) = length(unique(jump_chans));
                        end
                        badepoch=zeros(1, EEG.trials);
                        for ei=1:EEG.trials
                            if sum(badChans(:,ei)) > 0 % check if there are any bad chans
                                badepoch (ei)= 1;
                            end
                        end
                        badepoch=logical(badepoch);
                    end
                    % If all epochs are artifacted, save the dataset and ignore rest of the preprocessing for this subject.
                    if sum(EEG.reject.rejthresh)==EEG.trials || sum(EEG.reject.rejthresh)+1==EEG.trials || sum(badepoch)==EEG.trials || sum(badepoch)+1==EEG.trials
                        all_bad_epochs=1;
                        warning(['No usable data for datafile', datafile_names{subject}]);
                    else
                        if run_miniMADE == 0
                            EEG = pop_rejepoch(EEG,(EEG.reject.rejthresh), 0);
                            EEG = eeg_checkset(EEG);
                        elseif run_miniMADE == 1
                            EEG = pop_rejepoch(EEG, badepoch, 0);
                            EEG = eeg_checkset(EEG);
                        end
                    end
                end % end of epoch level channel interpolation if statement
            end % end of voltage threshold rejection if statement
            
        elseif allow_missing_chans == 1 % If advanced option to replace bad channels with NaNs is selected
            if rerefer_data==1
                % grab channels for rereferencing
                if iscell(reref)==1
                    reref_idx=zeros(1, length(reref));
                    for rr=1:length(reref)
                        reref_idx(rr)=find(strcmp({EEG.chanlocs.labels}, reref{rr}));
                    end
                    reref_chans = reref_idx;
                else
                    reref_chans = reref; 
                end
                % perform traditional artifact rejection for rereference channels
                EEG = pop_eegthresh(EEG, 1, reref_chans, volt_threshold(1), volt_threshold(2), EEG.xmin, EEG.xmax, 0, 0);
                EEG = eeg_checkset( EEG );
                EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
                
                if length(find(EEG.reject.rejthresh)) == EEG.trials || length(find(EEG.reject.rejthresh))+1 == EEG.trials % if all epochs marked for rejection
                    all_bad_epochs = 1; % set at 1 so we don't save the file below
                else
                    % reject epochs where rereference channels are bad
                    EEG = pop_rejepoch( EEG, (EEG.reject.rejthresh), 0);
                    EEG = eeg_checkset(EEG);
                    % re-reference to new reference channels
                    EEG = eeg_checkset(EEG);
                    EEG = pop_reref( EEG, reref_chans);
                end
            end
            
            % grab new channel labels in case some channels have been removed
            chans=[]; chans_labels2=[];
            chans_labels2=cell(1,EEG.nbchan);
            for i=1:EEG.nbchan
                chans_labels2{i}= EEG.chanlocs(i).labels;
            end
                
            if all_bad_epochs == 0 && blink_check == 1 % look at user entered frontal channels and remove artefacted epochs
                frontal_channels_idx=zeros(1, length(frontal_channels));
                for rr=1:length(frontal_channels)
                    frontal_channels_idx(rr)=find(strcmp({EEG.chanlocs.labels}, frontal_channels{rr}));
                end
                % perform traditional artifact rejection for rereference channels
                EEG = pop_eegthresh(EEG, 1, frontal_channels_idx, volt_threshold(1), volt_threshold(2), EEG.xmin, EEG.xmax, 0, 0);
                EEG = eeg_checkset( EEG );
                EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
                % only mark epochs for rejection if all frontal channels are bad
                blink_epochs = find(sum(EEG.reject.rejthreshE(frontal_channels_idx,:))==length(frontal_channels_idx));
                if length(blink_epochs) == EEG.trials || length(blink_epochs)+1 == EEG.trials  % if all epochs marked for rejection
                    all_bad_epochs = 1; % set at 1 so we don't save the file below
                else
                    % reject epochs where rereference channels are bad
                    EEG = pop_rejepoch( EEG, blink_epochs, 0);
                    EEG = eeg_checkset(EEG);
                end
            end
            
            if all_bad_epochs == 0
                % replace remaining bad channels with NaNs
                badChans = zeros(EEG.nbchan, EEG.trials);
                % Find artifacted epochs by detecting outlier voltage but don't remove
                for ch=1:EEG.nbchan
                    EEG = pop_eegthresh(EEG,1, ch, volt_threshold(1), volt_threshold(2), EEG.xmin, EEG.xmax,0,0);
                    EEG = eeg_checkset(EEG);
                    EEG = eeg_rejsuperpose(EEG, 1, 1, 1, 1, 1, 1, 1, 1);
                    badChans(ch,:) = EEG.reject.rejglobal;
                end
                tmpData = zeros(EEG.nbchan, EEG.pnts, EEG.trials);
                for e = 1:EEG.trials
                    EEGe = []; badChanNum = []; % Initialize variables EEGe and EEGe_interp;
                    %select only this epoch (e)
                    EEGe = pop_selectevent( EEG, 'epoch',e,'deleteevents','off','deleteepochs','on','invertepochs','off');
                    badChanNum = find(badChans(:,e)==1); %find which channels are bad for this epoch
                    % find and add flat chans to the bad chans list
                    flatChanNum = [find(range(EEGe.data(:,1:(EEG.pnts/2)),2) < 1); find(range(EEGe.data(:,(EEG.pnts/2):EEG.pnts),2) < 1)];
                    badChanNum  = unique([badChanNum; unique(flatChanNum)]);
                    % find chans with a large jump/deflection (also bad chans) by taking 1st derivative
                    [jump_chans, ~] = find( abs(diff(EEGe.data,1,2) ./ repmat(diff(1:EEGe.pnts),EEGe.nbchan,1)) > 50);
                    badChanNum = unique([badChanNum; unique(jump_chans)]);
                    % add any new bad channels back to the bad chan list
                    badChans(badChanNum,e) = 1; % modify badChan list
                    % remove the bad chans for this epoch (replace with NaN)
                    EEGe = eeg_checkset( EEGe );
                    EEGe.data(badChanNum,:) = NaN;
                    tmpData(:,:,e) = EEGe.data; % store NaN replaced data into matrix
                    % keep track of flat and jump channel information
                    %flat_mat(e) = ~isempty(flatChanNum);
                    %jump_mat(e) = length(unique(jump_chans));
                end
                EEG.data = tmpData; % save data containing NaNs back to EEG struct
                
                % reject epochs where NaN channel numbers are above the user entered threshold
                badepoch=zeros(1, EEG.trials);
                for ei=1:EEG.trials
                    if sum(badChans(:,ei)) > chan_thresh*EEG.nbchan % check if there are any bad chans
                        badepoch(ei)= 1;
                    end
                end
                badepoch=logical(badepoch);
                if sum(badepoch)==EEG.trials || sum(badepoch)+1==EEG.trials
                    all_bad_epochs=1;
                    warning(['No usable data for datafile', datafile_names{subject}]);
                else
                    EEG = pop_rejepoch(EEG, badepoch, 0);
                    EEG = eeg_checkset(EEG);
                    badChans = badChans(:,~badepoch);
                end
                    
                % reformat badChans for output table
                total_epochs_by_chan_after_artifact_rejection = {'';[]}; badChansSum = [];
                total_epochs_by_chan_after_artifact_rejection(1,1:EEG.nbchan) = strcat(chans_labels2,'_total_epochs_after_artifact_rejection');
                badChansSum = sum((badChans-1)*-1,2)'; 
                for cch = 1:EEG.nbchan
                    total_epochs_by_chan_after_artifact_rejection{2,cch} = badChansSum(cch);
                end
            else
                total_epochs_by_chan_after_artifact_rejection = {'';[]}; badChansSum = [];
                total_epochs_by_chan_after_artifact_rejection(1,1:EEG.nbchan) = strcat(chans_labels2,'_total_epochs_after_artifact_rejection');
                for cch = 1:EEG.nbchan
                    total_epochs_by_chan_after_artifact_rejection{2,cch} = 0;
                end
            end
            total_epochs_by_chan_after_artifact_rejection{1,EEG.nbchan+1} = 'datafile_names';
            total_epochs_by_chan_after_artifact_rejection{2,EEG.nbchan+1} =  datafile_names{subject};
        end % end advanced option to use NaNs
        
        % if all epochs are found bad during artifact rejection
        if all_bad_epochs==1
            total_epochs_after_artifact_rejection(subject)=0;
            total_channels_interpolated(subject)=0;
            if output_format==1
                EEG = eeg_checkset(EEG);
                EEG = pop_editset(EEG, 'setname',  strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_epochs'));
                EEG = pop_saveset(EEG, 'filename', strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_epochs.set'),'filepath', [output_location filesep 'processed_data' filesep ]); % save .set format
            elseif output_format==2
                save([[output_location filesep 'processed_data' filesep ] strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_epochs.mat')], 'EEG'); % save .mat format
            elseif output_format==3
                EEG = eeg_checkset(EEG);
                EEG = pop_editset(EEG, 'setname',  [current_subject, '_task-', task_name, '_run-01_eeg', '_no_usable_data_all_bad_epochs']);
                EEG = pop_saveset(EEG, 'filename', [current_subject, '_task-', task_name, '_run-01_eeg', '_no_usable_data_all_bad_epochs.set'],'filepath', [output_location_derivatives filesep 'eegpreprocess' filesep current_subject]); % save BIDS format
            end
            
            savedDataPath = strcat(output_location,'/processed_data/', strrep(datafile_names{subject}, ext, '_no_usable_data_all_bad_epochs.set'));

            cd(output_location)
                if exist('miniMADE_preprocessing_report.mat','file')
                    load miniMADE_preprocessing_report.mat
                end
                report_table_newrow=table(datafile_names(subject)', {curdate}, {savedDataPath}, total_epochs_before_artifact_rejection(subject)', total_epochs_after_artifact_rejection(subject)',...
                    LS_check, {N_TaskOnsetMarkers}, {AR1_badFrontCh}, {AR1_blinks}, Neps_postAR1, {AR2a_thresholds}, {AR2b_flat}, {AR2c_jumps}, {AR2_thr_flat_jump}, {Eps_BAD_InvalidInterp}, Neps_postAR2, {Chan_labels_curr});
                report_table_newrow.Properties.VariableNames={'datafile_names', 'CurDate',  'DataPath', 'total_epochs_before_artifact_rejection', 'total_epochs_after_artifact_rejection',...
                    'LS','EEGevents','AR1_badfrontCh','AR1_blinks','Neps_postAR1','AR2a_thres','AR2b_flat','AR2c_jump','AR2_summary','Eps_toomuchinterpolation','Neps_postAR2','Chan_labels'};
                miniMADE_report_table(subject,:) = report_table_newrow;
                save('miniMADE_preprocessing_report.mat','miniMADE_report_table')


            continue % ignore rest of the processing and go to next datafile
        else
            total_epochs_after_artifact_rejection(subject)=EEG.trials;
        end
        
        %% STEP 15: Interpolate deleted channels - not for low-density
        
        %% STEP 16: Rereference data
        if allow_missing_chans == 0
            if rerefer_data==1
                if iscell(reref)==1
                    reref_idx=zeros(1, length(reref));
                    for rr=1:length(reref)
                        reref_idx(rr)=find(strcmp({EEG.chanlocs.labels}, reref{rr}));
                    end
                    EEG = eeg_checkset(EEG);
                    EEG = pop_reref( EEG, reref_idx);
                else
                    EEG = eeg_checkset(EEG);
                    EEG = pop_reref(EEG, reref);
                end
            end
        end
        
        %% Save processed data
        if output_format==1
            EEG = eeg_checkset(EEG);
            EEG = pop_editset(EEG, 'setname',  strrep(datafile_names{subject}, ext, '_processed_data'));
            EEG = pop_saveset(EEG, 'filename', strrep(datafile_names{subject}, ext, '_processed_data.set'),'filepath', [output_location filesep 'processed_data' filesep ]); % save .set format
        elseif output_format==2
            save([[output_location filesep 'processed_data' filesep ] strrep(datafile_names{subject}, ext, '_processed_data.mat')], 'EEG'); % save .mat format
        elseif output_format==3
            EEG = eeg_checkset(EEG);
            EEG = pop_editset(EEG, 'setname',  [current_subject, '_task-', task_name, '_run-01', '_desc-processed_eeg']);
            EEG = pop_saveset(EEG, 'filename', [current_subject, '_task-', task_name, '_run-01', '_desc-processed_eeg.set'],'filepath', [output_location_derivatives filesep 'eegpreprocess' filesep current_subject]); % save BIDS format
        end

        savedDataPath = strcat(output_location,'/processed_data/',strrep(datafile_names{subject}, ext, '_processed_data.set'));

    catch
        if length(total_epochs_before_artifact_rejection) < subject || isempty(total_epochs_before_artifact_rejection(subject))
            total_epochs_before_artifact_rejection(subject) = NaN;
        end 
        if length(total_epochs_after_artifact_rejection) < subject || isempty(total_epochs_after_artifact_rejection(subject))
            total_epochs_after_artifact_rejection(subject) = NaN;
        end 
    end
    
    %% Create the report table for all the data files with relevant preprocessing outputs.
    cd(output_location)
    if exist('miniMADE_preprocessing_report.mat','file')
        load miniMADE_preprocessing_report.mat
    end
    report_table_newrow=table(datafile_names(subject)', {curdate}, {savedDataPath}, total_epochs_before_artifact_rejection(subject)', total_epochs_after_artifact_rejection(subject)',...
        LS_check, {N_TaskOnsetMarkers}, {AR1_badFrontCh}, {AR1_blinks}, Neps_postAR1, {AR2a_thresholds}, {AR2b_flat}, {AR2c_jumps}, {AR2_thr_flat_jump}, {Eps_BAD_InvalidInterp}, Neps_postAR2, {Chan_labels_curr});
    report_table_newrow.Properties.VariableNames={'datafile_names', 'CurDate',  'DataPath', 'total_epochs_before_artifact_rejection', 'total_epochs_after_artifact_rejection',...
        'LS','EEGevents','AR1_badfrontCh','AR1_blinks','Neps_postAR1','AR2a_thres','AR2b_flat','AR2c_jump','AR2_summary','Eps_toomuchinterpolation','Neps_postAR2','Chan_labels'};
    miniMADE_report_table(subject,:) = report_table_newrow;
    save('miniMADE_preprocessing_report.mat','miniMADE_report_table')

    % end
end % end of subject loop

disp('Done')