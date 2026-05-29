%% RCT1_01: initial preprocessing and filtering of data

% This script does the first step of preprocessing by reading in the EEGlab
% formatted data. Then it identifies the channel with the light sensor
% signal, and applies filtering to the data (band pass filtering and line
% noise reduction.
% This step is consistent for all the different EEG metrics calculated in
% the dynamic and resting state videos in RCT1. 

% Created by dr. Rianne Haartsen, 2024

% Additions/ changes to miniMADE:
% - Adding in reduction for line noise (cleanLine from Tim Mullen)


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
addpath('xxx/eeglab2022.0_old/'); % use this path for most of the processing
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
output_location = 'xxx/DATA/01_Preproc_Filtered_Continuous';

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

% light sensor check
LS_check = [];

% channels interpolated
% total_channels_interpolated=[]; % total_channels_interpolated=faster_bad_channels+ica_preparation_bad_channels
curdate=datestr(now,'dd-mm-yyyy'); % set current date here so that table won't crash if date changes

%% Loop over all data files
for subject = 1:length(datafile_names)
    EEG=[];
    
    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', subject, datafile_names{subject});

        %% Initialize output variables

        % light sensor check
        LS_check = NaN;
        %
        % channels interpolated
        % total_channels_interpolated=[]; % total_channels_interpolated=faster_bad_channels+ica_preparation_bad_channels
        curdate=datestr(now,'dd-mm-yyyy'); % set current date here so that table won't crash if date changes

        savedDataPath = NaN; 

        cd(output_location)
        if exist('miniMADE_preprocessing_report.mat','file')
            load miniMADE_preprocessing_report.mat
        end

        report_table_newrow=table(datafile_names(subject)', {curdate}, {savedDataPath}, LS_check);
        report_table_newrow.Properties.VariableNames={'datafile_names', 'CurDate',  'DataPath', 'LS_check'};
        miniMADE_filtered_table(subject,:) = report_table_newrow;
        save('miniMADE_filtered_report.mat','miniMADE_filtered_table')

    try
        %% STEP 1: Import EGI data file and relevant information
        % EEG = mff_import([rawdata_location filesep datafile_names{subject}]);
        EEG = pop_loadset([rawdata_location filesep datafile_names{subject}]);
        EEG = eeg_checkset(EEG);

        % check if the channel labels are as expected
        if strcmp(EEG.urchanlocs(1).labels,'Ch1') && size(EEG.urchanlocs,2) == 20 ...
                && size(EEG.chanlocs,2) == 20
            % rename labels
            Chan_names = {'P7','P4','Cz','Pz','P3','P8','Oz','O2','T8','PO8','C4','F4','AF8','Fz','C3','F3','AF7','T7','PO7','Fpz'};
            for cc = 1:size(EEG.urchanlocs,2)
                EEG.urchanlocs(cc).labels = Chan_names{1,cc};
                EEG.chanlocs(cc).labels = Chan_names{1,cc};
            end
        end
    
        %% STEP 1.b: Check for light sensor activity - RCT1 specific

        % a) check for light sensor across all channels
        % set minimum number of contiguous samples above threshold that we
        % consider to be indicative of a light sensor turning on 
        min_time = 0.020;
        min_samps = min_time * EEG.srate;
        thresh = 1000;
        % loop through channels...
        numChan = size(EEG.data, 1);
        found = zeros(numChan, 1);
        for c = 1:numChan
            % threshold and find runs
            idx = EEG.data(c, :) >= thresh;
            ct = findcontig2(idx);
            % if not voltages above threshold, move to next channel
            if isempty(ct)
                continue
            end
            % remove runs below minimum duration
            idx_tooShort = ct(:, 3) < min_samps;
            ct(idx_tooShort, :) = [];
            % count number of runs left
            found(c) = length(ct);
        end
        
        % check at last one channel had some runs
        anyFound = ~all(found == 0);
        if anyFound
            % find channel with most runs
            idx_chan = found == max(found);
        else
            idx_chan = false(numChan, 1);
        end

        if anyFound ~= 0
            Suggested_LS = EEG.chanlocs(idx_chan).labels;
        else
            Suggested_LS = NaN;
        end

        clear min_time min_samps thresh numChan found c idx ct idx_tooShort anyFound idx_chan

        % b) check for activity in LS according to set-up = O2
        % select data with the lightsensor only
        EEG_lightsensor = pop_select(EEG, 'nochannel', find(~ismember(lower({EEG.chanlocs.labels}),lower('O2')))); 
        % check if there is a signal
        Diffs_between_samples = diff(EEG_lightsensor.data);
        if var(Diffs_between_samples,[],2) == 0
            warning('Flat signal for light sensor')
            LS_check = 0;
        elseif strcmp(Suggested_LS, 'O2') % suggested LS consistent with protocol
            LS_check = 1;
        elseif isnan(Suggested_LS) % no LS identified
            LS_check = 2;
        else 
            LS_check = 3;
        end

        clear EEG_lightsensor Diffs_between_samples Suggested_LS

        % Delete the LS channel 
        chans_labels=cell(1,EEG.nbchan);
        for i=1:EEG.nbchan
            chans_labels{i}= EEG.chanlocs(i).labels;
        end
        [chans,chansidx] = ismember(outerlayer_channel, chans_labels);
        outerlayer_channel_idx = chansidx(chansidx ~= 0);
        if delete_outerlayer==1
            if isempty(outerlayer_channel_idx)==1
                error(['None of the outer layer channels present in channel locations of data.'...
                    ' Make sure outer layer channels are present in channel labels of data (EEG.chanlocs.labels).']);
            else
                EEG = pop_select( EEG,'nochannel', outerlayer_channel_idx);
                EEG = eeg_checkset( EEG );
            end
        end

        clear chans_labels chans chansidx outlayer_channel_idx i

       %% STEP 1.5: Delete discontinuous data from the raw data file (OPTIONAL, but necessary for most EGI files)
        % NA here
        
        %% STEP 2: Import channel locations

        EEG_x.chanlocs = EEG.chanlocs;

        % RH: change the eeglab function path to read in the layout
        rmpath('xxx/eeglab2022.0_old/functions/sigprocfunc/');
        addpath('xxx/eeglab2024.0/functions/sigprocfunc');
        
        EEG = pop_chanedit(EEG, 'load',{channel_locations 'filetype' 'autodetect'});
        EEG = eeg_checkset( EEG );

        % RH: change the eeglab function path back
        rmpath('xxx/eeglab2024.0/functions/sigprocfunc');
        addpath('xxx/eeglab2022.0_old/functions/sigprocfunc/');

        % RH: correct for chan locs containing labels '1'-'20'
        for cc = 1:size(EEG.chanlocs,2)
            EEG.chanlocs(cc).labels = EEG_x.chanlocs(cc).labels;
        end
        clear cc EEG_x

        % Check whether the channel locations were properly imported. The EEG signals and channel numbers should be same.
        if size(EEG.data, 1) ~= length(EEG.chanlocs)
            error('The size of the data does not match with channel numbers.');
        end
        % Throw a warning if system is low density and user has selected to run the full MADE pipeline or interpolation
        if length(EEG.chanlocs) < 32 % if a low density system
            if run_miniMADE == 0
                warning('Running MADE is not recommended for low-density systems. To run miniMADE instead set run_miniMADE equal to 0 at the top of the script');
            end
            if interp_channels == 1 && length(EEG.chanlocs) < 20
                warning('Channel interpolation is not recommended for low-density systems with fewer than 20 channels');
            end
        end
        % if reref is a vector (type double) instead of a cell array, grab the electrode name(s)
        if ~iscell(reref); reref = {EEG.chanlocs(reref).labels}; end % this will prevent using the wrong channel in cases where we remove channels
        
        %% STEP 2.5: Label the task (OPTIONAL)
        % NA here - use own task segmenting scripts
        
        %% STEP 3: Adjust anti-aliasing and task related time offset
        % NA here

        %% STEP 4: Change sampling rate
        % NA here; sampling rate all the same due to same system in RCT1
        
        %% STEP 5: Delete outer layer of channels 
        % NA here; already taken out the LS signal
        
        %% STEP 6: Filter data
        % Calculate filter order using the formula: m = dF / (df / fs), where m = filter order,
        % df = transition band width, dF = normalized transition width, fs = sampling rate
        % dF is specific for the window type. Hamming window dF = 3.3
        
        high_transband = highpass; % high pass transition band
        low_transband = 10; % low pass transition band
        
        hp_fl_order = 3.3 / (high_transband / EEG.srate);
        lp_fl_order = 3.3 / (low_transband / EEG.srate);
        
        % Round filter order to next higher even integer. Filter order is always even integer.
        if mod(floor(hp_fl_order),2) == 0
            hp_fl_order=floor(hp_fl_order);
        elseif mod(floor(hp_fl_order),2) == 1
            hp_fl_order=floor(hp_fl_order)+1;
        end
        
        if mod(floor(lp_fl_order),2) == 0
            lp_fl_order=floor(lp_fl_order)+2;
        elseif mod(floor(lp_fl_order),2) == 1
            lp_fl_order=floor(lp_fl_order)+1;
        end
        
        % Calculate cutoff frequency
        high_cutoff = highpass/2;
        low_cutoff = lowpass + (low_transband/2);
        
        % Performing high pass filtering
        EEG = eeg_checkset( EEG );
        EEG = pop_firws(EEG, 'fcutoff', high_cutoff, 'ftype', 'highpass', 'wtype', 'hamming', 'forder', hp_fl_order, 'minphase', 0);
        EEG = eeg_checkset( EEG );
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        
        % pop_firws() - filter window type hamming ('wtype', 'hamming')
        % pop_firws() - applying zero-phase (non-causal) filter ('minphase', 0)
        
        % Performing low pass filtering
        EEG = eeg_checkset( EEG );
        EEG = pop_firws(EEG, 'fcutoff', low_cutoff, 'ftype', 'lowpass', 'wtype', 'hamming', 'forder', lp_fl_order, 'minphase', 0);
        EEG = eeg_checkset( EEG );
        
        % pop_firws() - transition band width: 10 Hz
        % pop_firws() - filter window type hamming ('wtype', 'hamming')
        % pop_firws() - applying zero-phase (non-causal) filter ('minphase', 0)
    
    % RH: try outs to remove line noise
        % remove Line noise and harmonics at 50 and 100Hz
        % From Github from Tim Mullen
        % add folder
        addpath(genpath('/Users/riannehaartsen/Documents/MATLAB/eeglab2024.0/plugins/cleanline-master'))
        % reduce line noise
        EEG = pop_cleanline(EEG, 'Bandwidth',4,'ChanCompIndices',[1:EEG.nbchan] ,...
            'SignalType','Channels','ComputeSpectralPower',true,'LineFrequencies',[50 100] ,...
            'NormalizeSpectrum',false,'LineAlpha',0.01,'PaddingFactor',2,'PlotFigures',false,...
            'ScanForLines',true,'SmoothingFactor',100,'VerbosityLevel',1,'SlidingWinLength',4,...
            'SlidingWinStep',2);
    
   %% Save the filtered data for later segmenting
     
        %% Save processed data
        if output_format==1
            EEG = eeg_checkset(EEG);
            EEG = pop_editset(EEG, 'setname',  strrep(datafile_names{subject}, ext, '_filtered_data'));
            EEG = pop_saveset(EEG, 'filename', strrep(datafile_names{subject}, ext, '_filtered_data.set'),'filepath', [output_location filesep 'processed_data' filesep ]); % save .set format
        elseif output_format==2
            save([[output_location filesep 'filtered_data' filesep ] strrep(datafile_names{subject}, ext, '_filtered_data.mat')], 'EEG'); % save .mat format
        elseif output_format==3
            EEG = eeg_checkset(EEG);
            EEG = pop_editset(EEG, 'setname',  [current_subject, '_task-', task_name, '_run-01', '_desc-filtered_eeg']);
            EEG = pop_saveset(EEG, 'filename', [current_subject, '_task-', task_name, '_run-01', '_desc-filteredd_eeg.set'],'filepath', [output_location_derivatives filesep 'eegpreprocess' filesep current_subject]); % save BIDS format
        end

        savedDataPath = strcat(output_location,'/filtered_data/',strrep(datafile_names{subject}, ext, '_filtered_data.set'));

    end % end of try
    
    %% Create the report table for all the data files with relevant preprocessing outputs.
    cd(output_location)
    if exist('miniMADE_filtered_report.mat','file')
        load miniMADE_filtered_report.mat
    end
    report_table_newrow=table(datafile_names(subject)', {curdate}, {savedDataPath}, LS_check);
    report_table_newrow.Properties.VariableNames={'datafile_names', 'CurDate',  'DataPath', 'LS_check'};
    miniMADE_filtered_table(subject,:) = report_table_newrow;
    save('miniMADE_filtered_report.mat','miniMADE_filtered_table')

    % end
end % end of subject loop

disp('Done')