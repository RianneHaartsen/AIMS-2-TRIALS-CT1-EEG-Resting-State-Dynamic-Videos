%% RCT1_02B: version B) Power values from miniMADE pipeline data and preparation for FOOOF

% This first section of the script takes the data processed with script RCT1_01A and converges
% it to a fieldtrip structure. Epochs are split into conditions. 
% Then power is calculated for 1 to 120 Hz for each epoch (log power and 
% absolute power) and averagered across epochs within conditions. Power
% across all channels is calculated and plotted into a figure for visual
% inspection.
% Processing is looped through participants with paths from the
% miniMADE_preprocessing_report.mat file from the previous script
% (RCT1_01B). Power data are saved and summaries are reported in
% miniMADE_MeasuresOfInterest_report matlab file. 

% The second section of the script averages the absolute power per region
% and saves it as a .csv file for further FOOOF analyses in Python. 

% Rianne Haartsen, PhD.; 2024
% Birkbeck University of London

%%

clear % clear matlab workspace
clc % clear matlab command window
addpath('xxx/eeglab2022.0');% enter the path of the EEGLAB folder in this line
eeglab nogui;% call eeglab to set up the plugin
addpath('xxx/fieldtrip-20230118/')
ft_defaults

% add path to RCT1 scripts
addpath('xxx/Rbaclofen_RCT1')

InfoPath = 'xxx/DATA/01B_PreprocOOF';

% create new folders
if exist(strcat(InfoPath,'/A_Power_data'),'dir') == 0
    mkdir(strcat(InfoPath,'/A_Power_data'))
end
if exist(strcat(InfoPath,'/FOOOF_powerspectra'),'dir') == 0
    mkdir(strcat(InfoPath,'/FOOOF_powerspectra'))
end

Pow_location = strcat(InfoPath,'/A_Power_data');
output_location = InfoPath;

output4FOOOF = strcat(InfoPath,'/FOOOF_powerspectra');

%%
% load miniMADE report
load xxx/DATA/01B_PreprocOOF/miniMADE_preprocessing_report.mat

for subject = 1:height(miniMADE_report_table)
    
    CurrFileName = miniMADE_report_table.datafile_names{subject};
    CurrID = extractBefore(CurrFileName, '_rawEEGlab.set');
    DateCur{subject} = datestr(now,'dd-mm-yyyy');
    
    fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', subject, CurrID);
    
    if miniMADE_report_table.total_epochs_after_artifact_rejection(subject) > 0

        try
            % load preprocessed data
            EEG = pop_loadset(miniMADE_report_table.DataPath{subject});
            EEG = eeg_checkset(EEG);
    
            %% Convert preprocessed eeglab data to fieldtrip data
            ft_EEG = eeglab2fieldtrip(EEG,'raw');
    
            % create trialinfo and sampleinfo fields in accordance with fieldtrip
                ft_EEG.trialinfo_eeglab = ft_EEG.trialinfo;
                % go through epochs
                trialinfo_cur = zeros(size(ft_EEG.trialinfo_eeglab,1),1);
                sampleinfo_cur = zeros(size(ft_EEG.trialinfo_eeglab,1),2);
                eps_start_urevent = zeros(size(ft_EEG.trialinfo_eeglab,1),1);
                for eps = 1:size(ft_EEG.trialinfo_eeglab,1)
                    trialinfo_cur(eps) = ft_EEG.trialinfo_eeglab.type(eps);
                    % check for consistency of event marker
                    eps_start_urevent(eps) = EEG.epoch(eps).eventtype{1};
                    if ~isequal(trialinfo_cur(eps),EEG.urevent(EEG.epoch(eps).eventurevent{1}).type)
                        error('Inconsistency ft and eeglab event types')
                    end
                    % get sampleinfo
                    begsample = EEG.urevent(EEG.epoch(eps).eventurevent{1}).latency;
                    endsample = begsample + size(ft_EEG.time{eps},2)-1;
                    sampleinfo_cur(eps,:) = [begsample endsample];
                    clear ind_event begsample endsample
                end
                ft_EEG.trialinfo = trialinfo_cur;
                ft_EEG.sampleinfo = sampleinfo_cur;
                clear trialinfo_cur sampleinfo_cur eps eps_start_urevent
    
            % create time in accordance with fieldtrip
                ft_EEG.time_eeglab = ft_EEG.time;
                time = cell(1,size(ft_EEG.time_eeglab,2));
                % time vector
                time_end = ft_EEG.sampleinfo(end,2)*(1/ft_EEG.fsample);
                time_full = 0:(1/ft_EEG.fsample):time_end;
                for eps = 1:size(ft_EEG.time_eeglab,2)
                    time_cur = time_full(ft_EEG.sampleinfo(eps,1):ft_EEG.sampleinfo(eps,2));
                    time{1,eps} = time_cur;
                    clear time_cur
                end
                ft_EEG.time = time;
                clear time eps time_end time_full
    
            %% Split into conditions
            % markers: '910','911','917', '918', '901', '903'
            
            % social video
            if sum(ft_EEG.trialinfo == 910) > 0 % english
                cfg                 = [];
                cfg.trials          = ft_EEG.trialinfo == 910;
                ft_EEG_sc = ft_selectdata(cfg, ft_EEG);
                Ntrls_sc = sum(ft_EEG.trialinfo == 910);
            elseif sum(ft_EEG.trialinfo == 917) > 0 % spanish
                cfg                 = [];
                cfg.trials          = ft_EEG.trialinfo == 917;
                ft_EEG_sc = ft_selectdata(cfg, ft_EEG);
                Ntrls_sc = sum(ft_EEG.trialinfo == 917);
            elseif sum(ft_EEG.trialinfo == 918) > 0 % french
                cfg                 = [];
                cfg.trials          = ft_EEG.trialinfo == 918;
                ft_EEG_sc = ft_selectdata(cfg, ft_EEG);
                Ntrls_sc = sum(ft_EEG.trialinfo == 918);
            else 
                ft_EEG_sc = []; Ntrls_sc = 0;
            end
            DATA(1).Cond = {'Soc'};
            DATA(1).Ntrials = Ntrls_sc;
            DATA(1).ft_EEG_timeseries = ft_EEG_sc;
            clear Ntrls_sc ft_EEG_sc

            % toy video
            if sum(ft_EEG.trialinfo == 911) > 0
                cfg                 = [];
                cfg.trials          = ft_EEG.trialinfo == 911;
                ft_EEG_sc = ft_selectdata(cfg, ft_EEG);
                Ntrls_sc = sum(ft_EEG.trialinfo == 911);
            else
                ft_EEG_sc = []; Ntrls_sc = 0;
            end
            DATA(2).Cond = {'Toy'};
            DATA(2).Ntrials = Ntrls_sc;
            DATA(2).ft_EEG_timeseries = ft_EEG_sc;
            clear Ntrls_sc ft_EEG_sc
            

            % abstract video
            if sum(ft_EEG.trialinfo == 901) > 0
                cfg                 = [];
                cfg.trials          = ft_EEG.trialinfo == 901;
                ft_EEG_sc = ft_selectdata(cfg, ft_EEG);
                Ntrls_sc = sum(ft_EEG.trialinfo == 901);
            else
                ft_EEG_sc = []; Ntrls_sc = 0;
            end
            DATA(3).Cond = {'Abs'};
            DATA(3).Ntrials = Ntrls_sc;
            DATA(3).ft_EEG_timeseries = ft_EEG_sc;
            clear Ntrls_sc ft_EEG_sc

            % fixation
            if sum(ft_EEG.trialinfo == 903) > 0
                cfg                 = [];
                cfg.trials          = ft_EEG.trialinfo == 903;
                ft_EEG_sc = ft_selectdata(cfg, ft_EEG);
                Ntrls_sc = sum(ft_EEG.trialinfo == 903);
            else
                ft_EEG_sc = []; Ntrls_sc = 0;
            end
            DATA(4).Cond = {'Fix'};
            DATA(4).Ntrials = Ntrls_sc;
            DATA(4).ft_EEG_timeseries = ft_EEG_sc;
            clear Ntrls_sc ft_EEG_sc
            
            %% Parameters across power analyses
            % for power analysis 
                cfgpow                 = [];
                cfgpow.output          = 'pow';
                cfgpow.method          = 'mtmfft';
                cfgpow.taper           = 'hanning';
                cfgpow.foi             = 0:0.5:120;  
                cfgpow.keeptrials      = 'yes';

            %% For each condition
            for cc = 1:length(DATA)
                if DATA(cc).Ntrials > 1
                    % Power analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        SpectralPower = ft_freqanalysis(cfgpow, DATA(cc).ft_EEG_timeseries);
                        % log power
                        SpectralPower.logpowspctrm = log(SpectralPower.powspctrm);
                        SpectralPower.logpowspctrm_avg = squeeze(mean(SpectralPower.logpowspctrm,1));
                        SpectralPower.logpowspctrm_avg_global = squeeze(mean(SpectralPower.logpowspctrm_avg,1));
                        % abs power for FOOOF
                        PowAbs_avg = squeeze(mean(SpectralPower.powspctrm,1));
                        PowAbs_avg_glob = squeeze(mean(PowAbs_avg,1));
                        % other vars
                        freqs_var = SpectralPower.freq;
                        logPow = SpectralPower.logpowspctrm_avg_global;
                else % no data for that condition
                        SpectralPower = [];
                        freqs_var = [];
                        logPow = [];
                        PowAbs_avg_glob = [];
                end
                % add data to DATA struct
                DATA(cc).ft_EEG_freqdata = SpectralPower;
                DATA(cc).freqs = freqs_var;
                DATA(cc).Logpow = logPow;
                DATA(cc).Abspow = PowAbs_avg_glob;
                clear SpectralPower freqs_var logPow PowAbs_avg_glob PowAbs_avg
            end
            clear cc

            %% Save absolute power for FOOOF if available
            % note: this saves the power across all channels
            for ss = 1:length(DATA)
                if DATA(ss).Ntrials > 20
                    % reformat and save
                    cd(output4FOOOF)
                    Tab = table(DATA(ss).Abspow');
                    writetable(Tab, strcat(CurrID,'_AbsPow_', char(DATA(ss).Cond) ,'_allchs.csv'))
                    clear Tab
                end
            end
            clear ss
            
        catch
            warning('Unable to calculate measures of interest')
        end % of try line 62

        
        %% save data and add to current subject to report table 
        % keep track of saved data location and subject preprocessed
        % fieldtrip data - all trials
        Ft_file{subject} = strcat(output_location, filesep, 'processed_data', filesep, CurrID, '_fieldtrip_data.mat');
        save(Ft_file{subject},'ft_EEG')
        % Power
        Power_file{subject} = strcat(Pow_location, filesep, CurrID, '_Power_data.mat');
        save(Power_file{subject}, 'DATA')

    else % if there are 0 trials
        Ft_file{subject} = [];
        Power_file{subject} = [];
        DATA(1).Cond = {'Soc'}; DATA(1).Ntrials = 0; DATA(1).freqs = []; DATA(1).Logpow = []; DATA(1).Abspow = [];
        DATA(2).Cond = {'Toy'}; DATA(2).Ntrials = 0; DATA(2).freqs = []; DATA(2).Logpow = []; DATA(2).Abspow = [];
        DATA(3).Cond = {'Abs'}; DATA(3).Ntrials = 0; DATA(3).freqs = []; DATA(3).Logpow = []; DATA(3).Abspow = [];
        DATA(4).Cond = {'Fix'}; DATA(4).Ntrials = 0; DATA(4).freqs = []; DATA(4).Logpow = []; DATA(4).Abspow = [];
    end % end of if file exists

        % add to tracking table
        cd(output_location)
        if exist(strcat(InfoPath, '/miniMADE_MeasuresOfInterest_report.mat'),'file')
            load miniMADE_MeasuresOfInterest_report.mat
        end
        data_table_newrow = table({CurrID}, DateCur(1,subject), Ft_file(1,subject), Power_file(1,subject), ...
            DATA(1).Cond, DATA(1).Ntrials, {DATA(1).freqs}, {DATA(1).Logpow}, {DATA(1).Abspow}, ...
            DATA(2).Cond, DATA(2).Ntrials, {DATA(2).freqs}, {DATA(2).Logpow}, {DATA(2).Abspow}, ...
            DATA(3).Cond, DATA(3).Ntrials, {DATA(3).freqs}, {DATA(3).Logpow}, {DATA(3).Abspow}, ...
            DATA(4).Cond, DATA(4).Ntrials, {DATA(4).freqs}, {DATA(4).Logpow}, {DATA(4).Abspow});
        data_table_newrow.Properties.VariableNames={'ID','Date','Fieldtrip_file','Power_file',  ...
            'SCond','Seps' ,'Sfreqs', 'SLogPow','SAbsPow', ...
            'TCond','Teps' ,'Tfreqs', 'TLogPow','TAbsPow', ...
            'ACond','Aeps' ,'Afreqs', 'ALogPow','AAbsPow', ...
            'FCond','Feps' ,'Ffreqs', 'FLogPow','FAbsPow'};
        miniMADE_MOI_table(subject,:) = data_table_newrow;
        save('miniMADE_MeasuresOfInterest_report.mat','miniMADE_MOI_table');
        clear miniMADE_MOI_table

    % clear up variables
    clear ft_EEG DATA EEG
    clear SpectralPower 
    clear cfgpow cfg
    clear Ft_file Power_file CurrID CurrFileName data_table_newrow DateCur

end



%% Extract spower spectra per region for FOOOF

clear % clear matlab workspace
clc % clear matlab command window

% add path to RCT1 scripts
addpath('xxx/Rbaclofen_RCT1')
InfoPath = 'xxx/DATA/01B_PreprocOOF';

% create new folders
if exist(strcat(InfoPath,'/A_Power_data'),'dir') == 0
    mkdir(strcat(InfoPath,'/A_Power_data'))
end
if exist(strcat(InfoPath,'/FOOOF_powerspectra_Region'),'dir') == 0
    mkdir(strcat(InfoPath,'/FOOOF_powerspectra_Region'))
end

Pow_location = strcat(InfoPath,'/A_Power_data');
output_location = InfoPath;
output4FOOOF = strcat(InfoPath,'/FOOOF_powerspectra_Region');

% load miniMADE measures of interest report
cd(output_location)
load miniMADE_MeasuresOfInterest_report.mat

for subject = 1:height(miniMADE_MOI_table)

    % info
    CurrID = miniMADE_MOI_table.ID{subject};

    % Load power data if available
    if ~isempty(miniMADE_MOI_table.Power_file{subject})
        load(miniMADE_MOI_table.Power_file{subject})

            % Loop through conditions
            for Cond = 1:4
                if DATA(Cond).Ntrials >= 20
                    % get absolute powspctrm
                    AbsPow = DATA(Cond).ft_EEG_freqdata;
                    % average across trials
                    AbsPow_mntrls = shiftdim(mean(AbsPow.powspctrm,1)); 
        
                    % average across regions
                    Ch_labels = DATA(Cond).ft_EEG_freqdata.label;
                    Chind_front = find(ismember(Ch_labels,{'Fz','F3','F4'}));
                    Chind_cent = find(ismember(Ch_labels,{'Cz','C3','C4'}));
                    Chind_pari = find(ismember(Ch_labels,{'Pz','P3','P4'}));
                    Chind_occip = find(ismember(Ch_labels,{'PO7','Oz','PO8'}));
        
                    AbsPow_Fr = mean(AbsPow_mntrls(Chind_front,:),1);
                    AbsPow_Ce = mean(AbsPow_mntrls(Chind_cent,:),1);
                    AbsPow_Pa = mean(AbsPow_mntrls(Chind_pari,:),1);
                    AbsPow_Oc = mean(AbsPow_mntrls(Chind_occip,:),1);
        
                    % save data into csv files
                    cd(output4FOOOF)
                    writetable(array2table(AbsPow_Fr), strcat(CurrID,'_', DATA(Cond).Cond{1,1},'_Frontal.csv'))
                    writetable(array2table(AbsPow_Ce), strcat(CurrID,'_', DATA(Cond).Cond{1,1},'_Central.csv'))
                    writetable(array2table(AbsPow_Pa), strcat(CurrID,'_', DATA(Cond).Cond{1,1},'_Parietal.csv'))
                    writetable(array2table(AbsPow_Oc), strcat(CurrID,'_', DATA(Cond).Cond{1,1},'_Occipital.csv'))
                    
                    % clean up
                    clear AbsPow AbsPow_mntrls Ch_labels Chind_front Chind_cent Chind_pari Chind_occip
                    clear AbsPow_Fr AbsPow_Ce AbsPow_Pa AbsPow_Oc
                
                end
            end
        
            % clean up 
            clear DATA Cond CurrID

    end

   
end
