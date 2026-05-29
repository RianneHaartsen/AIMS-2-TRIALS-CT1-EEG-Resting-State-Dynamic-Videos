%% Rbaclofen RCT1: FC values from miniMADE pipeline data 

% This script takes the fieldtrip structure. Data are split
% into conditions, then an FFT is calculated
% for 1 to 32 Hz for each trial. The dbWPLI is calculated across the frequencies. 
% Finally, dbWPLI values are averaged across frequencies for the theta and alpha band. 
% For theta: 6-7Hz
% For alpha: 8-12Hz
% dbWPLI across all channels is calculated and plotted into a figure for visual
% inspection.
% Processing is looped through participants with paths from the
% miniMADE_preprocessing_report.mat file from the previous script
% (RCT1_01C). Connectivity data are saved and summaries are reported in
% miniMADE_MeasuresOfInterest_report matlab file. 

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

InfoPath = 'xxx/DATA/01C_PreprocFC';

% create new folders
if exist(strcat(InfoPath,'/FC_data'),'dir') == 0
    mkdir(strcat(InfoPath,'/FC_data'))
end
if exist(strcat(InfoPath,'/FC_data/FC_figures'),'dir') == 0
    mkdir(strcat(InfoPath,'/FC_data/FC_figures'))
end

FC_location = strcat(InfoPath,'/FC_data');
output_location = InfoPath;

%%
% load miniMADE report
cd(InfoPath)
load miniMADE_preprocessing_report.mat

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
                clear EEG
    
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
            clear Ntrls_sc ft_EEG_sc cfg

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
            clear Ntrls_sc ft_EEG_sc cfg
            

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
            clear Ntrls_sc ft_EEG_sc cfg

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
            clear Ntrls_sc ft_EEG_sc cfg

            %% For each condition
            % Parameters across FFT analyses for FC
            % for power analysis 
                cfgfft             = [];
                cfgfft.method      = 'mtmfft';
                cfgfft.taper       =  'hanning'; 
                cfgfft.output      = 'fourier'; 
                cfgfft.tapsmofrq   = 1;
                cfgfft.foi         = 1:1:32; 
                cfgfft.keeptrials  = 'yes'; 
                % frequencies (in points) to examine
                F = 1:length(cfgfft.foi);

            for cc = 1:length(DATA)
                if DATA(cc).Ntrials >= 90
                    % Power analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    SpectralPower = ft_freqanalysis(cfgfft, DATA(cc).ft_EEG_timeseries);
                    % dbWPLI analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % find segments with overlap
                    if size(DATA(cc).ft_EEG_timeseries.trialinfo,1) > 1
                        OVERLAP = zeros(size(DATA(cc).ft_EEG_timeseries.trialinfo,1), size(DATA(cc).ft_EEG_timeseries.trialinfo,1));
                        for i=1: size(DATA(cc).ft_EEG_timeseries.trialinfo,1)-1
                            if DATA(cc).ft_EEG_timeseries.sampleinfo(i,2) > DATA(cc).ft_EEG_timeseries.sampleinfo((i+1),1) % if overlaps
                               OVERLAP (i, i+1)=1;  OVERLAP (i+1, i)=1; 
                            end
                        end
                        clear i
                    end
                    OVERLAP(find(eye(size(DATA(cc).ft_EEG_timeseries.trialinfo,1))==1)) = 1;
                    % dbWPLI
                    [PLI, WPLI, ubPLI, dbWPLI] = PLIbasic_and_directional2(SpectralPower.fourierspctrm, size(SpectralPower.fourierspctrm,1), F, OVERLAP); %%dPLIbasic
                    FC = struct();
                    FC.PLI = PLI;
                    FC.WPLI = WPLI;
                    FC.ubPLI = ubPLI;
                    FC.dbWPLI = dbWPLI;
                    FC.Freqs_Hz = SpectralPower.freq(F);
                    FC.Chan_label = SpectralPower.label;
                    FC.overlap = OVERLAP;
                    clear PLI WPLI ubPLI dbWPLI OVERLAP
                    % global dbWPLI
                    GlobdbWPLI = zeros(size(FC.dbWPLI,1),1);
                    for ff = 1:size(FC.dbWPLI,1)
                        FCmat = squeeze(FC.dbWPLI(ff,:,:));
                        Mat1 = tril(FCmat,-1);
                        Mat1(Mat1 == 0) = NaN;
                        GlobdbWPLI(ff,1) = mean(Mat1, 'all', 'omitnan'); 
                        clear Mat1 FCmat
                    end
                    FC.dbWPLI_glob = GlobdbWPLI';  
                    clear ff
                else % no data for that condition
                        SpectralPower = [];
                        FC = [];
                        FC.Freqs_Hz = [];
                        GlobdbWPLI = [];
                end
                % add data to DATA struct
                DATA_fc(cc).Cond = DATA(cc).Cond;
                DATA_fc(cc).Ntrials = DATA(cc).Ntrials;
                DATA_fc(cc).ft_EEG_FFTdata = SpectralPower;
                DATA_fc(cc).FCdata = FC;
                DATA_fc(cc).FC_freqs = FC.Freqs_Hz;
                DATA_fc(cc).FC_globdbWPLI = GlobdbWPLI;
                clear SpectralPower FC GlobdbWPLI
            end
            clear cc


             %% Create figure of FC spectra and Ns
             if any([DATA_fc.Ntrials] >= 90) 
                FunConCheck = figure('Visible','off');
                FunConVals = zeros(4,32);
                PowerFreq = 1:1:32;
                for cc = 1:length(DATA_fc)
                    if DATA_fc(cc).Ntrials >= 90
                        FunConVals(cc,:) = DATA_fc(cc).FC_globdbWPLI;
                    end
                end
                plot(PowerFreq,FunConVals)
                legend({strcat(DATA_fc(1).Cond{1,1}, ' - N=', num2str(DATA_fc(1).Ntrials)),...
                    strcat(DATA_fc(2).Cond{1,1}, ' - N=', num2str(DATA_fc(2).Ntrials)),...
                    strcat(DATA_fc(3).Cond{1,1}, ' - N=', num2str(DATA_fc(3).Ntrials)),...
                    strcat(DATA_fc(4).Cond{1,1}, ' - N=', num2str(DATA_fc(4).Ntrials))});
                ylabel('Global FC (dbWPLI)'); xlabel('Frequency (Hz)')
                title(CurrID)
    
                % save the figure
                cd([InfoPath filesep 'FC_data/FC_figures'])
                filename = sprintf(strcat(extractBefore(miniMADE_report_table.datafile_names{subject}, '_rawEEGlab.set'),'_FCSpectrum.png'));
                saveas(FunConCheck, filename);
            
                close(FunConCheck)
                clear cc FunConVals PowerFreq FunConCheck
             end
            
        catch
            warning('Unable to calculate connectivity measures of interest')
        end % of try line 62

        
        %% save data and add to current subject to report table 
        % keep track of saved data location and subject preprocessed
        % FC
        FC_file{subject} = strcat(FC_location, filesep, CurrID, '_FunCon_data.mat');
        save(FC_file{subject}, 'DATA_fc', 'DATA')

    else % if there are 0 trials
        FC_file{subject} = [];
        DATA_fc(1).Cond = {'Soc'}; DATA_fc(1).Ntrials = 0; DATA_fc(1).ft_EEG_FFTdata = []; DATA_fc(1).FCdata = []; DATA_fc(1).FC_freqs = []; DATA_fc(1).FC_globdbWPLI = [];
        DATA_fc(2).Cond = {'Toy'}; DATA_fc(2).Ntrials = 0; DATA_fc(2).ft_EEG_FFTdata = []; DATA_fc(2).FCdata = []; DATA_fc(2).FC_freqs = []; DATA_fc(2).FC_globdbWPLI = [];
        DATA_fc(3).Cond = {'Abs'}; DATA_fc(3).Ntrials = 0; DATA_fc(3).ft_EEG_FFTdata = []; DATA_fc(3).FCdata = []; DATA_fc(3).FC_freqs = []; DATA_fc(3).FC_globdbWPLI = [];
        DATA_fc(4).Cond = {'Fix'}; DATA_fc(4).Ntrials = 0; DATA_fc(4).ft_EEG_FFTdata = []; DATA_fc(4).FCdata = []; DATA_fc(4).FC_freqs = []; DATA_fc(4).FC_globdbWPLI = [];
    end % end of if file exists

        % add to tracking table
        cd(output_location)
        if exist(strcat(InfoPath, '/miniMADE_MeasuresOfInterestFC_report.mat'),'file')
            load miniMADE_MeasuresOfInterestFC_report.mat
        end
        data_table_newrow = table({CurrID}, DateCur(1,subject), FC_file(1,subject), ...
            DATA_fc(1).Cond, DATA_fc(1).Ntrials, {DATA_fc(1).FC_freqs}, {DATA_fc(1).FC_globdbWPLI}, ...
            DATA_fc(2).Cond, DATA_fc(2).Ntrials, {DATA_fc(2).FC_freqs}, {DATA_fc(2).FC_globdbWPLI}, ...
            DATA_fc(3).Cond, DATA_fc(3).Ntrials, {DATA_fc(3).FC_freqs}, {DATA_fc(3).FC_globdbWPLI}, ...
            DATA_fc(4).Cond, DATA_fc(4).Ntrials, {DATA_fc(4).FC_freqs}, {DATA_fc(4).FC_globdbWPLI});
        data_table_newrow.Properties.VariableNames={'ID','Date','FC_file',  ...
            'SCond','Seps' ,'Sfreqs', 'SGlobdbWPLI',...
            'TCond','Teps' ,'Tfreqs', 'TGlobdbWPLI',...
            'ACond','Aeps' ,'Afreqs', 'AGlobdbWPLI',...
            'FCond','Feps' ,'Ffreqs', 'FGlobdbWPLI'};
        miniMADE_MOIfc_table(subject,:) = data_table_newrow;
        save('miniMADE_MeasuresOfInterestFC_report.mat','miniMADE_MOIfc_table');
        clear miniMADE_MOIfc_table

    % clear up variables
    clear DATA_fc ft_EEG DATA 
    clear cfgfft F
    clear FC_file CurrID data_table_newrow DateCur

end

disp('Done')

