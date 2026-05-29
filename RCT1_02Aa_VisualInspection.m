%% RCT_02Aa: Visual inspection of QC power spectrum whole session pre pre-proc vs cleaned for task

% This script takes the raw data of the whole session and the pre-processed
% data and creates a reporting figure for visual inspection. 
% The following are visualised for reporting for visual inspection:
% - channel at which the light sensor has been identified
% - 50Hz line noise per channel (raw data - session quality)
% - variance of the amplitude differences (raw data - session quality)
% - correlation and variance across channels across the session (raw data - session quality)
% - power spectra for each channel (raw data - session quality)
% - number of EEG events in the data for each condition
% - whole scalp power (across all clean data) for each condition with the
% number of clean epochs per condition

% Rianne Haartsen, PhD.; 2024
% Birkbeck University of London

%%
% add paths with scripts, data, and output location
addpath('xxx/Rbaclofen_RCT1')

output_location = 'xxx/RCT1_EEG/DATA/01A_PreprocPower';

cd(output_location)
if exist([output_location filesep 'PowSpreport_figures'], 'dir') == 0
    mkdir([output_location filesep 'PowSpreport_figures'])
end
load xxx/RCT1_EEG/DATA/01A_PreprocPower/miniMADE_MeasuresOfInterest_report.mat
load xxx/RCT1_EEG/DATA/01A_PreprocPower/miniMADE_preprocessing_report.mat

ftrawdata_location = 'xxx';

    %add code from Luke Mason's TaskEngine2 and Rianne's code
    addpath(genpath('xxx/TaskEngine2'))
    addpath(genpath('xxx/lm_tools'));
    addpath(genpath('xxx/braintools_code_LM'));  
    addpath(genpath('xxx/eegtools'));
    addpath('xxx/BrT_Arb_scripts');
    addpath('xxx/braintools_code_RH');
    addpath('xxx/braintools_code_RH/Arbaclofen');
    % add fieldtrip path
    addpath('xxx/fieldtrip-20230118/')
    ft_defaults


for ss = 1:height(miniMADE_MOI_table)

        % read in data
        cd(ftrawdata_location)
        load(strcat('fieldtrip_',miniMADE_MOI_table.ID{ss},'.mat'))
        if exist('ft_data')
            RAWdata = ft_data;
            clear ft_data
        else
            RAWdata = ft;
            clear ft 
        end

        % create eeg data quality summary struct
            eegQC = struct();
            eegQC.QC_date = datestr(now,'dd-mm-yyyy'); 
        
        % QC with continous data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if strcmp(RAWdata.label{1},'Ch1') && size(RAWdata.label,2) == 20
                % rename labels
                Chan_names = {'P7','P4','Cz','Pz','P3','P8','Oz','O2','T8','PO8','C4','F4','AF8','Fz','C3','F3','AF7','T7','PO7','Fpz'};
                Nchs = size(RAWdata.label,2);
                RAWdata.label = [];
                for cc = 1:Nchs
                    RAWdata.label{cc,1} = Chan_names{1,cc};
                end
            end
            eegQC.ChanLabels = RAWdata.label;
            
        % Check for LS data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            [anyFound, idx_chan] = eegFT_findLightSensorChannel(RAWdata);
                [anyFound, idx_chan] = eegFT_findLightSensorChannel(RAWdata);
            if isequal(anyFound,1)
                disp(strcat('Light sensor channel identified:    ', RAWdata.label{idx_chan}))
                eegQC.LS_found = 1;
                eegQC.LS_idx_chan = idx_chan;
            else
                disp({'Light sensor channel not identified'; 'Using enobio markers to segment'})
                eegQC.LS_found = 0;
                eegQC.LS_idx_chan = [];
            end
    
    
        % Amplitude for 50Hz noise %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            eegQC.Noise50Hz = zeros(size(eegQC.ChanLabels,1),1); 
        % segment it into 1-second pieces with 50% overlap
            cfg = [];
            cfg.length      = 1;
            cfg.overlap     = .5;
            RawSegmented = ft_redefinetrial(cfg, RAWdata);
        % frequency analysis for 49 - 51 Hz
            cfg = [];
            cfg.method      = 'mtmfft';
            cfg.output      = 'pow';
            cfg.taper       = 'hanning';
            cfg.foilim      = [49 51];
            FreqDataS = ft_freqanalysis(cfg,RawSegmented);
            eegQC.Noise50Hz = log(mean(FreqDataS.powspctrm,2)); % calculate average natural log power across 49-51 Hz for each channel
            clear FreqDataS cfg 
        % power spectra 
            cfgpow                 = [];
            cfgpow.output          = 'pow';
            cfgpow.method          = 'mtmfft';
            cfgpow.taper           = 'hanning';
            cfgpow.foi             = 1:120;  
            cfgpow.keeptrials      = 'no';
            FreqData_Ch = ft_freqanalysis(cfgpow, RawSegmented);
            eegQC.Pow_perChannel = FreqData_Ch; % calculate average natural log power across 49-51 Hz for each channel
            clear FreqDataS cfg RawSegmented FreqData_Ch;
               
        % Channel variance and correlations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % segment into 1-sec epochs without overlap
            cfg = [];
            cfg.length      = 1;
            cfg.overlap     = 0;
            cfg.channel     = {'all','-O2'};
            RawSeconds = ft_redefinetrial(cfg, RAWdata);
            clear cfg
            
        % calculate differences in amplitudes between samples
            eegQC.CorrDiff = zeros((size(RawSeconds.label,1)*(size(RawSeconds.label,1)-1))/2,size(RawSeconds.trial,2));
            eegQC.VarChan = zeros(size(RawSeconds.label,1),size(RawSeconds.trial,2));
            
            for ii = 1:size(RawSeconds.trial,2) % for each 1-sec segment:
                
                % calculate differences in amplitude between subsequent samples
                % within a 1-sec segment
                    DiffAmpl = diff(RawSeconds.trial{1,ii},1,2);
                % calculate the correlation between electrodes
                    x = corr(DiffAmpl'); 
                    eegQC.CorrDiff(:,ii) = squareform(tril(x,-1));
                    clear x 
                % calculate the variance for the differences in amplitudes across
                % the tsegment for each channel
                    eegQC.VarChan(:,ii) = var(RawSeconds.trial{1,ii},0,2);
                    clear DiffAmpl
            end
            clear ii
            clear RawContinuous RawSeconds 


        % Create figure
        
            QCwhole_powspctrmcleaned = figure('position',[1 84 1512 782],'MenuBar','none','Visible','off');
            Titles_size = 12;
            AxLabel_size = 10;
            Colours = jet(20);
            % create string on LS data
            if eegQC.LS_found == 1
                LSstr = strcat('Light sensor: ', eegQC.ChanLabels(eegQC.LS_idx_chan));
            else
                LSstr = 'No LS signal found';
            end
            % plot the 50Hz noise for each channel
                subplot(3,3,1) 
                plot(1:1:size(eegQC.ChanLabels,1), eegQC.Noise50Hz, 'Color',[.8, .8, .8])
                hold on
                scatter(1:1:size(eegQC.ChanLabels,1),eegQC.Noise50Hz, 30, Colours, 'filled');
                ylabel('Log Power (Hz)', 'Fontsize', AxLabel_size);
                ylim([0 20])
                xlim([0 (size(eegQC.ChanLabels,1)+1)])
                xticks(1:1:size(eegQC.ChanLabels,1))
                xticklabels(eegQC.ChanLabels)
                xlabel('Channels', 'Fontsize', AxLabel_size);
                title({['50 Hz noise (mean: ' num2str(mean(eegQC.Noise50Hz(~isinf(eegQC.Noise50Hz)),1)) ' Hz)'], [LSstr{1,1}]}, 'Fontsize', Titles_size);
            % plot the channel variance for each channel
                subplot(3,3,4) 
                x = ones(1,size(eegQC.VarChan,2));
                scatter(x,eegQC.VarChan(1,:)', 20, Colours(1,:));
                hold on
                for ii = 2:size(eegQC.VarChan,1)
                    scatter(x*ii, eegQC.VarChan(ii,:)', 20, Colours(ii,:))
                end
                ylabel('Variance (signal)', 'Fontsize', AxLabel_size);
                set(gca, 'YScale', 'log')
                xlim([0 (size(eegQC.ChanLabels,1)+1)])
                xticks(1:1:size(eegQC.ChanLabels,1))
                xticklabels(eegQC.ChanLabels)
                xlabel('Channels', 'Fontsize', AxLabel_size);
                title({'Variance for 1-sec segments (ampl diff)'}, 'Fontsize', Titles_size);           
            % plot the channel correlations for differences in amplitudes between subsequent samples
                subplot(3,3,7)
                yyaxis left 
                plot(1:1:size(eegQC.CorrDiff,2),nanmean(eegQC.CorrDiff,1));
                ylabel('Mean', 'Fontsize', AxLabel_size); ylim([0 1]);
                hold on
            % plot the channel correlations for differences in amplitudes between subsequent samples
                yyaxis right
                plot(1:1:size(eegQC.CorrDiff,2),nanvar(eegQC.CorrDiff,0,1));
                ylabel('Variance', 'Fontsize', AxLabel_size); ylim([0 1]);
                xlim([0 (size(eegQC.CorrDiff,2))])
                xlabel('Time (sec)', 'Fontsize', AxLabel_size);
                title({'Correlations across all channels (ampl diff)'}, 'Fontsize', Titles_size);          
            % plot power spectra for all channels
                subplot(3,3,[2 3])
                plot(eegQC.Pow_perChannel.freq, log(eegQC.Pow_perChannel.powspctrm(1,:)),'Color', Colours(1,:))
                hold on
                for ii = 2:size(eegQC.VarChan,1)
                    plot(eegQC.Pow_perChannel.freq, log(eegQC.Pow_perChannel.powspctrm(ii,:)),'Color', Colours(ii,:))
                end
                ylabel('Log power'); xlabel('Frequency (Hz)')
                legend(eegQC.Pow_perChannel.label, 'Location','eastoutside')
                title({'Power spectra for each channel (raw data)'}, 'Fontsize', Titles_size);
            % plot N markers
                subplot(3,3,5)
                bar([1:6],miniMADE_report_table.EEGevents{ss}(2,:))
                xticklabels({'Soc En', 'Toy', 'Soc Sp', 'Soc Fr', 'Abs','Fix'})
                xtickangle(45)
                title({'EEG events for videos task'}, 'Fontsize', Titles_size)
                xlabel('Stimulus '); ylabel('N presentations')                
            % plot power spectrum for cleaned data
                subplot(3,3,[8 9])
                if ~isempty(miniMADE_MOI_table.Power_file{ss})
                    load(miniMADE_MOI_table.Power_file{ss})
                    PowerVals = zeros(4,120);
                    for cc = 1:length(DATA)
                        if DATA(cc).Ntrials > 1
                            PowerFreq = DATA(cc).freqs;
                            PowerVals(cc,:) = DATA(cc).Logpow;
                        end
                    end
                    plot(PowerFreq,PowerVals)
                    legend({strcat(DATA(1).Cond{1,1}, ' - N=', num2str(DATA(1).Ntrials)),...
                        strcat(DATA(2).Cond{1,1}, ' - N=', num2str(DATA(2).Ntrials)),...
                        strcat(DATA(3).Cond{1,1}, ' - N=', num2str(DATA(3).Ntrials)),...
                        strcat(DATA(4).Cond{1,1}, ' - N=', num2str(DATA(4).Ntrials))});
                    ylabel('Global log power'); xlabel('Frequency (Hz)')
                    title({'Power spectra across all channels (clean data)'},'Fontsize', Titles_size)
                else
                    text(0.5, 0.5, 'No power data available', ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'middle', ...
                        'Units', 'normalized', ...
                        'FontSize', 16, ...
                        'FontWeight', 'bold', ...
                        'Color', [0.5 0.5 0.5]);
                    title({'Power spectra across all channels (clean data)'},'Fontsize', Titles_size)
                end

            % general title
            IDparts = strsplit(miniMADE_MOI_table.ID{ss}, '_');
            sgtitle({strcat('Subject:',IDparts{1,1}); strcat('Session: ',IDparts{1,3}, ' Site: ', IDparts{1,2})}, 'Fontsize',16);

            % save the figure and clear up
            cd([output_location filesep 'PowSpreport_figures'])
            filename = sprintf(strcat(extractBefore(miniMADE_report_table.datafile_names{ss}, '_rawEEGlab.set'),'_PowSpreport.png'));
            saveas(QCwhole_powspctrmcleaned, filename);
            close(QCwhole_powspctrmcleaned)
            clear anyFound AxLabel_size cc cfgpow Chan_names Channel_labels Colours DATA eegQC ft filename 
            clear IDparts idx_chan ii Info LSstr PowerFreq PowerVals RAWdata Titles_size QCwhole_powspctrmcleaned
end
