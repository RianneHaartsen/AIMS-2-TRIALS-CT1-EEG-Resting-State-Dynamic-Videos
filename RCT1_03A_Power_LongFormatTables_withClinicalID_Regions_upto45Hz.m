%% RCT_03A: Spectral power for statistical analyses

% This script takes the tables with Metrics of Interest from power (script
% RCT1_02A) and assembles them into long format for statistical analysis in 
% R studio. 
% The script averages power data within region rather than across the whole
% scalp. Log power, absolute power, and relative power are extracted.
% Here we are taking frequencies up to 45Hz instead of 95Hz to include the
% 28 dataset with bump artefacts in the frequencies above 50Hz.

% The first section of the script extracts the data for the low and high
% frequencies (2-12Hz and 15-45Hz). The second section of the script
% extracts the data for the canonical frequency bands for follow up
% analyses of the wider frequency ranges. 
% Data are saved in a .csv file for further statistical analyses in
% Rstudio.

% Rianne Haartsen, PhD.; Aug - 2025
% Birkbeck University of London

%% Spectral power: a) absolute, log, and relative power %%%%%%%%%%%%%%%%%%%

MOI_location = 'xxx/DATA/01A_PreprocPower';
cd(MOI_location)
load miniMADE_MeasuresOfInterest_report.mat

PSCcode_conversion = readtable('xxx/rct1_psc1_psc2.xlsx');


% Columns: IDfull, ID, session, condition (4), region (4), frequency band (2), N trials, absolute power, log power, relative power
% frequency bands: 
% low: 2-12, high: 15-45 (based on feedback Paul Wang & Sarah
% Lippé during meeting in Dec 2024)
% 32 rows per dataset

% Categories:
Conditions_N = 4;
Regions_N = 4;
Freqs_N = 2;

Condition = [repmat({'Social'}, Regions_N*Freqs_N, 1); ...
    repmat({'Toy'}, Regions_N*Freqs_N, 1);...
    repmat({'Abstract'}, Regions_N*Freqs_N, 1);...
    repmat({'Fix'}, Regions_N*Freqs_N, 1)];
Regions = repmat([repmat({'Frontal'}, Freqs_N, 1); ...
    repmat({'Central'}, Freqs_N, 1);...
    repmat({'Parietal'}, Freqs_N, 1);...
    repmat({'Occipital'}, Freqs_N, 1)],Conditions_N,1);
Frequencies = repmat({'Low';'High'},Conditions_N*Regions_N,1);

% Low vs high frequency bands
Freqs = miniMADE_MOI_table.Sfreqs{1};
LFreq_boundaries = [2 12];
nHFreq_boundaries = [15 45];
LowFreqs_inds = [find(Freqs == LFreq_boundaries(1,1)), find(Freqs == LFreq_boundaries(1,2))]; 
HighFreqs_inds = [find(Freqs == nHFreq_boundaries(1,1)), find(Freqs == nHFreq_boundaries(1,2))];
RelPow_boundaries = [1 45];
RelPow_inds = [find(Freqs == RelPow_boundaries(1,1)), find(Freqs == RelPow_boundaries(1,2))]; 


% Loop over all data files
for ss = 1:height(miniMADE_MOI_table)

    % info
    CurrID = miniMADE_MOI_table.ID{ss};
    IDfull = repmat({CurrID}, size(Condition,1), 1);
    IDparts = strsplit(CurrID, '_');
    IDcode = repmat({IDparts{1,1}}, size(Condition,1), 1);
    Site = repmat({IDparts{1,2}}, size(Condition,1), 1);
    Sessions = repmat({IDparts{1,3}}, size(Condition,1), 1);
    if strcmp(IDparts{1,3},'aphp')
        Sessions = repmat({IDparts{1,4}}, size(Condition,1), 1);
    end
    % find clinical id
    % check if format is PSC2, and if not change it to PSC2
    ClinicalID1 = [];
    for tt = 1:height(PSCcode_conversion)
        if strcmp(num2str(PSCcode_conversion.PSC2(tt)),IDparts{1,1})
            ClinicalID1 = num2str(PSCcode_conversion.PSC2(tt));
        end
    end
    clear tt
    if isempty(ClinicalID1)
        for tt = 1:height(PSCcode_conversion)
            if strcmp(PSCcode_conversion.PSC1{tt},IDparts{1,1})
                ClinicalID1 = num2str(PSCcode_conversion.PSC2(tt));
            end
        end
    end
    clear tt
    ClinicalID = repmat({ClinicalID1}, size(Condition,1), 1);

    % Load power data if available
    if ~isempty(miniMADE_MOI_table.Power_file{ss})
        load(miniMADE_MOI_table.Power_file{ss})
    end

    % Social condition (condition x region x frequency band) %%%%%%%%%%%%%%
    S_Neps = miniMADE_MOI_table.Seps(ss);
    if S_Neps >= 20
        
        Ch_labels = DATA(1).ft_EEG_freqdata.label;
        Chind_front = find(ismember(Ch_labels,{'Fz','F3','F4'}));
        Chind_cent = find(ismember(Ch_labels,{'Cz','C3','C4'}));
        Chind_pari = find(ismember(Ch_labels,{'Pz','P3','P4'}));
        Chind_occip = find(ismember(Ch_labels,{'PO7','Oz','PO8'}));

        % calculate power for different frequencies
        % get absolute power per channel
        PowAbs_Ch = DATA(1).Abspow_chs;
            % Frontal region:
                % low freqs
                F_LF_pow = mean(PowAbs_Ch(Chind_front,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowAbs_Ch(Chind_front,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    F_HF_pow = mean(Vals,'all'); clear Vals
            % Central region:
                % low freqs
                C_LF_pow = mean(PowAbs_Ch(Chind_cent,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowAbs_Ch(Chind_cent,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    C_HF_pow = mean(Vals,'all'); clear Vals
            % Parietal region:
                % low freqs
                P_LF_pow = mean(PowAbs_Ch(Chind_pari,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowAbs_Ch(Chind_pari,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    P_HF_pow = mean(Vals,'all'); clear Vals
            % Occipital region:
                % low freqs
                O_LF_pow = mean(PowAbs_Ch(Chind_occip,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowAbs_Ch(Chind_occip,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    O_HF_pow = mean(Vals,'all'); clear Vals

            Apow = [F_LF_pow; F_HF_pow; C_LF_pow; C_HF_pow; P_LF_pow; P_HF_pow; O_LF_pow; O_HF_pow];
            clear PowAbs_Ch F_LF_pow F_HF_pow C_LF_pow C_HF_pow P_LF_pow P_HF_pow O_LF_pow O_HF_pow
 
        % log power
        % get log power per channel
        PowLog_Ch = DATA(1).ft_EEG_freqdata.logpowspctrm_avg;
            % Frontal region:
                % low freqs
                F_LF_pow = mean(PowLog_Ch(Chind_front,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowLog_Ch(Chind_front,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    F_HF_pow = mean(Vals,'all'); clear Vals
            % Central region:
                % low freqs
                C_LF_pow = mean(PowLog_Ch(Chind_cent,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowLog_Ch(Chind_cent,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    C_HF_pow = mean(Vals,'all'); clear Vals
            % Parietal region:
                % low freqs
                P_LF_pow = mean(PowLog_Ch(Chind_pari,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowLog_Ch(Chind_pari,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    P_HF_pow = mean(Vals,'all'); clear Vals
            % Occipital region:
                % low freqs
                O_LF_pow = mean(PowLog_Ch(Chind_occip,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowLog_Ch(Chind_occip,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    O_HF_pow = mean(Vals,'all'); clear Vals

            Lpow = [F_LF_pow; F_HF_pow; C_LF_pow; C_HF_pow; P_LF_pow; P_HF_pow; O_LF_pow; O_HF_pow];
            clear PowLog_Ch F_LF_pow F_HF_pow C_LF_pow C_HF_pow P_LF_pow P_HF_pow O_LF_pow O_HF_pow


        % relative power: sum(band of interest) / sum (1-45Hz) => not
        % including higher gamma
        PowAbs_Ch = DATA(1).Abspow_chs;
            % Frontal region:
                PowAbs_FrChs = mean(PowAbs_Ch(Chind_front,:),1);
                Deno = sum(PowAbs_FrChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                % low freqs
                LF_pow_sum = sum(PowAbs_FrChs(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);
                F_LF_pow = LF_pow_sum/Deno; clear LF_pow_sum
                % high freqs
                    HF_pow_sum = sum(PowAbs_FrChs(1,HighFreqs_inds(1,1):HighFreqs_inds(1,2)),2);
                    F_HF_pow = HF_pow_sum/Deno; clear HF_pow_sum
                clear Deno PowAbs_FrChs 
            % Central region:
                PowAbs_CeChs = mean(PowAbs_Ch(Chind_cent,:),1);
                Deno = sum(PowAbs_CeChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                % low freqs
                LF_pow_sum = sum(PowAbs_CeChs(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);
                C_LF_pow = LF_pow_sum/Deno; clear LF_pow_sum
                % high freqs
                    HF_pow_sum = sum(PowAbs_CeChs(1,HighFreqs_inds(1,1):HighFreqs_inds(1,2)),2);
                    C_HF_pow = HF_pow_sum/Deno; clear HF_pow_sum
                clear Deno PowAbs_CeChs
            % Parietal region:
                PowAbs_PaChs = mean(PowAbs_Ch(Chind_pari,:),1);
                Deno = sum(PowAbs_PaChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                % low freqs
                LF_pow_sum = sum(PowAbs_PaChs(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);
                P_LF_pow = LF_pow_sum/Deno; clear LF_pow_sum
                % high freqs
                    HF_pow_sum = sum(PowAbs_PaChs(1,HighFreqs_inds(1,1):HighFreqs_inds(1,2)),2);
                    P_HF_pow = HF_pow_sum/Deno; clear HF_pow_sum
                clear Deno PowAbs_PaChs
            % Occipital region:
                PowAbs_OcChs = mean(PowAbs_Ch(Chind_occip,:),1);
                Deno = sum(PowAbs_OcChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                % low freqs
                LF_pow_sum = sum(PowAbs_OcChs(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);
                O_LF_pow = LF_pow_sum/Deno; clear LF_pow_sum
                % high freqs
                    HF_pow_sum = sum(PowAbs_OcChs(1,HighFreqs_inds(1,1):HighFreqs_inds(1,2)),2);
                    O_HF_pow = HF_pow_sum/Deno; clear HF_pow_sum
                clear Deno PowAbs_OcChs
            Rpow = [F_LF_pow; F_HF_pow; C_LF_pow; C_HF_pow; P_LF_pow; P_HF_pow; O_LF_pow; O_HF_pow];
            clear PowAbs_Ch F_LF_pow F_HF_pow C_LF_pow C_HF_pow P_LF_pow P_HF_pow O_LF_pow O_HF_pow

    else % no sufficient data for this condition
        % absolute power
        Apow = nan(Regions_N*Freqs_N,1);
        % log power
        Lpow = nan(Regions_N*Freqs_N,1);
        % relative power
        Rpow = nan(Regions_N*Freqs_N,1);
    end
    % Rename vars
    S_Apow = Apow; clear Apow
    S_Lpow = Lpow; clear Lpow
    S_Rpow = Rpow; clear Rpow





    % Toy condition (condition x region x frequency band) %%%%%%%%%%%%%%%%%
    T_Neps = miniMADE_MOI_table.Teps(ss);
    if T_Neps >= 20
        
        Ch_labels = DATA(2).ft_EEG_freqdata.label;
        Chind_front = find(ismember(Ch_labels,{'Fz','F3','F4'}));
        Chind_cent = find(ismember(Ch_labels,{'Cz','C3','C4'}));
        Chind_pari = find(ismember(Ch_labels,{'Pz','P3','P4'}));
        Chind_occip = find(ismember(Ch_labels,{'PO7','Oz','PO8'}));

        % calculate power for different frequencies
        % get absolute power per channel
        PowAbs_Ch = DATA(2).Abspow_chs;
            % Frontal region:
                % low freqs
                F_LF_pow = mean(PowAbs_Ch(Chind_front,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowAbs_Ch(Chind_front,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    F_HF_pow = mean(Vals,'all'); clear Vals
            % Central region:
                % low freqs
                C_LF_pow = mean(PowAbs_Ch(Chind_cent,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowAbs_Ch(Chind_cent,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    C_HF_pow = mean(Vals,'all'); clear Vals
            % Parietal region:
                % low freqs
                P_LF_pow = mean(PowAbs_Ch(Chind_pari,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowAbs_Ch(Chind_pari,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    P_HF_pow = mean(Vals,'all'); clear Vals
            % Occipital region:
                % low freqs
                O_LF_pow = mean(PowAbs_Ch(Chind_occip,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowAbs_Ch(Chind_occip,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    O_HF_pow = mean(Vals,'all'); clear Vals

            Apow = [F_LF_pow; F_HF_pow; C_LF_pow; C_HF_pow; P_LF_pow; P_HF_pow; O_LF_pow; O_HF_pow];
            clear PowAbs_Ch F_LF_pow F_HF_pow C_LF_pow C_HF_pow P_LF_pow P_HF_pow O_LF_pow O_HF_pow
 
        % log power
        % get log power per channel
        PowLog_Ch = DATA(2).ft_EEG_freqdata.logpowspctrm_avg;
            % Frontal region:
                % low freqs
                F_LF_pow = mean(PowLog_Ch(Chind_front,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowLog_Ch(Chind_front,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    F_HF_pow = mean(Vals,'all'); clear Vals
            % Central region:
                % low freqs
                C_LF_pow = mean(PowLog_Ch(Chind_cent,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowLog_Ch(Chind_cent,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    C_HF_pow = mean(Vals,'all'); clear Vals
            % Parietal region:
                % low freqs
                P_LF_pow = mean(PowLog_Ch(Chind_pari,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowLog_Ch(Chind_pari,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    P_HF_pow = mean(Vals,'all'); clear Vals
            % Occipital region:
                % low freqs
                O_LF_pow = mean(PowLog_Ch(Chind_occip,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowLog_Ch(Chind_occip,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    O_HF_pow = mean(Vals,'all'); clear Vals

            Lpow = [F_LF_pow; F_HF_pow; C_LF_pow; C_HF_pow; P_LF_pow; P_HF_pow; O_LF_pow; O_HF_pow];
            clear PowLog_Ch F_LF_pow F_HF_pow C_LF_pow C_HF_pow P_LF_pow P_HF_pow O_LF_pow O_HF_pow


        % relative power: sum(band of interest) / sum (1-45Hz) => not
        % including higher gamma
        PowAbs_Ch = DATA(2).Abspow_chs;
            % Frontal region:
                PowAbs_FrChs = mean(PowAbs_Ch(Chind_front,:),1);
                Deno = sum(PowAbs_FrChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                % low freqs
                LF_pow_sum = sum(PowAbs_FrChs(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);
                F_LF_pow = LF_pow_sum/Deno; clear LF_pow_sum
                % high freqs
                    HF_pow_sum = sum(PowAbs_FrChs(1,HighFreqs_inds(1,1):HighFreqs_inds(1,2)),2);
                    F_HF_pow = HF_pow_sum/Deno; clear HF_pow_sum
                clear Deno PowAbs_FrChs 
            % Central region:
                PowAbs_CeChs = mean(PowAbs_Ch(Chind_cent,:),1);
                Deno = sum(PowAbs_CeChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                % low freqs
                LF_pow_sum = sum(PowAbs_CeChs(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);
                C_LF_pow = LF_pow_sum/Deno; clear LF_pow_sum
                % high freqs
                    HF_pow_sum = sum(PowAbs_CeChs(1,HighFreqs_inds(1,1):HighFreqs_inds(1,2)),2);
                    C_HF_pow = HF_pow_sum/Deno; clear HF_pow_sum
                clear Deno PowAbs_CeChs
            % Parietal region:
                PowAbs_PaChs = mean(PowAbs_Ch(Chind_pari,:),1);
                Deno = sum(PowAbs_PaChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                % low freqs
                LF_pow_sum = sum(PowAbs_PaChs(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);
                P_LF_pow = LF_pow_sum/Deno; clear LF_pow_sum
                % high freqs
                    HF_pow_sum = sum(PowAbs_PaChs(1,HighFreqs_inds(1,1):HighFreqs_inds(1,2)),2);
                    P_HF_pow = HF_pow_sum/Deno; clear HF_pow_sum
                clear Deno PowAbs_PaChs
            % Occipital region:
                PowAbs_OcChs = mean(PowAbs_Ch(Chind_occip,:),1);
                Deno = sum(PowAbs_OcChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                % low freqs
                LF_pow_sum = sum(PowAbs_OcChs(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);
                O_LF_pow = LF_pow_sum/Deno; clear LF_pow_sum
                % high freqs
                    HF_pow_sum = sum(PowAbs_OcChs(1,HighFreqs_inds(1,1):HighFreqs_inds(1,2)),2);
                    O_HF_pow = HF_pow_sum/Deno; clear HF_pow_sum
                clear Deno PowAbs_OcChs
            Rpow = [F_LF_pow; F_HF_pow; C_LF_pow; C_HF_pow; P_LF_pow; P_HF_pow; O_LF_pow; O_HF_pow];
            clear PowAbs_Ch F_LF_pow F_HF_pow C_LF_pow C_HF_pow P_LF_pow P_HF_pow O_LF_pow O_HF_pow

    else % no sufficient data for this condition
        % absolute power
        Apow = nan(Regions_N*Freqs_N,1);
        % log power
        Lpow = nan(Regions_N*Freqs_N,1);
        % relative power
        Rpow = nan(Regions_N*Freqs_N,1);
    end
    % Rename vars
    T_Apow = Apow; clear Apow
    T_Lpow = Lpow; clear Lpow
    T_Rpow = Rpow; clear Rpow






    % Abstract condition (condition x region x frequency band) %%%%%%%%%%%%
    A_Neps = miniMADE_MOI_table.Aeps(ss);
    if A_Neps >= 20

        Ch_labels = DATA(3).ft_EEG_freqdata.label;
        Chind_front = find(ismember(Ch_labels,{'Fz','F3','F4'}));
        Chind_cent = find(ismember(Ch_labels,{'Cz','C3','C4'}));
        Chind_pari = find(ismember(Ch_labels,{'Pz','P3','P4'}));
        Chind_occip = find(ismember(Ch_labels,{'PO7','Oz','PO8'}));

        % calculate power for different frequencies
        % get absolute power per channel
        PowAbs_Ch = DATA(3).Abspow_chs;
            % Frontal region:
                % low freqs
                F_LF_pow = mean(PowAbs_Ch(Chind_front,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowAbs_Ch(Chind_front,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    F_HF_pow = mean(Vals,'all'); clear Vals
            % Central region:
                % low freqs
                C_LF_pow = mean(PowAbs_Ch(Chind_cent,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowAbs_Ch(Chind_cent,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    C_HF_pow = mean(Vals,'all'); clear Vals
            % Parietal region:
                % low freqs
                P_LF_pow = mean(PowAbs_Ch(Chind_pari,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowAbs_Ch(Chind_pari,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    P_HF_pow = mean(Vals,'all'); clear Vals
            % Occipital region:
                % low freqs
                O_LF_pow = mean(PowAbs_Ch(Chind_occip,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowAbs_Ch(Chind_occip,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    O_HF_pow = mean(Vals,'all'); clear Vals

            Apow = [F_LF_pow; F_HF_pow; C_LF_pow; C_HF_pow; P_LF_pow; P_HF_pow; O_LF_pow; O_HF_pow];
            clear PowAbs_Ch F_LF_pow F_HF_pow C_LF_pow C_HF_pow P_LF_pow P_HF_pow O_LF_pow O_HF_pow
 
        % log power
        % get log power per channel
        PowLog_Ch = DATA(3).ft_EEG_freqdata.logpowspctrm_avg;
            % Frontal region:
                % low freqs
                F_LF_pow = mean(PowLog_Ch(Chind_front,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowLog_Ch(Chind_front,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    F_HF_pow = mean(Vals,'all'); clear Vals
            % Central region:
                % low freqs
                C_LF_pow = mean(PowLog_Ch(Chind_cent,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowLog_Ch(Chind_cent,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    C_HF_pow = mean(Vals,'all'); clear Vals
            % Parietal region:
                % low freqs
                P_LF_pow = mean(PowLog_Ch(Chind_pari,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowLog_Ch(Chind_pari,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    P_HF_pow = mean(Vals,'all'); clear Vals
            % Occipital region:
                % low freqs
                O_LF_pow = mean(PowLog_Ch(Chind_occip,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowLog_Ch(Chind_occip,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    O_HF_pow = mean(Vals,'all'); clear Vals

            Lpow = [F_LF_pow; F_HF_pow; C_LF_pow; C_HF_pow; P_LF_pow; P_HF_pow; O_LF_pow; O_HF_pow];
            clear PowLog_Ch F_LF_pow F_HF_pow C_LF_pow C_HF_pow P_LF_pow P_HF_pow O_LF_pow O_HF_pow


        % relative power: sum(band of interest) / sum (1-45Hz) => not
        % including higher gamma
        PowAbs_Ch = DATA(3).Abspow_chs;
            % Frontal region:
                PowAbs_FrChs = mean(PowAbs_Ch(Chind_front,:),1);
                Deno = sum(PowAbs_FrChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                % low freqs
                LF_pow_sum = sum(PowAbs_FrChs(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);
                F_LF_pow = LF_pow_sum/Deno; clear LF_pow_sum
                % high freqs
                    HF_pow_sum = sum(PowAbs_FrChs(1,HighFreqs_inds(1,1):HighFreqs_inds(1,2)),2);
                    F_HF_pow = HF_pow_sum/Deno; clear HF_pow_sum
                clear Deno PowAbs_FrChs 
            % Central region:
                PowAbs_CeChs = mean(PowAbs_Ch(Chind_cent,:),1);
                Deno = sum(PowAbs_CeChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                % low freqs
                LF_pow_sum = sum(PowAbs_CeChs(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);
                C_LF_pow = LF_pow_sum/Deno; clear LF_pow_sum
                % high freqs
                    HF_pow_sum = sum(PowAbs_CeChs(1,HighFreqs_inds(1,1):HighFreqs_inds(1,2)),2);
                    C_HF_pow = HF_pow_sum/Deno; clear HF_pow_sum
                clear Deno PowAbs_CeChs
            % Parietal region:
                PowAbs_PaChs = mean(PowAbs_Ch(Chind_pari,:),1);
                Deno = sum(PowAbs_PaChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                % low freqs
                LF_pow_sum = sum(PowAbs_PaChs(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);
                P_LF_pow = LF_pow_sum/Deno; clear LF_pow_sum
                % high freqs
                    HF_pow_sum = sum(PowAbs_PaChs(1,HighFreqs_inds(1,1):HighFreqs_inds(1,2)),2);
                    P_HF_pow = HF_pow_sum/Deno; clear HF_pow_sum
                clear Deno PowAbs_PaChs
            % Occipital region:
                PowAbs_OcChs = mean(PowAbs_Ch(Chind_occip,:),1);
                Deno = sum(PowAbs_OcChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                % low freqs
                LF_pow_sum = sum(PowAbs_OcChs(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);
                O_LF_pow = LF_pow_sum/Deno; clear LF_pow_sum
                % high freqs
                    HF_pow_sum = sum(PowAbs_OcChs(1,HighFreqs_inds(1,1):HighFreqs_inds(1,2)),2);
                    O_HF_pow = HF_pow_sum/Deno; clear HF_pow_sum
                clear Deno PowAbs_OcChs
            Rpow = [F_LF_pow; F_HF_pow; C_LF_pow; C_HF_pow; P_LF_pow; P_HF_pow; O_LF_pow; O_HF_pow];
            clear PowAbs_Ch F_LF_pow F_HF_pow C_LF_pow C_HF_pow P_LF_pow P_HF_pow O_LF_pow O_HF_pow

    else % no sufficient data for this condition
        % absolute power
        Apow = nan(Regions_N*Freqs_N,1);
        % log power
        Lpow = nan(Regions_N*Freqs_N,1);
        % relative power
        Rpow = nan(Regions_N*Freqs_N,1);
    end
    % Rename vars
    A_Apow = Apow; clear Apow
    A_Lpow = Lpow; clear Lpow
    A_Rpow = Rpow; clear Rpow






    % Fixation condition (condition x region x frequency band) %%%%%%%%%%%%
    F_Neps = miniMADE_MOI_table.Feps(ss);
    if F_Neps >= 20
        
        Ch_labels = DATA(4).ft_EEG_freqdata.label;
        Chind_front = find(ismember(Ch_labels,{'Fz','F3','F4'}));
        Chind_cent = find(ismember(Ch_labels,{'Cz','C3','C4'}));
        Chind_pari = find(ismember(Ch_labels,{'Pz','P3','P4'}));
        Chind_occip = find(ismember(Ch_labels,{'PO7','Oz','PO8'}));

        % calculate power for different frequencies
        % get absolute power per channel
        PowAbs_Ch = DATA(4).Abspow_chs;
            % Frontal region:
                % low freqs
                F_LF_pow = mean(PowAbs_Ch(Chind_front,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowAbs_Ch(Chind_front,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    F_HF_pow = mean(Vals,'all'); clear Vals
            % Central region:
                % low freqs
                C_LF_pow = mean(PowAbs_Ch(Chind_cent,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowAbs_Ch(Chind_cent,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    C_HF_pow = mean(Vals,'all'); clear Vals
            % Parietal region:
                % low freqs
                P_LF_pow = mean(PowAbs_Ch(Chind_pari,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowAbs_Ch(Chind_pari,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    P_HF_pow = mean(Vals,'all'); clear Vals
            % Occipital region:
                % low freqs
                O_LF_pow = mean(PowAbs_Ch(Chind_occip,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowAbs_Ch(Chind_occip,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    O_HF_pow = mean(Vals,'all'); clear Vals

            Apow = [F_LF_pow; F_HF_pow; C_LF_pow; C_HF_pow; P_LF_pow; P_HF_pow; O_LF_pow; O_HF_pow];
            clear PowAbs_Ch F_LF_pow F_HF_pow C_LF_pow C_HF_pow P_LF_pow P_HF_pow O_LF_pow O_HF_pow
 
        % log power
        % get log power per channel
        PowLog_Ch = DATA(4).ft_EEG_freqdata.logpowspctrm_avg;
            % Frontal region:
                % low freqs
                F_LF_pow = mean(PowLog_Ch(Chind_front,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowLog_Ch(Chind_front,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    F_HF_pow = mean(Vals,'all'); clear Vals
            % Central region:
                % low freqs
                C_LF_pow = mean(PowLog_Ch(Chind_cent,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowLog_Ch(Chind_cent,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    C_HF_pow = mean(Vals,'all'); clear Vals
            % Parietal region:
                % low freqs
                P_LF_pow = mean(PowLog_Ch(Chind_pari,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowLog_Ch(Chind_pari,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    P_HF_pow = mean(Vals,'all'); clear Vals
            % Occipital region:
                % low freqs
                O_LF_pow = mean(PowLog_Ch(Chind_occip,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),'all');
                % high freqs
                    Vals = [PowLog_Ch(Chind_occip,HighFreqs_inds(1,1):HighFreqs_inds(1,2))];
                    O_HF_pow = mean(Vals,'all'); clear Vals

            Lpow = [F_LF_pow; F_HF_pow; C_LF_pow; C_HF_pow; P_LF_pow; P_HF_pow; O_LF_pow; O_HF_pow];
            clear PowLog_Ch F_LF_pow F_HF_pow C_LF_pow C_HF_pow P_LF_pow P_HF_pow O_LF_pow O_HF_pow


        % relative power: sum(band of interest) / sum (1-45Hz) => not
        % including higher gamma
        PowAbs_Ch = DATA(4).Abspow_chs;
            % Frontal region:
                PowAbs_FrChs = mean(PowAbs_Ch(Chind_front,:),1);
                Deno = sum(PowAbs_FrChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                % low freqs
                LF_pow_sum = sum(PowAbs_FrChs(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);
                F_LF_pow = LF_pow_sum/Deno; clear LF_pow_sum
                % high freqs
                    HF_pow_sum = sum(PowAbs_FrChs(1,HighFreqs_inds(1,1):HighFreqs_inds(1,2)),2);
                    F_HF_pow = HF_pow_sum/Deno; clear HF_pow_sum
                clear Deno PowAbs_FrChs 
            % Central region:
                PowAbs_CeChs = mean(PowAbs_Ch(Chind_cent,:),1);
                Deno = sum(PowAbs_CeChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                % low freqs
                LF_pow_sum = sum(PowAbs_CeChs(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);
                C_LF_pow = LF_pow_sum/Deno; clear LF_pow_sum
                % high freqs
                    HF_pow_sum = sum(PowAbs_CeChs(1,HighFreqs_inds(1,1):HighFreqs_inds(1,2)),2);
                    C_HF_pow = HF_pow_sum/Deno; clear HF_pow_sum
                clear Deno PowAbs_CeChs
            % Parietal region:
                PowAbs_PaChs = mean(PowAbs_Ch(Chind_pari,:),1);
                Deno = sum(PowAbs_PaChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                % low freqs
                LF_pow_sum = sum(PowAbs_PaChs(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);
                P_LF_pow = LF_pow_sum/Deno; clear LF_pow_sum
                % high freqs
                    HF_pow_sum = sum(PowAbs_PaChs(1,HighFreqs_inds(1,1):HighFreqs_inds(1,2)),2);
                    P_HF_pow = HF_pow_sum/Deno; clear HF_pow_sum
                clear Deno PowAbs_PaChs
            % Occipital region:
                PowAbs_OcChs = mean(PowAbs_Ch(Chind_occip,:),1);
                Deno = sum(PowAbs_OcChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                % low freqs
                LF_pow_sum = sum(PowAbs_OcChs(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);
                O_LF_pow = LF_pow_sum/Deno; clear LF_pow_sum
                % high freqs
                    HF_pow_sum = sum(PowAbs_OcChs(1,HighFreqs_inds(1,1):HighFreqs_inds(1,2)),2);
                    O_HF_pow = HF_pow_sum/Deno; clear HF_pow_sum
                clear Deno PowAbs_OcChs
            Rpow = [F_LF_pow; F_HF_pow; C_LF_pow; C_HF_pow; P_LF_pow; P_HF_pow; O_LF_pow; O_HF_pow];
            clear PowAbs_Ch F_LF_pow F_HF_pow C_LF_pow C_HF_pow P_LF_pow P_HF_pow O_LF_pow O_HF_pow

            clear DATA Ch_labels

    else % no sufficient data for this condition
        % absolute power
        Apow = nan(Regions_N*Freqs_N,1);
        % log power
        Lpow = nan(Regions_N*Freqs_N,1);
        % relative power
        Rpow = nan(Regions_N*Freqs_N,1);
    end
    % Rename vars
    F_Apow = Apow; clear Apow
    F_Lpow = Lpow; clear Lpow
    F_Rpow = Rpow; clear Rpow
   

    % Create table with new data
    % Columns: IDfull, ID, session, condition (4), region (4) frequency band (2), N trials, absolute power, log power, relative power
    % Ntrials
    N_trials = [repmat(S_Neps, Regions_N*Freqs_N,1); repmat(T_Neps, Regions_N*Freqs_N,1); ...
        repmat(A_Neps, Regions_N*Freqs_N,1); repmat(F_Neps, Regions_N*Freqs_N,1)];
    % Power
    Abs_Power = [S_Apow; T_Apow; A_Apow; F_Apow];
    Log_Power = [S_Lpow; T_Lpow; A_Lpow; F_Lpow];
    Rel_Power = [S_Rpow; T_Rpow; A_Rpow; F_Rpow];
    % assemble table
    dataTable_currID = table(ClinicalID, IDfull, IDcode, Site, Sessions, Condition, Regions, Frequencies, N_trials, Abs_Power, Log_Power, Rel_Power);
    % add current ID data to table
    cd(MOI_location)
    if exist('RCT1_Power_LvsHfreqs_longformat_CxRxFB_upto45Hz.mat','file') == 2
        load RCT1_Power_LvsHfreqs_longformat_CxRxFB_upto45Hz.mat
        NewRowInds = [(ss-1)*height(dataTable_currID)+1 ss*height(dataTable_currID)];
        RCT1_Power_longformat(NewRowInds(1,1):NewRowInds(1,2),:) = dataTable_currID;
        clear NewRowInds
    else % first ID
        RCT1_Power_longformat = dataTable_currID;
    end

    % save the table
    save('RCT1_Power_LvsHfreqs_longformat_CxRxFB_upto45Hz.mat','RCT1_Power_longformat')
    
    % clear up
    clear ClinicalID ClinicalID1 IDfull IDparts IDcode Site Sessions N_trials Abs_Power Log_Power Rel_Power 
    clear S_Apow T_Apow A_Apow F_Apow S_Lpow T_Lpow A_Lpow F_Lpow S_Rpow T_Rpow A_Rpow F_Rpow
    clear CurrID dataTable_currID S_Neps T_Neps A_Neps F_Neps PowCurr_abs
    clear RCT1_Power_longformat
    clear Chind_front Chind_cent Chind_pari Chind_occip


end


clear Freqs Freq_boundaries Freqs_inds Condition Frequencies
clear Freqs_N Regions Regions_N Conditions_N nHFreq_boundaries HighFreqs_inds LFreq_boundaries LowFreqs_inds ss RelPow_boundaries RelPow_inds
clear HFBUMPS_present PSCcode_conversion Subjects_with_HFbumps miniMADE_MOI_table

%% Save table into csv format for Rstudio
cd(MOI_location)
load RCT1_Power_LvsHfreqs_longformat_CxRxFB.mat
% save as csv file
writetable(RCT1_Power_longformat, 'RCT1_Power_LvsHfreqs_longformat_CxRxFB.csv')
clear RCT1_Power_longformat 


%% For canonical frequency bands %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Columns: IDfull, ID, session, condition (4), region (4), frequency band (6), N trials, absolute power, log power, relative power
% frequency bands: 
% delta band: 2-5Hz, theta band: 6-7Hz, alpha band: 8-12Hz, beta band: 15-25Hz, lower gamma band: 35-45Hz
% 96 rows per ID

clear all 
MOI_location = 'xxx/DATA/01A_PreprocPower';
cd(MOI_location)
load miniMADE_MeasuresOfInterest_report.mat

PSCcode_conversion = readtable('xxx/rct1_psc1_psc2.xlsx');

% Categories:
Conditions_N = 4;
Regions_N = 4;
Freqs_N = 5;

% Settings
Freqs = miniMADE_MOI_table.Sfreqs{1};
Freq_boundaries = [2 5; 6 7; 8 12; 15 25; 35 45];
Freqs_inds = [find(Freqs == Freq_boundaries(1,1)), find(Freqs == Freq_boundaries(1,2)); ...
    find(Freqs == Freq_boundaries(2,1)), find(Freqs == Freq_boundaries(2,2)); ...
    find(Freqs == Freq_boundaries(3,1)), find(Freqs == Freq_boundaries(3,2)); ...
    find(Freqs == Freq_boundaries(4,1)), find(Freqs == Freq_boundaries(4,2)); ...
    find(Freqs == Freq_boundaries(5,1)), find(Freqs == Freq_boundaries(5,2))];
RelPow_boundaries = [1 45];
RelPow_inds = [find(Freqs == RelPow_boundaries(1,1)), find(Freqs == RelPow_boundaries(1,2))]; 
Condition = [repmat({'Social'}, Regions_N*Freqs_N, 1); ...
    repmat({'Toy'}, Regions_N*Freqs_N, 1);...
    repmat({'Abstract'}, Regions_N*Freqs_N, 1);...
    repmat({'Fix'}, Regions_N*Freqs_N, 1)];
Regions = repmat([repmat({'Frontal'}, Freqs_N, 1); ...
    repmat({'Central'}, Freqs_N, 1);...
    repmat({'Parietal'}, Freqs_N, 1);...
    repmat({'Occipital'}, Freqs_N, 1)],Conditions_N,1);
Frequency_bands = {'Delta','Theta','Alpha','Beta','Lower Gamma'}';
Frequencies = repmat(Frequency_bands,Conditions_N*Regions_N,1);


% Loop through participants
for ss = 1:height(miniMADE_MOI_table)

    % info
    CurrID = miniMADE_MOI_table.ID{ss};
    IDfull = repmat({CurrID}, size(Condition,1), 1);
    IDparts = strsplit(CurrID, '_');
    IDcode = repmat({IDparts{1,1}}, size(Condition,1), 1);
    Site = repmat({IDparts{1,2}}, size(Condition,1), 1);
    Sessions = repmat({IDparts{1,3}}, size(Condition,1), 1);
    if strcmp(IDparts{1,3},'aphp')
        Sessions = repmat({IDparts{1,4}}, size(Condition,1), 1);
    end
    % find clinical id
    % check if format is PSC2, and if not change it to PSC2
    ClinicalID1 = [];
    for tt = 1:height(PSCcode_conversion)
        if strcmp(num2str(PSCcode_conversion.PSC2(tt)),IDparts{1,1})
            ClinicalID1 = num2str(PSCcode_conversion.PSC2(tt));
        end
    end
    clear tt
    if isempty(ClinicalID1)
        for tt = 1:height(PSCcode_conversion)
            if strcmp(PSCcode_conversion.PSC1{tt},IDparts{1,1})
                ClinicalID1 = num2str(PSCcode_conversion.PSC2(tt));
            end
        end
    end
    clear tt
    ClinicalID = repmat({ClinicalID1}, size(Condition,1), 1);

    % Load power data if available
    if ~isempty(miniMADE_MOI_table.Power_file{ss})
        load(miniMADE_MOI_table.Power_file{ss})
    end

    % Social condition (condition x region x frequency band) %%%%%%%%%%%%%%
    S_Neps = miniMADE_MOI_table.Seps(ss);
    if S_Neps >= 20

        Ch_labels = DATA(1).ft_EEG_freqdata.label;
        Chind_front = find(ismember(Ch_labels,{'Fz','F3','F4'}));
        Chind_cent = find(ismember(Ch_labels,{'Cz','C3','C4'}));
        Chind_pari = find(ismember(Ch_labels,{'Pz','P3','P4'}));
        Chind_occip = find(ismember(Ch_labels,{'PO7','Oz','PO8'}));

        % calculate power for different frequencies
        % get absolute power per channel
        PowAbs_Ch = DATA(1).Abspow_chs;
            % Frontal region:
                F_Apow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(F_Apow,1)
                    F_Apow(ff,1) = mean(PowAbs_Ch(Chind_front,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff
            % Central region:
                C_Apow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(C_Apow,1)
                    C_Apow(ff,1) = mean(PowAbs_Ch(Chind_cent,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff
            % Parietal region:
                P_Apow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(P_Apow,1)
                    P_Apow(ff,1) = mean(PowAbs_Ch(Chind_pari,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff
            % Occipital region:
                O_Apow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(O_Apow,1)
                    O_Apow(ff,1) = mean(PowAbs_Ch(Chind_occip,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff

            Apow = [F_Apow; C_Apow; P_Apow; O_Apow];
            clear PowAbs_Ch F_Apow C_Apow P_Apow O_Apow
 
        % log power
        % get log power per channel
        PowLog_Ch = DATA(1).ft_EEG_freqdata.logpowspctrm_avg;
            % Frontal region:
                F_Lpow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(F_Lpow,1)
                    F_Lpow(ff,1) = mean(PowLog_Ch(Chind_front,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff
            % Central region:
                C_Lpow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(C_Lpow,1)
                    C_Lpow(ff,1) = mean(PowLog_Ch(Chind_cent,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff
            % Parietal region:
                P_Lpow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(P_Lpow,1)
                    P_Lpow(ff,1) = mean(PowLog_Ch(Chind_pari,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff
            % Occipital region:
                O_Lpow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(O_Lpow,1)
                    O_Lpow(ff,1) = mean(PowLog_Ch(Chind_occip,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff

            Lpow = [F_Lpow; C_Lpow; P_Lpow; O_Lpow];
            clear PowLog_Ch F_Lpow C_Lpow P_Lpow O_Lpow

            

        % relative power: sum(band of interest) / sum (1-45Hz) => not
        % including higher gamma
        PowAbs_Ch = DATA(1).Abspow_chs;
            % Frontal region:
                PowAbs_FrChs = mean(PowAbs_Ch(Chind_front,:),1);
                Deno = sum(PowAbs_FrChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                F_Rpow = nan(Freqs_N,1);
                for ff = 1:5 
                    Pow_sum = sum(PowAbs_FrChs(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    F_Rpow(ff) = Pow_sum / Deno;
                    clear Pow_sum
                end
                clear ff PowAbs_FrChs Deno
            % Central region:
                PowAbs_CeChs = mean(PowAbs_Ch(Chind_cent,:),1);
                Deno = sum(PowAbs_CeChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                C_Rpow = nan(Freqs_N,1);
                for ff = 1:5 
                    Pow_sum = sum(PowAbs_CeChs(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    C_Rpow(ff) = Pow_sum / Deno;
                    clear Pow_sum
                end
                clear ff PowAbs_CeChs Deno
            % Parietal region:
                PowAbs_PaChs = mean(PowAbs_Ch(Chind_pari,:),1);
                Deno = sum(PowAbs_PaChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                P_Rpow = nan(Freqs_N,1);
                for ff = 1:5 
                    Pow_sum = sum(PowAbs_PaChs(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    P_Rpow(ff) = Pow_sum / Deno;
                    clear Pow_sum
                end
                clear ff PowAbs_PaChs Deno
            % Occipital region:
                PowAbs_OcChs = mean(PowAbs_Ch(Chind_occip,:),1);
                Deno = sum(PowAbs_OcChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                O_Rpow = nan(Freqs_N,1);
                for ff = 1:5 
                    Pow_sum = sum(PowAbs_OcChs(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    O_Rpow(ff) = Pow_sum / Deno;
                    clear Pow_sum
                end
                clear ff PowAbs_OcChs Deno

            Rpow = [F_Rpow; C_Rpow; P_Rpow; O_Rpow];
            clear PowAbs_Ch F_LF_pow F_HF_pow C_LF_pow C_HF_pow P_LF_pow P_HF_pow O_LF_pow O_HF_pow

            clear Ch_labels Chind_front Chind_cent Chind_pari Chind_occip

    else % no sufficient data for this condition
        % absolute power
        Apow = nan(Regions_N*Freqs_N,1);
        % log power
        Lpow = nan(Regions_N*Freqs_N,1);
        % relative power
        Rpow = nan(Regions_N*Freqs_N,1);
    end
    % Rename vars
    S_Apow = Apow; clear Apow
    S_Lpow = Lpow; clear Lpow
    S_Rpow = Rpow; clear Rpow







    % Toy condition (condition x region x frequency band) %%%%%%%%%%%%%%%%%
    T_Neps = miniMADE_MOI_table.Teps(ss);
    if T_Neps >= 20

        Ch_labels = DATA(2).ft_EEG_freqdata.label;
        Chind_front = find(ismember(Ch_labels,{'Fz','F3','F4'}));
        Chind_cent = find(ismember(Ch_labels,{'Cz','C3','C4'}));
        Chind_pari = find(ismember(Ch_labels,{'Pz','P3','P4'}));
        Chind_occip = find(ismember(Ch_labels,{'PO7','Oz','PO8'}));

        % calculate power for different frequencies
        % get absolute power per channel
        PowAbs_Ch = DATA(2).Abspow_chs;
            % Frontal region:
                F_Apow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(F_Apow,1)
                    F_Apow(ff,1) = mean(PowAbs_Ch(Chind_front,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff
            % Central region:
                C_Apow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(C_Apow,1)
                    C_Apow(ff,1) = mean(PowAbs_Ch(Chind_cent,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff
            % Parietal region:
                P_Apow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(P_Apow,1)
                    P_Apow(ff,1) = mean(PowAbs_Ch(Chind_pari,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff
            % Occipital region:
                O_Apow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(O_Apow,1)
                    O_Apow(ff,1) = mean(PowAbs_Ch(Chind_occip,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff

            Apow = [F_Apow; C_Apow; P_Apow; O_Apow];
            clear PowAbs_Ch F_Apow C_Apow P_Apow O_Apow
 
        % log power
        % get log power per channel
        PowLog_Ch = DATA(2).ft_EEG_freqdata.logpowspctrm_avg;
            % Frontal region:
                F_Lpow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(F_Lpow,1)
                    F_Lpow(ff,1) = mean(PowLog_Ch(Chind_front,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff
            % Central region:
                C_Lpow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(C_Lpow,1)
                    C_Lpow(ff,1) = mean(PowLog_Ch(Chind_cent,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff
            % Parietal region:
                P_Lpow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(P_Lpow,1)
                    P_Lpow(ff,1) = mean(PowLog_Ch(Chind_pari,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff
            % Occipital region:
                O_Lpow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(O_Lpow,1)
                    O_Lpow(ff,1) = mean(PowLog_Ch(Chind_occip,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff

            Lpow = [F_Lpow; C_Lpow; P_Lpow; O_Lpow];
            clear PowLog_Ch F_Lpow C_Lpow P_Lpow O_Lpow

            

        % relative power: sum(band of interest) / sum (1-45Hz) => not
        % including higher gamma
        PowAbs_Ch = DATA(2).Abspow_chs;
            % Frontal region:
                PowAbs_FrChs = mean(PowAbs_Ch(Chind_front,:),1);
                Deno = sum(PowAbs_FrChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                F_Rpow = nan(Freqs_N,1);
                for ff = 1:5 %exclude higher gamma at end
                    Pow_sum = sum(PowAbs_FrChs(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    F_Rpow(ff) = Pow_sum / Deno;
                    clear Pow_sum
                end
                clear ff PowAbs_FrChs Deno
            % Central region:
                PowAbs_CeChs = mean(PowAbs_Ch(Chind_cent,:),1);
                Deno = sum(PowAbs_CeChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                C_Rpow = nan(Freqs_N,1);
                for ff = 1:5 %exclude higher gamma at end
                    Pow_sum = sum(PowAbs_CeChs(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    C_Rpow(ff) = Pow_sum / Deno;
                    clear Pow_sum
                end
                clear ff PowAbs_CeChs Deno
            % Parietal region:
                PowAbs_PaChs = mean(PowAbs_Ch(Chind_pari,:),1);
                Deno = sum(PowAbs_PaChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                P_Rpow = nan(Freqs_N,1);
                for ff = 1:5 %exclude higher gamma at end
                    Pow_sum = sum(PowAbs_PaChs(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    P_Rpow(ff) = Pow_sum / Deno;
                    clear Pow_sum
                end
                clear ff PowAbs_PaChs Deno
            % Occipital region:
                PowAbs_OcChs = mean(PowAbs_Ch(Chind_occip,:),1);
                Deno = sum(PowAbs_OcChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                O_Rpow = nan(Freqs_N,1);
                for ff = 1:5 %exclude higher gamma at end
                    Pow_sum = sum(PowAbs_OcChs(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    O_Rpow(ff) = Pow_sum / Deno;
                    clear Pow_sum
                end
                clear ff PowAbs_OcChs Deno

            Rpow = [F_Rpow; C_Rpow; P_Rpow; O_Rpow];
            clear PowAbs_Ch F_LF_pow F_HF_pow C_LF_pow C_HF_pow P_LF_pow P_HF_pow O_LF_pow O_HF_pow

            clear Ch_labels Chind_front Chind_cent Chind_pari Chind_occip

    else % no sufficient data for this condition
        % absolute power
        Apow = nan(Regions_N*Freqs_N,1);
        % log power
        Lpow = nan(Regions_N*Freqs_N,1);
        % relative power
        Rpow = nan(Regions_N*Freqs_N,1);
    end
    % Rename vars
    T_Apow = Apow; clear Apow
    T_Lpow = Lpow; clear Lpow
    T_Rpow = Rpow; clear Rpow



    % Abstract condition (condition x region x frequency band) %%%%%%%%%%%%
    A_Neps = miniMADE_MOI_table.Aeps(ss);
    if A_Neps >= 20

        Ch_labels = DATA(3).ft_EEG_freqdata.label;
        Chind_front = find(ismember(Ch_labels,{'Fz','F3','F4'}));
        Chind_cent = find(ismember(Ch_labels,{'Cz','C3','C4'}));
        Chind_pari = find(ismember(Ch_labels,{'Pz','P3','P4'}));
        Chind_occip = find(ismember(Ch_labels,{'PO7','Oz','PO8'}));

        % calculate power for different frequencies
        % get absolute power per channel
        PowAbs_Ch = DATA(3).Abspow_chs;
            % Frontal region:
                F_Apow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(F_Apow,1)
                    F_Apow(ff,1) = mean(PowAbs_Ch(Chind_front,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff
            % Central region:
                C_Apow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(C_Apow,1)
                    C_Apow(ff,1) = mean(PowAbs_Ch(Chind_cent,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff
            % Parietal region:
                P_Apow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(P_Apow,1)
                    P_Apow(ff,1) = mean(PowAbs_Ch(Chind_pari,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff
            % Occipital region:
                O_Apow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(O_Apow,1)
                    O_Apow(ff,1) = mean(PowAbs_Ch(Chind_occip,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff

            Apow = [F_Apow; C_Apow; P_Apow; O_Apow];
            clear PowAbs_Ch F_Apow C_Apow P_Apow O_Apow
 
        % log power
        % get log power per channel
        PowLog_Ch = DATA(3).ft_EEG_freqdata.logpowspctrm_avg;
            % Frontal region:
                F_Lpow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(F_Lpow,1)
                    F_Lpow(ff,1) = mean(PowLog_Ch(Chind_front,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff
            % Central region:
                C_Lpow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(C_Lpow,1)
                    C_Lpow(ff,1) = mean(PowLog_Ch(Chind_cent,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff
            % Parietal region:
                P_Lpow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(P_Lpow,1)
                    P_Lpow(ff,1) = mean(PowLog_Ch(Chind_pari,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff
            % Occipital region:
                O_Lpow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(O_Lpow,1)
                    O_Lpow(ff,1) = mean(PowLog_Ch(Chind_occip,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff

            Lpow = [F_Lpow; C_Lpow; P_Lpow; O_Lpow];
            clear PowLog_Ch F_Lpow C_Lpow P_Lpow O_Lpow

            

        % relative power: sum(band of interest) / sum (1-45Hz) => not
        % including higher gamma
        PowAbs_Ch = DATA(3).Abspow_chs;
            % Frontal region:
                PowAbs_FrChs = mean(PowAbs_Ch(Chind_front,:),1);
                Deno = sum(PowAbs_FrChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                F_Rpow = nan(Freqs_N,1);
                for ff = 1:5 %exclude higher gamma at end
                    Pow_sum = sum(PowAbs_FrChs(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    F_Rpow(ff) = Pow_sum / Deno;
                    clear Pow_sum
                end
                clear ff PowAbs_FrChs Deno
            % Central region:
                PowAbs_CeChs = mean(PowAbs_Ch(Chind_cent,:),1);
                Deno = sum(PowAbs_CeChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                C_Rpow = nan(Freqs_N,1);
                for ff = 1:5 %exclude higher gamma at end
                    Pow_sum = sum(PowAbs_CeChs(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    C_Rpow(ff) = Pow_sum / Deno;
                    clear Pow_sum
                end
                clear ff PowAbs_CeChs Deno
            % Parietal region:
                PowAbs_PaChs = mean(PowAbs_Ch(Chind_pari,:),1);
                Deno = sum(PowAbs_PaChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                P_Rpow = nan(Freqs_N,1);
                for ff = 1:5 %exclude higher gamma at end
                    Pow_sum = sum(PowAbs_PaChs(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    P_Rpow(ff) = Pow_sum / Deno;
                    clear Pow_sum
                end
                clear ff PowAbs_PaChs Deno
            % Occipital region:
                PowAbs_OcChs = mean(PowAbs_Ch(Chind_occip,:),1);
                Deno = sum(PowAbs_OcChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                O_Rpow = nan(Freqs_N,1);
                for ff = 1:5 %exclude higher gamma at end
                    Pow_sum = sum(PowAbs_OcChs(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    O_Rpow(ff) = Pow_sum / Deno;
                    clear Pow_sum
                end
                clear ff PowAbs_OcChs Deno

            Rpow = [F_Rpow; C_Rpow; P_Rpow; O_Rpow];
            clear PowAbs_Ch F_Rpow C_Rpow P_Rpow O_Rpow

            clear Ch_labels Chind_front Chind_cent Chind_pari Chind_occip

    else % no sufficient data for this condition
        % absolute power
        Apow = nan(Regions_N*Freqs_N,1);
        % log power
        Lpow = nan(Regions_N*Freqs_N,1);
        % relative power
        Rpow = nan(Regions_N*Freqs_N,1);
    end
    % Rename vars
    A_Apow = Apow; clear Apow
    A_Lpow = Lpow; clear Lpow
    A_Rpow = Rpow; clear Rpow




    % Fixation condition (condition x region x frequency band) %%%%%%%%%%%%
    F_Neps = miniMADE_MOI_table.Feps(ss);
    if F_Neps >= 20
        
        Ch_labels = DATA(4).ft_EEG_freqdata.label;
        Chind_front = find(ismember(Ch_labels,{'Fz','F3','F4'}));
        Chind_cent = find(ismember(Ch_labels,{'Cz','C3','C4'}));
        Chind_pari = find(ismember(Ch_labels,{'Pz','P3','P4'}));
        Chind_occip = find(ismember(Ch_labels,{'PO7','Oz','PO8'}));

        % calculate power for different frequencies
        % get absolute power per channel
        PowAbs_Ch = DATA(4).Abspow_chs;
            % Frontal region:
                F_Apow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(F_Apow,1)
                    F_Apow(ff,1) = mean(PowAbs_Ch(Chind_front,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff
            % Central region:
                C_Apow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(C_Apow,1)
                    C_Apow(ff,1) = mean(PowAbs_Ch(Chind_cent,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff
            % Parietal region:
                P_Apow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(P_Apow,1)
                    P_Apow(ff,1) = mean(PowAbs_Ch(Chind_pari,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff
            % Occipital region:
                O_Apow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(O_Apow,1)
                    O_Apow(ff,1) = mean(PowAbs_Ch(Chind_occip,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff

            Apow = [F_Apow; C_Apow; P_Apow; O_Apow];
            clear PowAbs_Ch F_Apow C_Apow P_Apow O_Apow
 
        % log power
        % get log power per channel
        PowLog_Ch = DATA(4).ft_EEG_freqdata.logpowspctrm_avg;
            % Frontal region:
                F_Lpow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(F_Lpow,1)
                    F_Lpow(ff,1) = mean(PowLog_Ch(Chind_front,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff
            % Central region:
                C_Lpow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(C_Lpow,1)
                    C_Lpow(ff,1) = mean(PowLog_Ch(Chind_cent,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff
            % Parietal region:
                P_Lpow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(P_Lpow,1)
                    P_Lpow(ff,1) = mean(PowLog_Ch(Chind_pari,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff
            % Occipital region:
                O_Lpow = nan(size(Freqs_inds,1),1);
                for ff = 1:size(O_Lpow,1)
                    O_Lpow(ff,1) = mean(PowLog_Ch(Chind_occip,Freqs_inds(ff,1):Freqs_inds(ff,2)),'all');
                end
                clear ff

            Lpow = [F_Lpow; C_Lpow; P_Lpow; O_Lpow];
            clear PowLog_Ch F_Lpow C_Lpow P_Lpow O_Lpow

            

        % relative power: sum(band of interest) / sum (1-45Hz) => not
        % including higher gamma
        PowAbs_Ch = DATA(4).Abspow_chs;
            % Frontal region:
                PowAbs_FrChs = mean(PowAbs_Ch(Chind_front,:),1);
                Deno = sum(PowAbs_FrChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                F_Rpow = nan(Freqs_N,1);
                for ff = 1:5 %exclude higher gamma at end
                    Pow_sum = sum(PowAbs_FrChs(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    F_Rpow(ff) = Pow_sum / Deno;
                    clear Pow_sum
                end
                clear ff PowAbs_FrChs Deno
            % Central region:
                PowAbs_CeChs = mean(PowAbs_Ch(Chind_cent,:),1);
                Deno = sum(PowAbs_CeChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                C_Rpow = nan(Freqs_N,1);
                for ff = 1:5 %exclude higher gamma at end
                    Pow_sum = sum(PowAbs_CeChs(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    C_Rpow(ff) = Pow_sum / Deno;
                    clear Pow_sum
                end
                clear ff PowAbs_CeChs Deno
            % Parietal region:
                PowAbs_PaChs = mean(PowAbs_Ch(Chind_pari,:),1);
                Deno = sum(PowAbs_PaChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                P_Rpow = nan(Freqs_N,1);
                for ff = 1:5 %exclude higher gamma at end
                    Pow_sum = sum(PowAbs_PaChs(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    P_Rpow(ff) = Pow_sum / Deno;
                    clear Pow_sum
                end
                clear ff PowAbs_PaChs Deno
            % Occipital region:
                PowAbs_OcChs = mean(PowAbs_Ch(Chind_occip,:),1);
                Deno = sum(PowAbs_OcChs(1, RelPow_inds(1,1):RelPow_inds(1,2)),2);
                O_Rpow = nan(Freqs_N,1);
                for ff = 1:5 %exclude higher gamma at end
                    Pow_sum = sum(PowAbs_OcChs(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    O_Rpow(ff) = Pow_sum / Deno;
                    clear Pow_sum
                end
                clear ff PowAbs_OcChs Deno

            Rpow = [F_Rpow; C_Rpow; P_Rpow; O_Rpow];
            clear PowAbs_Ch F_Rpow C_Rpow P_Rpow O_Rpow

            clear Ch_labels Chind_front Chind_cent Chind_pari Chind_occip

    else % no sufficient data for this condition
        % absolute power
        Apow = nan(Regions_N*Freqs_N,1);
        % log power
        Lpow = nan(Regions_N*Freqs_N,1);
        % relative power
        Rpow = nan(Regions_N*Freqs_N,1);
    end
    % Rename vars
    F_Apow = Apow; clear Apow
    F_Lpow = Lpow; clear Lpow
    F_Rpow = Rpow; clear Rpow



   

    % Create table with new data
    % Columns: IDfull, ID, session, condition (4), regions (4) frequency band (5), N trials,  absolute power, log power, relative power
    N_trials = [repmat(S_Neps, Regions_N*Freqs_N,1); repmat(T_Neps, Regions_N*Freqs_N,1); ...
        repmat(A_Neps, Regions_N*Freqs_N,1); repmat(F_Neps, Regions_N*Freqs_N,1)];
    % Power
    Abs_Power = [S_Apow; T_Apow; A_Apow; F_Apow];
    Log_Power = [S_Lpow; T_Lpow; A_Lpow; F_Lpow];
    Rel_Power = [S_Rpow; T_Rpow; A_Rpow; F_Rpow];
    % assemble table
    dataTable_currID = table(ClinicalID, IDfull, IDcode, Site, Sessions, Condition, Regions, Frequencies, N_trials, Abs_Power, Log_Power, Rel_Power);

    
    
    
    % add current ID data to table
    cd(MOI_location)
    if exist('RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.mat','file') == 2
        load RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.mat
        NewRowInds = [(ss-1)*height(dataTable_currID)+1 ss*height(dataTable_currID)];
        RCT1_Power_longformat(NewRowInds(1,1):NewRowInds(1,2),:) = dataTable_currID;
        clear NewRowInds
    else % first ID
        RCT1_Power_longformat = dataTable_currID;
    end

    % save the table
    save('RCT1_Power_CanBands_longformat_CxRxFB_upto45Hz.mat','RCT1_Power_longformat')
    
    % clear up
    clear ClinicalID ClinicalID1 IDfull IDcode Site Sessions N_trials Abs_Power Log_Power Rel_Power IDparts
    clear S_Apow T_Apow A_Apow F_Apow S_Lpow T_Lpow A_Lpow F_Lpow S_Rpow T_Rpow A_Rpow F_Rpow
    clear CurrID dataTable_currID S_Neps T_Neps A_Neps F_Neps PowCurr_abs
    clear RCT1_Power_longformat
    clear DATA

end


clear Freqs Freq_boundaries Freqs_inds Condition Frequency_bands Frequencies

%% Save table into csv format for Rstudio
cd(MOI_location)
load RCT1_Power_CanBands_longformat_CxRxFB.mat
% save as csv file
writetable(RCT1_Power_longformat, 'RCT1_Power_CanBands_longformat_CxRxFB.csv')
clear RCT1_Power_longformat 

