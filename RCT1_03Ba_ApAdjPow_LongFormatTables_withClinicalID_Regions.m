%% RCT_03Ba: Aperiodic-adjusted power for statistical analyses

% This script takes csv files from Python that fitted one-over-f models
% onto the absolute power data from script RCT1_02B. It takes the 1/f
% offset, 1/f slope, and Rsq and error fit of the 1/f model per condition
% per region, for the narrow range (3-28Hz). 
% Aperiodic-adjusted power is only extracted on power spectra based on 20 or more epochs
% and with a 1/f Rsq higher than .95 are included in the tables, otherwise
% the value is set to NaN. 

% The first section of the script extracts aperiodic-adjusted power for the
% low (6-12Hz) and high (15-25Hz) frequencies within the narrow fitting
% range (3-28Hz). 
% The second section of the script extracts aperiodic-adjusted power for
% the theta (6-7Hz), alpha (8-12Hz), and beta (15-25Hz) frequency bands. 

% Data are saved in a .csv file for further statistical analyses in
% Rstudio.

% Rianne Haartsen, PhD.; Feb - 2025
% Birkbeck University of London

%% Low vs high frequencies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Columns: CilinicalID, IDfull, IDcode for EEG, site, session, condition 
% (4), region (4), frequencies  (2), Rsq and Aperiodic Adjusted power for 
% 1/f fitted across narrow range
% frequency bands: 
% low: 6-12, high: 15-25 for narrow 1/f range
% 32 rows per ID (4*4*2)

Path_FOOOFspectra = 'xxx/DATA/01B_PreprocOOF/FOOOF_powerspectra_Region/';
FOOOF_Narrow = readtable(strcat(Path_FOOOFspectra, 'Narrow_FBoI/RCT1_Region_Narrow_data.csv'),'VariableNamingRule','preserve'); 

MOI_location = 'xxx/DATA/01B_PreprocOOF/';
cd(MOI_location)
load miniMADE_MeasuresOfInterest_report.mat

PSCcode_conversion = readtable('xxx/rct1_psc1_psc2.xlsx');



% Low vs high frequency bands
Freqs = miniMADE_MOI_table.Sfreqs{1};
LFreq_boundaries = [6 12];
LowFreqs_inds = [find(Freqs == LFreq_boundaries(1,1)), find(Freqs == LFreq_boundaries(1,2))]; 
nHFreq_boundaries = [15 25]; % for narrow range
nHighFreqs_inds = [find(Freqs == nHFreq_boundaries(1,1)), find(Freqs == nHFreq_boundaries(1,2))];

% Settings
Conditions_N = 4;
Regions_N = 4;
Frequencies_N = 2;

Condition = [repmat({'Social'}, Regions_N*Frequencies_N, 1); ...
    repmat({'Toy'}, Regions_N*Frequencies_N, 1);...
    repmat({'Abstract'}, Regions_N*Frequencies_N, 1);...
    repmat({'Fix'}, Regions_N*Frequencies_N, 1)];
Regions = repmat([repmat({'Frontal'}, Frequencies_N, 1); ...
    repmat({'Central'}, Frequencies_N, 1);...
    repmat({'Parietal'}, Frequencies_N, 1);...
    repmat({'Occipital'}, Frequencies_N, 1)],Conditions_N,1);
Frequencies = repmat({'Low';'High'},Conditions_N*Regions_N,1);


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
        % channel labels
        Ch_labels = DATA(1).ft_EEG_freqdata.label;
        Chind_front = find(ismember(Ch_labels,{'Fz','F3','F4'}));
        Chind_cent = find(ismember(Ch_labels,{'Cz','C3','C4'}));
        Chind_pari = find(ismember(Ch_labels,{'Pz','P3','P4'}));
        Chind_occip = find(ismember(Ch_labels,{'PO7','Oz','PO8'}));
        % get power for this condition and average across trials
        PowAbs_Ch = shiftdim(mean(DATA(1).ft_EEG_freqdata.powspctrm,1));

        % Frontal region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Soc_Frontal');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && S_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    fRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_front,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';

                    % figure; plot(Freqs, Aperiodic,'LineWidth',2)
                    % hold on
                    % plot(Freqs, Region_pow,'LineWidth',2)
                    % plot(Freqs, log(Region_pow),'LineWidth',2)
                    % plot(Freqs, log10(Region_pow),'LineWidth',2)
                    % legend({'1/f model','Absolute power (raw)','Log (base ln) power','Log (base 10) power'})
                    % xlim([3 28])


                    fApAdj_LFreq = mean(Region_logpow(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2) - mean(Aperiodic(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);    
                    fApAdj_HFreq = mean(Region_logpow(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2) - mean(Aperiodic(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2);    
                    clear Aperiodic Region_pow Region_logpow
                else
                    fRsq = FOOOF_Narrow.("r squared")(Nrow);
                    fApAdj_LFreq = NaN;
                    fApAdj_HFreq = NaN;
                end
            else
                fRsq = NaN;
                fApAdj_LFreq = NaN;
                fApAdj_HFreq = NaN;
            end
            clear Nrow Name Nr_offset Nr_slope


        % Central region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Soc_Central');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && S_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    cRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_cent,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    cApAdj_LFreq = mean(Region_logpow(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2) - mean(Aperiodic(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);    
                    cApAdj_HFreq = mean(Region_logpow(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2) - mean(Aperiodic(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2);    
                    clear Aperiodic Region_pow Region_logpow
                else
                    cRsq = FOOOF_Narrow.("r squared")(Nrow);
                    cApAdj_LFreq = NaN;
                    cApAdj_HFreq = NaN;
                end
            else
                cRsq = NaN;
                cApAdj_LFreq = NaN;
                cApAdj_HFreq = NaN;
            end
            clear Nrow Name Nr_offset Nr_slope

        % Parietal region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Soc_Parietal');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && S_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    pRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_pari,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    pApAdj_LFreq = mean(Region_logpow(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2) - mean(Aperiodic(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);    
                    pApAdj_HFreq = mean(Region_logpow(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2) - mean(Aperiodic(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2);    
                    clear Aperiodic Region_pow Region_logpow
                else
                    pRsq = FOOOF_Narrow.("r squared")(Nrow);
                    pApAdj_LFreq = NaN;
                    pApAdj_HFreq = NaN;
                end
            else
                pRsq = NaN;
                pApAdj_LFreq = NaN;
                pApAdj_HFreq = NaN;
            end
            clear Nrow Name Nr_offset Nr_slope


        % Occipital region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Soc_Occipital');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && S_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    oRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_occip,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    oApAdj_LFreq = mean(Region_logpow(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2) - mean(Aperiodic(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);    
                    oApAdj_HFreq = mean(Region_logpow(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2) - mean(Aperiodic(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2);    
                    clear Aperiodic Region_pow Region_logpow
                else
                    oRsq = FOOOF_Narrow.("r squared")(Nrow);
                    oApAdj_LFreq = NaN;
                    oApAdj_HFreq = NaN;
                end
            else
                oRsq = NaN;
                oApAdj_LFreq = NaN;
                oApAdj_HFreq = NaN;
            end
            clear Nrow Name Nr_offset Nr_slope

        % clear up
        clear Ch_labels Chind_front Chind_cent Chind_pari Chind_occip PowAbs_Ch 

    else % no sufficient data for this condition
        fRsq = NaN; cRsq = NaN; pRsq = NaN; oRsq = NaN;
        fApAdj_LFreq = NaN; cApAdj_LFreq = NaN; pApAdj_LFreq = NaN; oApAdj_LFreq = NaN;
        fApAdj_HFreq = NaN; cApAdj_HFreq = NaN; pApAdj_HFreq = NaN; oApAdj_HFreq = NaN;
    end
    % Rename vars
    S_nRsq = [fRsq; fRsq ; cRsq; cRsq; pRsq; pRsq; oRsq; oRsq];
    S_aapow = [fApAdj_LFreq; fApAdj_HFreq; cApAdj_LFreq; cApAdj_HFreq; pApAdj_LFreq; pApAdj_HFreq; oApAdj_LFreq; oApAdj_HFreq];
    clear fRsq cRsq pRsq oRsq 
    clear fApAdj_LFreq fApAdj_HFreq cApAdj_LFreq cApAdj_HFreq pApAdj_LFreq pApAdj_HFreq oApAdj_LFreq oApAdj_HFreq
    




    % Toy condition (condition x region x frequency band) %%%%%%%%%%%%%%%%%%%%%
    T_Neps = miniMADE_MOI_table.Teps(ss);
    if T_Neps >= 20 
        % channel labels
        Ch_labels = DATA(2).ft_EEG_freqdata.label;
        Chind_front = find(ismember(Ch_labels,{'Fz','F3','F4'}));
        Chind_cent = find(ismember(Ch_labels,{'Cz','C3','C4'}));
        Chind_pari = find(ismember(Ch_labels,{'Pz','P3','P4'}));
        Chind_occip = find(ismember(Ch_labels,{'PO7','Oz','PO8'}));
        % get power for this condition and average across trials
        PowAbs_Ch = shiftdim(mean(DATA(2).ft_EEG_freqdata.powspctrm,1));

        % Frontal region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Toy_Frontal');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && T_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    fRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_front,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    fApAdj_LFreq = mean(Region_logpow(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2) - mean(Aperiodic(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);    
                    fApAdj_HFreq = mean(Region_logpow(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2) - mean(Aperiodic(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2);    
                    clear Aperiodic Region_pow Region_logpow
                else
                    fRsq = FOOOF_Narrow.("r squared")(Nrow);
                    fApAdj_LFreq = NaN;
                    fApAdj_HFreq = NaN;
                end
            else
                fRsq = NaN;
                fApAdj_LFreq = NaN;
                fApAdj_HFreq = NaN;
            end
            clear Nrow Name Nr_offset Nr_slope


        % Central region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Toy_Central');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && T_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    cRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_cent,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    cApAdj_LFreq = mean(Region_logpow(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2) - mean(Aperiodic(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);    
                    cApAdj_HFreq = mean(Region_logpow(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2) - mean(Aperiodic(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2);    
                    clear Aperiodic Region_pow Region_logpow
                else
                    cRsq = FOOOF_Narrow.("r squared")(Nrow);
                    cApAdj_LFreq = NaN;
                    cApAdj_HFreq = NaN;
                end
            else
                cRsq = NaN;
                cApAdj_LFreq = NaN;
                cApAdj_HFreq = NaN;
            end
            clear Nrow Name Nr_offset Nr_slope

        % Parietal region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Toy_Parietal');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && T_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    pRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_pari,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    pApAdj_LFreq = mean(Region_logpow(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2) - mean(Aperiodic(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);    
                    pApAdj_HFreq = mean(Region_logpow(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2) - mean(Aperiodic(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2);    
                    clear Aperiodic Region_pow Region_logpow
                else
                    pRsq = FOOOF_Narrow.("r squared")(Nrow);
                    pApAdj_LFreq = NaN;
                    pApAdj_HFreq = NaN;
                end
            else
                pRsq = NaN;
                pApAdj_LFreq = NaN;
                pApAdj_HFreq = NaN;
            end
            clear Nrow Name Nr_offset Nr_slope


        % Occipital region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Toy_Occipital');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && T_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    oRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_occip,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    oApAdj_LFreq = mean(Region_logpow(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2) - mean(Aperiodic(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);    
                    oApAdj_HFreq = mean(Region_logpow(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2) - mean(Aperiodic(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2);    
                    clear Aperiodic Region_pow Region_logpow
                else
                    oRsq = FOOOF_Narrow.("r squared")(Nrow);
                    oApAdj_LFreq = NaN;
                    oApAdj_HFreq = NaN;
                end
            else
                oRsq = NaN;
                oApAdj_LFreq = NaN;
                oApAdj_HFreq = NaN;
            end
            clear Nrow Name Nr_offset Nr_slope

        % clear up
        clear Ch_labels Chind_front Chind_cent Chind_pari Chind_occip PowAbs_Ch 

    else % no sufficient data for this condition
        fRsq = NaN; cRsq = NaN; pRsq = NaN; oRsq = NaN;
        fApAdj_LFreq = NaN; cApAdj_LFreq = NaN; pApAdj_LFreq = NaN; oApAdj_LFreq = NaN;
        fApAdj_HFreq = NaN; cApAdj_HFreq = NaN; pApAdj_HFreq = NaN; oApAdj_HFreq = NaN;
    end
    % Rename vars
    T_nRsq = [fRsq; fRsq ; cRsq; cRsq; pRsq; pRsq; oRsq; oRsq];
    T_aapow = [fApAdj_LFreq; fApAdj_HFreq; cApAdj_LFreq; cApAdj_HFreq; pApAdj_LFreq; pApAdj_HFreq; oApAdj_LFreq; oApAdj_HFreq];
    clear fRsq cRsq pRsq oRsq 
    clear fApAdj_LFreq fApAdj_HFreq cApAdj_LFreq cApAdj_HFreq pApAdj_LFreq pApAdj_HFreq oApAdj_LFreq oApAdj_HFreq
    


% Abstract condition (condition x region x frequency band) %%%%%%%%%%%%%%%%
    A_Neps = miniMADE_MOI_table.Aeps(ss);
    if A_Neps >= 20 
        % channel labels
        Ch_labels = DATA(3).ft_EEG_freqdata.label;
        Chind_front = find(ismember(Ch_labels,{'Fz','F3','F4'}));
        Chind_cent = find(ismember(Ch_labels,{'Cz','C3','C4'}));
        Chind_pari = find(ismember(Ch_labels,{'Pz','P3','P4'}));
        Chind_occip = find(ismember(Ch_labels,{'PO7','Oz','PO8'}));
        % get power for this condition and average across trials
        PowAbs_Ch = shiftdim(mean(DATA(3).ft_EEG_freqdata.powspctrm,1));

        % Frontal region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Abs_Frontal');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && A_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    fRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_front,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    fApAdj_LFreq = mean(Region_logpow(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2) - mean(Aperiodic(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);    
                    fApAdj_HFreq = mean(Region_logpow(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2) - mean(Aperiodic(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2);    
                    clear Aperiodic Region_pow Region_logpow
                else
                    fRsq = FOOOF_Narrow.("r squared")(Nrow);
                    fApAdj_LFreq = NaN;
                    fApAdj_HFreq = NaN;
                end
            else
                fRsq = NaN;
                fApAdj_LFreq = NaN;
                fApAdj_HFreq = NaN;
            end
            clear Nrow Name Nr_offset Nr_slope


        % Central region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Abs_Central');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && A_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    cRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_cent,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    cApAdj_LFreq = mean(Region_logpow(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2) - mean(Aperiodic(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);    
                    cApAdj_HFreq = mean(Region_logpow(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2) - mean(Aperiodic(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2);    
                    clear Aperiodic Region_pow Region_logpow
                else
                    cRsq = FOOOF_Narrow.("r squared")(Nrow);
                    cApAdj_LFreq = NaN;
                    cApAdj_HFreq = NaN;
                end
            else
                cRsq = NaN;
                cApAdj_LFreq = NaN;
                cApAdj_HFreq = NaN;
            end
            clear Nrow Name Nr_offset Nr_slope

        % Parietal region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Abs_Parietal');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && A_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    pRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_pari,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    pApAdj_LFreq = mean(Region_logpow(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2) - mean(Aperiodic(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);    
                    pApAdj_HFreq = mean(Region_logpow(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2) - mean(Aperiodic(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2);    
                    clear Aperiodic Region_pow Region_logpow
                else
                    pRsq = FOOOF_Narrow.("r squared")(Nrow);
                    pApAdj_LFreq = NaN;
                    pApAdj_HFreq = NaN;
                end
            else
                pRsq = NaN;
                pApAdj_LFreq = NaN;
                pApAdj_HFreq = NaN;
            end
            clear Nrow Name Nr_offset Nr_slope


        % Occipital region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Abs_Occipital');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && A_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    oRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_occip,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    oApAdj_LFreq = mean(Region_logpow(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2) - mean(Aperiodic(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);    
                    oApAdj_HFreq = mean(Region_logpow(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2) - mean(Aperiodic(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2);    
                    clear Aperiodic Region_pow Region_logpow
                else
                    oRsq = FOOOF_Narrow.("r squared")(Nrow);
                    oApAdj_LFreq = NaN;
                    oApAdj_HFreq = NaN;
                end
            else
                oRsq = NaN;
                oApAdj_LFreq = NaN;
                oApAdj_HFreq = NaN;
            end
            clear Nrow Name Nr_offset Nr_slope

        % clear up
        clear Ch_labels Chind_front Chind_cent Chind_pari Chind_occip PowAbs_Ch 

    else % no sufficient data for this condition
        fRsq = NaN; cRsq = NaN; pRsq = NaN; oRsq = NaN;
        fApAdj_LFreq = NaN; cApAdj_LFreq = NaN; pApAdj_LFreq = NaN; oApAdj_LFreq = NaN;
        fApAdj_HFreq = NaN; cApAdj_HFreq = NaN; pApAdj_HFreq = NaN; oApAdj_HFreq = NaN;
    end
    % Rename vars
    A_nRsq = [fRsq; fRsq ; cRsq; cRsq; pRsq; pRsq; oRsq; oRsq];
    A_aapow = [fApAdj_LFreq; fApAdj_HFreq; cApAdj_LFreq; cApAdj_HFreq; pApAdj_LFreq; pApAdj_HFreq; oApAdj_LFreq; oApAdj_HFreq];
    clear fRsq cRsq pRsq oRsq 
    clear fApAdj_LFreq fApAdj_HFreq cApAdj_LFreq cApAdj_HFreq pApAdj_LFreq pApAdj_HFreq oApAdj_LFreq oApAdj_HFreq
    


    % Fixation condition (condition x region x frequency band) %%%%%%%%%%%%%%%%
    F_Neps = miniMADE_MOI_table.Feps(ss);
    if F_Neps >= 20 
        % channel labels
        Ch_labels = DATA(4).ft_EEG_freqdata.label;
        Chind_front = find(ismember(Ch_labels,{'Fz','F3','F4'}));
        Chind_cent = find(ismember(Ch_labels,{'Cz','C3','C4'}));
        Chind_pari = find(ismember(Ch_labels,{'Pz','P3','P4'}));
        Chind_occip = find(ismember(Ch_labels,{'PO7','Oz','PO8'}));
        % get power for this condition and average across trials
        PowAbs_Ch = shiftdim(mean(DATA(4).ft_EEG_freqdata.powspctrm,1));

        % Frontal region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Fix_Frontal');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && F_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    fRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_front,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    fApAdj_LFreq = mean(Region_logpow(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2) - mean(Aperiodic(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);    
                    fApAdj_HFreq = mean(Region_logpow(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2) - mean(Aperiodic(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2);    
                    clear Aperiodic Region_pow Region_logpow
                else
                    fRsq = FOOOF_Narrow.("r squared")(Nrow);
                    fApAdj_LFreq = NaN;
                    fApAdj_HFreq = NaN;
                end
            else
                fRsq = NaN;
                fApAdj_LFreq = NaN;
                fApAdj_HFreq = NaN;
            end
            clear Nrow Name Nr_offset Nr_slope


        % Central region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Fix_Central');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && F_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    cRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_cent,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    cApAdj_LFreq = mean(Region_logpow(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2) - mean(Aperiodic(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);    
                    cApAdj_HFreq = mean(Region_logpow(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2) - mean(Aperiodic(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2);    
                    clear Aperiodic Region_pow Region_logpow
                else
                    cRsq = FOOOF_Narrow.("r squared")(Nrow);
                    cApAdj_LFreq = NaN;
                    cApAdj_HFreq = NaN;
                end
            else
                cRsq = NaN;
                cApAdj_LFreq = NaN;
                cApAdj_HFreq = NaN;
            end
            clear Nrow Name Nr_offset Nr_slope

        % Parietal region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Fix_Parietal');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && F_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    pRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_pari,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    pApAdj_LFreq = mean(Region_logpow(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2) - mean(Aperiodic(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);    
                    pApAdj_HFreq = mean(Region_logpow(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2) - mean(Aperiodic(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2);    
                    clear Aperiodic Region_pow Region_logpow
                else
                    pRsq = FOOOF_Narrow.("r squared")(Nrow);
                    pApAdj_LFreq = NaN;
                    pApAdj_HFreq = NaN;
                end
            else
                pRsq = NaN;
                pApAdj_LFreq = NaN;
                pApAdj_HFreq = NaN;
            end
            clear Nrow Name Nr_offset Nr_slope


        % Occipital region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Fix_Occipital');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && F_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    oRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_occip,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    oApAdj_LFreq = mean(Region_logpow(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2) - mean(Aperiodic(1,LowFreqs_inds(1,1):LowFreqs_inds(1,2)),2);    
                    oApAdj_HFreq = mean(Region_logpow(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2) - mean(Aperiodic(1,nHighFreqs_inds(1,1):nHighFreqs_inds(1,2)),2);    
                    clear Aperiodic Region_pow Region_logpow
                else
                    oRsq = FOOOF_Narrow.("r squared")(Nrow);
                    oApAdj_LFreq = NaN;
                    oApAdj_HFreq = NaN;
                end
            else
                oRsq = NaN;
                oApAdj_LFreq = NaN;
                oApAdj_HFreq = NaN;
            end
            clear Nrow Name Nr_offset Nr_slope

        % clear up
        clear Ch_labels Chind_front Chind_cent Chind_pari Chind_occip PowAbs_Ch 

    else % no sufficient data for this condition
        fRsq = NaN; cRsq = NaN; pRsq = NaN; oRsq = NaN;
        fApAdj_LFreq = NaN; cApAdj_LFreq = NaN; pApAdj_LFreq = NaN; oApAdj_LFreq = NaN;
        fApAdj_HFreq = NaN; cApAdj_HFreq = NaN; pApAdj_HFreq = NaN; oApAdj_HFreq = NaN;
    end
    % Rename vars
    F_nRsq = [fRsq; fRsq ; cRsq; cRsq; pRsq; pRsq; oRsq; oRsq];
    F_aapow = [fApAdj_LFreq; fApAdj_HFreq; cApAdj_LFreq; cApAdj_HFreq; pApAdj_LFreq; pApAdj_HFreq; oApAdj_LFreq; oApAdj_HFreq];
    clear fRsq cRsq pRsq oRsq 
    clear fApAdj_LFreq fApAdj_HFreq cApAdj_LFreq cApAdj_HFreq pApAdj_LFreq pApAdj_HFreq oApAdj_LFreq oApAdj_HFreq
    



    % Create table with new data
    % Columns: IDfull, ID, session, condition (4), region (4), frequency
    % band (2), N trials, aperiodic adjusted power
    % Ntrials
    N_trials = [repmat(S_Neps, Regions_N*Frequencies_N,1); repmat(T_Neps, Regions_N*Frequencies_N,1); ...
        repmat(A_Neps, Regions_N*Frequencies_N,1); repmat(F_Neps, Regions_N*Frequencies_N,1)];
    % Power
    Nr_Rsq = [S_nRsq; T_nRsq; A_nRsq; F_nRsq];
    Nr_Aapow = [S_aapow; T_aapow; A_aapow; F_aapow];
    % assemble table
    dataTable_currID = table(ClinicalID, IDfull, IDcode, Site, Sessions, Condition, Regions, Frequencies, N_trials, Nr_Rsq, Nr_Aapow);
    
    % add current ID data to table
    cd(MOI_location)
    if exist('RCT1_ApAdjPower_LvsHfreqs_longformat_CxRxFB.mat','file') == 2
        load RCT1_ApAdjPower_LvsHfreqs_longformat_CxRxFB.mat
        NewRowInds = [(ss-1)*height(dataTable_currID)+1 ss*height(dataTable_currID)];
        RCT1_ApAdjPower_longformat(NewRowInds(1,1):NewRowInds(1,2),:) = dataTable_currID;
        clear NewRowInds
    else % first ID
        RCT1_ApAdjPower_longformat = dataTable_currID;
    end

    % save the table
    save('RCT1_ApAdjPower_LvsHfreqs_longformat_CxRxFB.mat','RCT1_ApAdjPower_longformat')
    
    % clear up
    clear ClinicalID ClinicalID1 IDfull IDparts IDcode Site Sessions N_trials 
    clear Nr_Rsq S_nRsq T_nRsq A_nRsq F_nRsq Nr_Aapow S_aapow T_aapow A_aapow F_aapow
    clear CurrID dataTable_currID S_Neps T_Neps A_Neps F_Neps 
    clear RCT1_ApAdjPower_longformat
    clear DATA

end

% clear up all
clear Freqs Frequencies_N Freqs_inds Condition Conditions_N Frequency_bands Frequencies Regions Regions_N
clear LFreq_boundaries LowFreqs_inds nHFreq_boundaries nHighFreqs_inds ss 
clear FOOOF_Narrow miniMADE_MOI_table MOI_location  Path_FOOOFspectra PSCcode_conversion


% save as csv file
load("RCT1_ApAdjPower_LvsHfreqs_longformat_CxRxFB.mat")
writetable(RCT1_ApAdjPower_longformat, 'RCT1_ApAdjPower_LvsHfreqs_longformat_CxRxFB.csv')
clear RCT1_ApAdjPower_longformat









%% Aperiodic adjusted power for canonical frequency bands
% based on the OOF data


% Columns: ClinicalID, IDfull, IDcode for EEG, site, session, condition 
% (4), region (4), frequencies  (4), Rsq and Aperiodic Adjusted power for 
% 1/f fitted across narrow range
% frequency bands: 
% delta band: 3-5Hz, theta band: 6-7Hz, alpha band: 8-12Hz, beta band: 15-25Hz
% 48 rows per ID (4*4*3)


Path_FOOOFspectra = 'xxx/DATA/01B_PreprocOOF/FOOOF_powerspectra_Region/';
FOOOF_Narrow = readtable(strcat(Path_FOOOFspectra, 'Narrow_FBoI/RCT1_Region_Narrow_data.csv'),'VariableNamingRule','preserve'); 

MOI_location = 'xxx/DATA/01B_PreprocOOF/';
cd(MOI_location)
load miniMADE_MeasuresOfInterest_report.mat

PSCcode_conversion = readtable('xxx/rct1_psc1_psc2.xlsx');

% Settings
Conditions_N = 4;
Regions_N = 4;
Frequencies_N = 3;

Freqs = miniMADE_MOI_table.Sfreqs{1};
Freq_boundaries = [6 7; 8 12; 15 25];
Freqs_inds = [find(Freqs == Freq_boundaries(1,1)), find(Freqs == Freq_boundaries(1,2)); ...
    find(Freqs == Freq_boundaries(2,1)), find(Freqs == Freq_boundaries(2,2)); ...
    find(Freqs == Freq_boundaries(3,1)), find(Freqs == Freq_boundaries(3,2))];
Condition = [repmat({'Social'}, Regions_N*Frequencies_N, 1); ...
    repmat({'Toy'}, Regions_N*Frequencies_N, 1);...
    repmat({'Abstract'}, Regions_N*Frequencies_N, 1);...
    repmat({'Fix'}, Regions_N*Frequencies_N, 1)];
Regions = repmat([repmat({'Frontal'}, Frequencies_N, 1); ...
    repmat({'Central'}, Frequencies_N, 1);...
    repmat({'Parietal'}, Frequencies_N, 1);...
    repmat({'Occipital'}, Frequencies_N, 1)],Conditions_N,1);
Frequencies = repmat({'Theta';'Alpha';'Beta'},Conditions_N*Regions_N,1);


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
        % channel labels
        Ch_labels = DATA(1).ft_EEG_freqdata.label;
        Chind_front = find(ismember(Ch_labels,{'Fz','F3','F4'}));
        Chind_cent = find(ismember(Ch_labels,{'Cz','C3','C4'}));
        Chind_pari = find(ismember(Ch_labels,{'Pz','P3','P4'}));
        Chind_occip = find(ismember(Ch_labels,{'PO7','Oz','PO8'}));
        % get power for this condition and average across trials
        PowAbs_Ch = shiftdim(mean(DATA(1).ft_EEG_freqdata.powspctrm,1));

        % Frontal region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Soc_Frontal');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && S_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    fRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_front,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';

                    % figure; plot(Freqs, Aperiodic,'LineWidth',2)
                    % hold on
                    % plot(Freqs, Region_pow,'LineWidth',2)
                    % plot(Freqs, log(Region_pow),'LineWidth',2)
                    % plot(Freqs, log10(Region_pow),'LineWidth',2)
                    % legend({'1/f model','Absolute power (raw)','Log (base ln) power','Log (base 10) power'})
                    % xlim([1 28])
                    % xline([2 3 5 6 7 8 12 15 25]); 
                    % grid on

                    F_ApAdjpow = zeros(size(Freqs_inds,1),1);
                    for ff = 1:size(Freqs_inds,1)
                        F_ApAdjpow(ff) = mean(Region_logpow(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2) - mean(Aperiodic(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    end
                    clear PowCurr Aperiodic Region_logpow Region_pow ff
                else
                    fRsq = FOOOF_Narrow.("r squared")(Nrow); F_ApAdjpow = nan(size(Freqs_inds,1),1);
                end
            else
                fRsq = NaN; F_ApAdjpow = nan(size(Freqs_inds,1),1);
            end
            clear Nrow Name Nr_offset Nr_slope
    
        % Central region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Soc_Central');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && S_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    cRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_cent,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    C_ApAdjpow = zeros(size(Freqs_inds,1),1);
                    for ff = 1:size(Freqs_inds,1)
                        C_ApAdjpow(ff) = mean(Region_logpow(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2) - mean(Aperiodic(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    end
                    clear PowCurr Aperiodic Region_logpow Region_pow ff
                else
                    cRsq = FOOOF_Narrow.("r squared")(Nrow); C_ApAdjpow = nan(size(Freqs_inds,1),1);
                end
            else
                cRsq = NaN; C_ApAdjpow = nan(size(Freqs_inds,1),1);
            end
            clear Nrow Name Nr_offset Nr_slope

        % Parietal region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Soc_Parietal');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && S_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    pRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_pari,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    P_ApAdjpow = zeros(size(Freqs_inds,1),1);
                    for ff = 1:size(Freqs_inds,1)
                        P_ApAdjpow(ff) = mean(Region_logpow(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2) - mean(Aperiodic(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    end
                    clear PowCurr Aperiodic Region_logpow Region_pow ff
                else
                    pRsq = FOOOF_Narrow.("r squared")(Nrow); P_ApAdjpow = nan(size(Freqs_inds,1),1);
                end
            else
                pRsq = NaN; P_ApAdjpow = nan(size(Freqs_inds,1),1);
            end
            clear Nrow Name Nr_offset Nr_slope


        % Occipital region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Soc_Occipital');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && S_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    oRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_occip,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    O_ApAdjpow = zeros(size(Freqs_inds,1),1);
                    for ff = 1:size(Freqs_inds,1)
                        O_ApAdjpow(ff) = mean(Region_logpow(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2) - mean(Aperiodic(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    end
                    clear PowCurr Aperiodic Region_logpow Region_pow ff
                else
                    oRsq = FOOOF_Narrow.("r squared")(Nrow); O_ApAdjpow = nan(size(Freqs_inds,1),1);
                end
            else
                oRsq = NaN; O_ApAdjpow = nan(size(Freqs_inds,1),1);
            end
            clear Nrow Name Nr_offset Nr_slope


        % clear up
        clear Ch_labels Chind_front Chind_cent Chind_pari Chind_occip PowAbs_Ch 

    else % no sufficient data for this condition
        fRsq = NaN; cRsq = NaN; pRsq = NaN; oRsq = NaN;
        F_ApAdjpow = nan(size(Freqs_inds,1),1);
        C_ApAdjpow = nan(size(Freqs_inds,1),1);
        P_ApAdjpow = nan(size(Freqs_inds,1),1);
        O_ApAdjpow = nan(size(Freqs_inds,1),1);
    end
    % Rename vars
    S_nRsq = [repmat(fRsq, size(F_ApAdjpow,1), 1); repmat(cRsq, size(C_ApAdjpow,1), 1); repmat(pRsq, size(P_ApAdjpow,1), 1); repmat(oRsq, size(O_ApAdjpow,1), 1)];
    S_aapow = [F_ApAdjpow; C_ApAdjpow; P_ApAdjpow; O_ApAdjpow];
    clear fRsq cRsq pRsq oRsq 
    clear F_ApAdjpow C_ApAdjpow P_ApAdjpow O_ApAdjpow
    

    % Toy condition (condition x region x frequency band) %%%%%%%%%%%%%%
    T_Neps = miniMADE_MOI_table.Teps(ss);
    if T_Neps >= 20 
        % channel labels
        Ch_labels = DATA(2).ft_EEG_freqdata.label;
        Chind_front = find(ismember(Ch_labels,{'Fz','F3','F4'}));
        Chind_cent = find(ismember(Ch_labels,{'Cz','C3','C4'}));
        Chind_pari = find(ismember(Ch_labels,{'Pz','P3','P4'}));
        Chind_occip = find(ismember(Ch_labels,{'PO7','Oz','PO8'}));
        % get power for this condition and average across trials
        PowAbs_Ch = shiftdim(mean(DATA(2).ft_EEG_freqdata.powspctrm,1));

        % Frontal region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Toy_Frontal');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && S_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    fRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_front,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    F_ApAdjpow = zeros(size(Freqs_inds,1),1);
                    for ff = 1:size(Freqs_inds,1)
                        F_ApAdjpow(ff) = mean(Region_logpow(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2) - mean(Aperiodic(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    end
                    clear PowCurr Aperiodic Region_logpow Region_pow ff
                else
                    fRsq = FOOOF_Narrow.("r squared")(Nrow); F_ApAdjpow = nan(size(Freqs_inds,1),1);
                end
            else
                fRsq = NaN; F_ApAdjpow = nan(size(Freqs_inds,1),1);
            end
            clear Nrow Name Nr_offset Nr_slope
    
        % Central region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Toy_Central');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && S_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    cRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_cent,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    C_ApAdjpow = zeros(size(Freqs_inds,1),1);
                    for ff = 1:size(Freqs_inds,1)
                        C_ApAdjpow(ff) = mean(Region_logpow(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2) - mean(Aperiodic(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    end
                    clear PowCurr Aperiodic Region_logpow Region_pow ff
                else
                    cRsq = FOOOF_Narrow.("r squared")(Nrow); C_ApAdjpow = nan(size(Freqs_inds,1),1);
                end
            else
                cRsq = NaN; C_ApAdjpow = nan(size(Freqs_inds,1),1);
            end
            clear Nrow Name Nr_offset Nr_slope

        % Parietal region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Toy_Parietal');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && S_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    pRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_pari,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    P_ApAdjpow = zeros(size(Freqs_inds,1),1);
                    for ff = 1:size(Freqs_inds,1)
                        P_ApAdjpow(ff) = mean(Region_logpow(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2) - mean(Aperiodic(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    end
                    clear PowCurr Aperiodic Region_logpow Region_pow ff
                else
                    pRsq = FOOOF_Narrow.("r squared")(Nrow); P_ApAdjpow = nan(size(Freqs_inds,1),1);
                end
            else
                pRsq = NaN; P_ApAdjpow = nan(size(Freqs_inds,1),1);
            end
            clear Nrow Name Nr_offset Nr_slope


        % Occipital region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Toy_Occipital');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && S_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    oRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_occip,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    O_ApAdjpow = zeros(size(Freqs_inds,1),1);
                    for ff = 1:size(Freqs_inds,1)
                        O_ApAdjpow(ff) = mean(Region_logpow(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2) - mean(Aperiodic(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    end
                    clear PowCurr Aperiodic Region_logpow Region_pow ff
                else
                    oRsq = FOOOF_Narrow.("r squared")(Nrow); O_ApAdjpow = nan(size(Freqs_inds,1),1);
                end
            else
                oRsq = NaN; O_ApAdjpow = nan(size(Freqs_inds,1),1);
            end
            clear Nrow Name Nr_offset Nr_slope


        % clear up
        clear Ch_labels Chind_front Chind_cent Chind_pari Chind_occip PowAbs_Ch 

    else % no sufficient data for this condition
        fRsq = NaN; cRsq = NaN; pRsq = NaN; oRsq = NaN;
        F_ApAdjpow = nan(size(Freqs_inds,1),1);
        C_ApAdjpow = nan(size(Freqs_inds,1),1);
        P_ApAdjpow = nan(size(Freqs_inds,1),1);
        O_ApAdjpow = nan(size(Freqs_inds,1),1);
    end
    % Rename vars
    T_nRsq = [repmat(fRsq, size(F_ApAdjpow,1), 1); repmat(cRsq, size(C_ApAdjpow,1), 1); repmat(pRsq, size(P_ApAdjpow,1), 1); repmat(oRsq, size(O_ApAdjpow,1), 1)];
    T_aapow = [F_ApAdjpow; C_ApAdjpow; P_ApAdjpow; O_ApAdjpow];
    clear fRsq cRsq pRsq oRsq 
    clear F_ApAdjpow C_ApAdjpow P_ApAdjpow O_ApAdjpow
    


    % Abstract condition (condition x region x frequency band) %%%%%%%%%%%%%%
    A_Neps = miniMADE_MOI_table.Aeps(ss);
    if A_Neps >= 20 
        % channel labels
        Ch_labels = DATA(3).ft_EEG_freqdata.label;
        Chind_front = find(ismember(Ch_labels,{'Fz','F3','F4'}));
        Chind_cent = find(ismember(Ch_labels,{'Cz','C3','C4'}));
        Chind_pari = find(ismember(Ch_labels,{'Pz','P3','P4'}));
        Chind_occip = find(ismember(Ch_labels,{'PO7','Oz','PO8'}));
        % get power for this condition and average across trials
        PowAbs_Ch = shiftdim(mean(DATA(3).ft_EEG_freqdata.powspctrm,1));

        % Frontal region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Abs_Frontal');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && S_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    fRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_front,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    F_ApAdjpow = zeros(size(Freqs_inds,1),1);
                    for ff = 1:size(Freqs_inds,1)
                        F_ApAdjpow(ff) = mean(Region_logpow(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2) - mean(Aperiodic(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    end
                    clear PowCurr Aperiodic Region_logpow Region_pow ff
                else
                    fRsq = FOOOF_Narrow.("r squared")(Nrow); F_ApAdjpow = nan(size(Freqs_inds,1),1);
                end
            else
                fRsq = NaN; F_ApAdjpow = nan(size(Freqs_inds,1),1);
            end
            clear Nrow Name Nr_offset Nr_slope
    
        % Central region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Abs_Central');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && S_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    cRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_cent,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    C_ApAdjpow = zeros(size(Freqs_inds,1),1);
                    for ff = 1:size(Freqs_inds,1)
                        C_ApAdjpow(ff) = mean(Region_logpow(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2) - mean(Aperiodic(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    end
                    clear PowCurr Aperiodic Region_logpow Region_pow ff
                else
                    cRsq = FOOOF_Narrow.("r squared")(Nrow); C_ApAdjpow = nan(size(Freqs_inds,1),1);
                end
            else
                cRsq = NaN; C_ApAdjpow = nan(size(Freqs_inds,1),1);
            end
            clear Nrow Name Nr_offset Nr_slope

        % Parietal region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Abs_Parietal');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && S_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    pRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_pari,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    P_ApAdjpow = zeros(size(Freqs_inds,1),1);
                    for ff = 1:size(Freqs_inds,1)
                        P_ApAdjpow(ff) = mean(Region_logpow(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2) - mean(Aperiodic(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    end
                    clear PowCurr Aperiodic Region_logpow Region_pow ff
                else
                    pRsq = FOOOF_Narrow.("r squared")(Nrow); P_ApAdjpow = nan(size(Freqs_inds,1),1);
                end
            else
                pRsq = NaN; P_ApAdjpow = nan(size(Freqs_inds,1),1);
            end
            clear Nrow Name Nr_offset Nr_slope


        % Occipital region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Abs_Occipital');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && S_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    oRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_occip,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    O_ApAdjpow = zeros(size(Freqs_inds,1),1);
                    for ff = 1:size(Freqs_inds,1)
                        O_ApAdjpow(ff) = mean(Region_logpow(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2) - mean(Aperiodic(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    end
                    clear PowCurr Aperiodic Region_logpow Region_pow ff
                else
                    oRsq = FOOOF_Narrow.("r squared")(Nrow); O_ApAdjpow = nan(size(Freqs_inds,1),1);
                end
            else
                oRsq = NaN; O_ApAdjpow = nan(size(Freqs_inds,1),1);
            end
            clear Nrow Name Nr_offset Nr_slope


        % clear up
        clear Ch_labels Chind_front Chind_cent Chind_pari Chind_occip PowAbs_Ch 

    else % no sufficient data for this condition
        fRsq = NaN; cRsq = NaN; pRsq = NaN; oRsq = NaN;
        F_ApAdjpow = nan(size(Freqs_inds,1),1);
        C_ApAdjpow = nan(size(Freqs_inds,1),1);
        P_ApAdjpow = nan(size(Freqs_inds,1),1);
        O_ApAdjpow = nan(size(Freqs_inds,1),1);
    end
    % Rename vars
    A_nRsq = [repmat(fRsq, size(F_ApAdjpow,1), 1); repmat(cRsq, size(C_ApAdjpow,1), 1); repmat(pRsq, size(P_ApAdjpow,1), 1); repmat(oRsq, size(O_ApAdjpow,1), 1)];
    A_aapow = [F_ApAdjpow; C_ApAdjpow; P_ApAdjpow; O_ApAdjpow];
    clear fRsq cRsq pRsq oRsq 
    clear F_ApAdjpow C_ApAdjpow P_ApAdjpow O_ApAdjpow

    
    % Fixation condition (condition x region x frequency band) %%%%%%%%%%%%%%
    F_Neps = miniMADE_MOI_table.Feps(ss);
    if F_Neps >= 20 
        % channel labels
        Ch_labels = DATA(4).ft_EEG_freqdata.label;
        Chind_front = find(ismember(Ch_labels,{'Fz','F3','F4'}));
        Chind_cent = find(ismember(Ch_labels,{'Cz','C3','C4'}));
        Chind_pari = find(ismember(Ch_labels,{'Pz','P3','P4'}));
        Chind_occip = find(ismember(Ch_labels,{'PO7','Oz','PO8'}));
        % get power for this condition and average across trials
        PowAbs_Ch = shiftdim(mean(DATA(4).ft_EEG_freqdata.powspctrm,1));

        % Frontal region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Fix_Frontal');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && S_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    fRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_front,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    F_ApAdjpow = zeros(size(Freqs_inds,1),1);
                    for ff = 1:size(Freqs_inds,1)
                        F_ApAdjpow(ff) = mean(Region_logpow(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2) - mean(Aperiodic(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    end
                    clear PowCurr Aperiodic Region_logpow Region_pow ff
                else
                    fRsq = FOOOF_Narrow.("r squared")(Nrow); F_ApAdjpow = nan(size(Freqs_inds,1),1);
                end
            else
                fRsq = NaN; F_ApAdjpow = nan(size(Freqs_inds,1),1);
            end
            clear Nrow Name Nr_offset Nr_slope
    
        % Central region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Fix_Central');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && S_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    cRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_cent,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    C_ApAdjpow = zeros(size(Freqs_inds,1),1);
                    for ff = 1:size(Freqs_inds,1)
                        C_ApAdjpow(ff) = mean(Region_logpow(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2) - mean(Aperiodic(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    end
                    clear PowCurr Aperiodic Region_logpow Region_pow ff
                else
                    cRsq = FOOOF_Narrow.("r squared")(Nrow); C_ApAdjpow = nan(size(Freqs_inds,1),1);
                end
            else
                cRsq = NaN; C_ApAdjpow = nan(size(Freqs_inds,1),1);
            end
            clear Nrow Name Nr_offset Nr_slope

        % Parietal region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Fix_Parietal');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && S_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    pRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_pari,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    P_ApAdjpow = zeros(size(Freqs_inds,1),1);
                    for ff = 1:size(Freqs_inds,1)
                        P_ApAdjpow(ff) = mean(Region_logpow(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2) - mean(Aperiodic(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    end
                    clear PowCurr Aperiodic Region_logpow Region_pow ff
                else
                    pRsq = FOOOF_Narrow.("r squared")(Nrow); P_ApAdjpow = nan(size(Freqs_inds,1),1);
                end
            else
                pRsq = NaN; P_ApAdjpow = nan(size(Freqs_inds,1),1);
            end
            clear Nrow Name Nr_offset Nr_slope


        % Occipital region:
            % find 1/f parameters and calculate 1/f spectrum
            Name = strcat(CurrID, '_Fix_Occipital');
            Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
            if ~isempty(Nrow) && S_Neps >= 20
                if FOOOF_Narrow.("r squared")(Nrow) > .95
                    oRsq = FOOOF_Narrow.("r squared")(Nrow);
                    Nr_offset = FOOOF_Narrow.intercept(Nrow);
                    Nr_slope = FOOOF_Narrow.slope(Nrow);
                    % get power spectrum for this region and condition
                    Region_pow = mean(PowAbs_Ch(Chind_occip,:),1);
                    Region_logpow = log10(Region_pow); % apply log base 10 to data as in FOOOF
                    % calculate power for different frequencies in the
                    % narrow range
                    Aperiodic = zeros(length(Freqs),1);
                    for ii = 1:length(Aperiodic)
                        freq_cur = Freqs(1,ii);
                        Aperiodic(ii,1) = Nr_offset - log10(0+freq_cur^(Nr_slope));
                    end
                    clear ii freq_cur
                    Aperiodic = Aperiodic';
                    % calculate aperiodic adjusted power
                    O_ApAdjpow = zeros(size(Freqs_inds,1),1);
                    for ff = 1:size(Freqs_inds,1)
                        O_ApAdjpow(ff) = mean(Region_logpow(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2) - mean(Aperiodic(1,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
                    end
                    clear PowCurr Aperiodic Region_logpow Region_pow ff
                else
                    oRsq = FOOOF_Narrow.("r squared")(Nrow); O_ApAdjpow = nan(size(Freqs_inds,1),1);
                end
            else
                oRsq = NaN; O_ApAdjpow = nan(size(Freqs_inds,1),1);
            end
            clear Nrow Name Nr_offset Nr_slope


        % clear up
        clear Ch_labels Chind_front Chind_cent Chind_pari Chind_occip PowAbs_Ch 

    else % no sufficient data for this condition
        fRsq = NaN; cRsq = NaN; pRsq = NaN; oRsq = NaN;
        F_ApAdjpow = nan(size(Freqs_inds,1),1);
        C_ApAdjpow = nan(size(Freqs_inds,1),1);
        P_ApAdjpow = nan(size(Freqs_inds,1),1);
        O_ApAdjpow = nan(size(Freqs_inds,1),1);
    end
    % Rename vars
    F_nRsq = [repmat(fRsq, size(F_ApAdjpow,1), 1); repmat(cRsq, size(C_ApAdjpow,1), 1); repmat(pRsq, size(P_ApAdjpow,1), 1); repmat(oRsq, size(O_ApAdjpow,1), 1)];
    F_aapow = [F_ApAdjpow; C_ApAdjpow; P_ApAdjpow; O_ApAdjpow];
    clear fRsq cRsq pRsq oRsq 
    clear F_ApAdjpow C_ApAdjpow P_ApAdjpow O_ApAdjpow







    % Create table with new data
    % Columns: IDfull, ID, session, condition (4), region (4), frequency
    % band (3), N trials, aperiodic adjusted power
    % Ntrials
    N_trials = [repmat(S_Neps, Regions_N*Frequencies_N,1); repmat(T_Neps, Regions_N*Frequencies_N,1); ...
        repmat(A_Neps, Regions_N*Frequencies_N,1); repmat(F_Neps, Regions_N*Frequencies_N,1)];
    % Power
    Nr_Rsq = [S_nRsq; T_nRsq; A_nRsq; F_nRsq];
    Nr_Aapow = [S_aapow; T_aapow; A_aapow; F_aapow];
    % assemble table
    dataTable_currID = table(ClinicalID, IDfull, IDcode, Site, Sessions, Condition, Regions, Frequencies, N_trials, Nr_Rsq, Nr_Aapow);
    
    % add current ID data to table
    cd(MOI_location)
    if exist('RCT1_ApAdjPower_CanBands_longformat_CxRxFB.mat','file') == 2
        load RCT1_ApAdjPower_CanBands_longformat_CxRxFB.mat
        NewRowInds = [(ss-1)*height(dataTable_currID)+1 ss*height(dataTable_currID)];
        RCT1_ApAdjPower_longformat(NewRowInds(1,1):NewRowInds(1,2),:) = dataTable_currID;
        clear NewRowInds
    else % first ID
        RCT1_ApAdjPower_longformat = dataTable_currID;
    end

    % save the table
    save('RCT1_ApAdjPower_CanBands_longformat_CxRxFB.mat','RCT1_ApAdjPower_longformat')
    
    % clear up
    clear ClinicalID ClinicalID1 IDfull IDparts IDcode Site Sessions N_trials 
    clear Nr_Rsq S_nRsq T_nRsq A_nRsq F_nRsq Nr_Aapow S_aapow T_aapow A_aapow F_aapow
    clear CurrID dataTable_currID S_Neps T_Neps A_Neps F_Neps 
    clear RCT1_ApAdjPower_longformat
    clear DATA

end

% clear up all
clear Freqs Frequencies_N Freqs_inds Condition Conditions_N Frequency_bands Frequencies Regions Regions_N
clear Freq_boundaries ss 
clear FOOOF_Narrow miniMADE_MOI_table MOI_location  Path_FOOOFspectra PSCcode_conversion


% save as csv file
load("RCT1_ApAdjPower_CanBands_longformat_CxRxFB.mat")
writetable(RCT1_ApAdjPower_longformat, 'RCT1_ApAdjPower_CanBands_longformat_CxRxFB.csv')
clear RCT1_ApAdjPower_longformat


