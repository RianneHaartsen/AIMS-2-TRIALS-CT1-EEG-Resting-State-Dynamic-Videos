%% RCT_03B: 1/f features for statistical analyses

% This script takes the tables with Metrics of Interest from FC (script
% RCT1_02C) and assembles them into long format for statistical analysis in 
% R studio. 
% The script averages FC values from 1 region to the other EEG channels
% (e.g. frontal to other channels = Frontal-R) for each condition. 

% Data are saved in a .csv file for further statistical analyses in
% Rstudio.

% Rianne Haartsen, PhD.; Feb - 2025
% Birkbeck University of London

%% Functional connectivity - dbWPLI
% Columns: ID, condition (4), Region (4) frequency band (2), N trials, dbWPLI
% frequency bands: 
% theta band: 6-7Hz, alpha band: 8-12Hz
% 8 rows per ID

MOI_location = 'xxx/DATA/01C_PreprocFC';
cd(MOI_location)
load miniMADE_MeasuresOfInterestFC_report.mat

PSCcode_conversion = readtable('xxx/rct1_psc1_psc2.xlsx');


% Settings
Conditions_N = 4;
Regions_N = 4;
Freqs_N = 2;

Condition = [repmat({'Social'}, Regions_N*Freqs_N, 1); ...
    repmat({'Toy'}, Regions_N*Freqs_N, 1);...
    repmat({'Abstract'}, Regions_N*Freqs_N, 1);...
    repmat({'Fix'}, Regions_N*Freqs_N, 1)];
Regions = repmat([repmat({'Frontal-R'}, Freqs_N, 1); ...
    repmat({'Central-R'}, Freqs_N, 1);...
    repmat({'Parietal-R'}, Freqs_N, 1);...
    repmat({'Occipital-R'}, Freqs_N, 1)],Conditions_N,1);
Frequencies = repmat({'Theta';'Alpha'},Conditions_N*Regions_N,1);

Freqs = miniMADE_MOIfc_table.Sfreqs{1};
Freq_boundaries = [6 7; 8 12];
Freqs_inds = [find(Freqs == Freq_boundaries(1,1)), find(Freqs == Freq_boundaries(1,2)); ...
    find(Freqs == Freq_boundaries(2,1)), find(Freqs == Freq_boundaries(2,2))];


% Loop through participants
for ss = 1:height(miniMADE_MOIfc_table)

    % info
    CurrID = miniMADE_MOIfc_table.ID{ss};
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
    if ~isempty(miniMADE_MOIfc_table.FC_file{ss})
        load(miniMADE_MOIfc_table.FC_file{ss})
    end

    % Social condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S_Neps = miniMADE_MOIfc_table.Seps(ss);
    if S_Neps >= 90 && strcmp(DATA_fc(1).Cond{1,1},'Soc')

        % calculate FC for different frequencies
        FCCurr = DATA_fc(1).FCdata.dbWPLI;
        FC_mnFCBand = zeros(size(Freqs_inds,1),19,19);
        for ff = 1:size(Freqs_inds,1)
            FC_mnFCBand(ff,:,:) = shiftdim(mean(FCCurr(Freqs_inds(ff,1):Freqs_inds(ff,2),:,:),1));
        end
        clear FCCurr ff
        
        % get indices for channels
        Ch_labels = DATA_fc(1).FCdata.Chan_label;
        Chind_front = find(ismember(Ch_labels,{'Fz','F3','F4'}));
        Chind_frontR = find(ismember(Ch_labels,{'Cz','C3','C4','Pz','P3','P4','PO7','Oz','PO8'}));
        Chind_cent = find(ismember(Ch_labels,{'Cz','C3','C4'}));
        Chind_centR = find(ismember(Ch_labels,{'Fz','F3','F4','Pz','P3','P4','PO7','Oz','PO8'}));
        Chind_pari = find(ismember(Ch_labels,{'Pz','P3','P4'}));
        Chind_pariR = find(ismember(Ch_labels,{'Fz','F3','F4','Cz','C3','C4','PO7','Oz','PO8'}));
        Chind_occip = find(ismember(Ch_labels,{'PO7','Oz','PO8'}));
        Chind_occipR = find(ismember(Ch_labels,{'Fz','F3','F4','Cz','C3','C4','Pz','P3','P4'}));


        % frontal
        FrR_theta = mean(shiftdim(FC_mnFCBand(1,Chind_front, Chind_frontR)),'all');
        FrR_alpha = mean(shiftdim(FC_mnFCBand(2,Chind_front, Chind_frontR)),'all');
        % central
        CeR_theta = mean(shiftdim(FC_mnFCBand(1,Chind_cent, Chind_centR)),'all');
        CeR_alpha = mean(shiftdim(FC_mnFCBand(2,Chind_cent, Chind_centR)),'all');
        % parietal
        PaR_theta = mean(shiftdim(FC_mnFCBand(1,Chind_pari, Chind_pariR)),'all');
        PaR_alpha = mean(shiftdim(FC_mnFCBand(2,Chind_pari, Chind_pariR)),'all');
        % occipital
        OcR_theta = mean(shiftdim(FC_mnFCBand(1,Chind_occip, Chind_occipR)),'all');
        OcR_alpha = mean(shiftdim(FC_mnFCBand(2,Chind_occip, Chind_occipR)),'all');
        
        FCmn = [FrR_theta; FrR_alpha; CeR_theta; CeR_alpha; PaR_theta; PaR_alpha; OcR_theta; OcR_alpha];

        clear FrR_theta FrR_alpha CeR_theta CeR_alpha PaR_theta PaR_alpha OcR_theta OcR_alpha FC_mnFCBand
        clear Ch_labels Chind_front Chind_frontR Chind_cent Chind_centR Chind_pari Chind_pariR Chind_occip Chind_occipR
        

    else % no sufficient data for this condition
        FCmn = nan(Regions_N*Freqs_N,1);
    end
    % Rename vars
    S_fc = FCmn; clear FCmn


    % Toy condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T_Neps = miniMADE_MOIfc_table.Teps(ss);
    if T_Neps >= 90 && strcmp(DATA_fc(2).Cond{1,1},'Toy')

        % calculate FC for different frequencies
        FCCurr = DATA_fc(2).FCdata.dbWPLI;
        FC_mnFCBand = zeros(size(Freqs_inds,1),19,19);
        for ff = 1:size(Freqs_inds,1)
            FC_mnFCBand(ff,:,:) = shiftdim(mean(FCCurr(Freqs_inds(ff,1):Freqs_inds(ff,2),:,:),1));
        end
        clear FCCurr ff
        
        % get indices for channels
        Ch_labels = DATA_fc(2).FCdata.Chan_label;
        Chind_front = find(ismember(Ch_labels,{'Fz','F3','F4'}));
        Chind_frontR = find(ismember(Ch_labels,{'Cz','C3','C4','Pz','P3','P4','PO7','Oz','PO8'}));
        Chind_cent = find(ismember(Ch_labels,{'Cz','C3','C4'}));
        Chind_centR = find(ismember(Ch_labels,{'Fz','F3','F4','Pz','P3','P4','PO7','Oz','PO8'}));
        Chind_pari = find(ismember(Ch_labels,{'Pz','P3','P4'}));
        Chind_pariR = find(ismember(Ch_labels,{'Fz','F3','F4','Cz','C3','C4','PO7','Oz','PO8'}));
        Chind_occip = find(ismember(Ch_labels,{'PO7','Oz','PO8'}));
        Chind_occipR = find(ismember(Ch_labels,{'Fz','F3','F4','Cz','C3','C4','Pz','P3','P4'}));


        % frontal
        FrR_theta = mean(shiftdim(FC_mnFCBand(1,Chind_front, Chind_frontR)),'all');
        FrR_alpha = mean(shiftdim(FC_mnFCBand(2,Chind_front, Chind_frontR)),'all');
        % central
        CeR_theta = mean(shiftdim(FC_mnFCBand(1,Chind_cent, Chind_centR)),'all');
        CeR_alpha = mean(shiftdim(FC_mnFCBand(2,Chind_cent, Chind_centR)),'all');
        % parietal
        PaR_theta = mean(shiftdim(FC_mnFCBand(1,Chind_pari, Chind_pariR)),'all');
        PaR_alpha = mean(shiftdim(FC_mnFCBand(2,Chind_pari, Chind_pariR)),'all');
        % occipital
        OcR_theta = mean(shiftdim(FC_mnFCBand(1,Chind_occip, Chind_occipR)),'all');
        OcR_alpha = mean(shiftdim(FC_mnFCBand(2,Chind_occip, Chind_occipR)),'all');
        
        FCmn = [FrR_theta; FrR_alpha; CeR_theta; CeR_alpha; PaR_theta; PaR_alpha; OcR_theta; OcR_alpha];

        clear FrR_theta FrR_alpha CeR_theta CeR_alpha PaR_theta PaR_alpha OcR_theta OcR_alpha FC_mnFCBand
        clear Ch_labels Chind_front Chind_frontR Chind_cent Chind_centR Chind_pari Chind_pariR Chind_occip Chind_occipR
        

    else % no sufficient data for this condition
        FCmn = nan(Regions_N*Freqs_N,1);
    end
    % Rename vars
    T_fc = FCmn; clear FCmn


    % Abstract condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A_Neps = miniMADE_MOIfc_table.Aeps(ss);
    if A_Neps >= 90 && strcmp(DATA_fc(3).Cond{1,1},'Abs')

        % calculate FC for different frequencies
        FCCurr = DATA_fc(3).FCdata.dbWPLI;
        FC_mnFCBand = zeros(size(Freqs_inds,1),19,19);
        for ff = 1:size(Freqs_inds,1)
            FC_mnFCBand(ff,:,:) = shiftdim(mean(FCCurr(Freqs_inds(ff,1):Freqs_inds(ff,2),:,:),1));
        end
        clear FCCurr ff
        
        % get indices for channels
        Ch_labels = DATA_fc(3).FCdata.Chan_label;
        Chind_front = find(ismember(Ch_labels,{'Fz','F3','F4'}));
        Chind_frontR = find(ismember(Ch_labels,{'Cz','C3','C4','Pz','P3','P4','PO7','Oz','PO8'}));
        Chind_cent = find(ismember(Ch_labels,{'Cz','C3','C4'}));
        Chind_centR = find(ismember(Ch_labels,{'Fz','F3','F4','Pz','P3','P4','PO7','Oz','PO8'}));
        Chind_pari = find(ismember(Ch_labels,{'Pz','P3','P4'}));
        Chind_pariR = find(ismember(Ch_labels,{'Fz','F3','F4','Cz','C3','C4','PO7','Oz','PO8'}));
        Chind_occip = find(ismember(Ch_labels,{'PO7','Oz','PO8'}));
        Chind_occipR = find(ismember(Ch_labels,{'Fz','F3','F4','Cz','C3','C4','Pz','P3','P4'}));


        % frontal
        FrR_theta = mean(shiftdim(FC_mnFCBand(1,Chind_front, Chind_frontR)),'all');
        FrR_alpha = mean(shiftdim(FC_mnFCBand(2,Chind_front, Chind_frontR)),'all');
        % central
        CeR_theta = mean(shiftdim(FC_mnFCBand(1,Chind_cent, Chind_centR)),'all');
        CeR_alpha = mean(shiftdim(FC_mnFCBand(2,Chind_cent, Chind_centR)),'all');
        % parietal
        PaR_theta = mean(shiftdim(FC_mnFCBand(1,Chind_pari, Chind_pariR)),'all');
        PaR_alpha = mean(shiftdim(FC_mnFCBand(2,Chind_pari, Chind_pariR)),'all');
        % occipital
        OcR_theta = mean(shiftdim(FC_mnFCBand(1,Chind_occip, Chind_occipR)),'all');
        OcR_alpha = mean(shiftdim(FC_mnFCBand(2,Chind_occip, Chind_occipR)),'all');
        
        FCmn = [FrR_theta; FrR_alpha; CeR_theta; CeR_alpha; PaR_theta; PaR_alpha; OcR_theta; OcR_alpha];

        clear FrR_theta FrR_alpha CeR_theta CeR_alpha PaR_theta PaR_alpha OcR_theta OcR_alpha FC_mnFCBand
        clear Ch_labels Chind_front Chind_frontR Chind_cent Chind_centR Chind_pari Chind_pariR Chind_occip Chind_occipR
        

    else % no sufficient data for this condition
        FCmn = nan(Regions_N*Freqs_N,1);
    end
    % Rename vars
    A_fc = FCmn; clear FCmn

    
    % Fixation condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F_Neps = miniMADE_MOIfc_table.Feps(ss);
    if F_Neps >= 90 && strcmp(DATA_fc(4).Cond{1,1},'Fix')

        % calculate FC for different frequencies
        FCCurr = DATA_fc(4).FCdata.dbWPLI;
        FC_mnFCBand = zeros(size(Freqs_inds,1),19,19);
        for ff = 1:size(Freqs_inds,1)
            FC_mnFCBand(ff,:,:) = shiftdim(mean(FCCurr(Freqs_inds(ff,1):Freqs_inds(ff,2),:,:),1));
        end
        clear FCCurr ff
        
        % get indices for channels
        Ch_labels = DATA_fc(4).FCdata.Chan_label;
        Chind_front = find(ismember(Ch_labels,{'Fz','F3','F4'}));
        Chind_frontR = find(ismember(Ch_labels,{'Cz','C3','C4','Pz','P3','P4','PO7','Oz','PO8'}));
        Chind_cent = find(ismember(Ch_labels,{'Cz','C3','C4'}));
        Chind_centR = find(ismember(Ch_labels,{'Fz','F3','F4','Pz','P3','P4','PO7','Oz','PO8'}));
        Chind_pari = find(ismember(Ch_labels,{'Pz','P3','P4'}));
        Chind_pariR = find(ismember(Ch_labels,{'Fz','F3','F4','Cz','C3','C4','PO7','Oz','PO8'}));
        Chind_occip = find(ismember(Ch_labels,{'PO7','Oz','PO8'}));
        Chind_occipR = find(ismember(Ch_labels,{'Fz','F3','F4','Cz','C3','C4','Pz','P3','P4'}));


        % frontal
        FrR_theta = mean(shiftdim(FC_mnFCBand(1,Chind_front, Chind_frontR)),'all');
        FrR_alpha = mean(shiftdim(FC_mnFCBand(2,Chind_front, Chind_frontR)),'all');
        % central
        CeR_theta = mean(shiftdim(FC_mnFCBand(1,Chind_cent, Chind_centR)),'all');
        CeR_alpha = mean(shiftdim(FC_mnFCBand(2,Chind_cent, Chind_centR)),'all');
        % parietal
        PaR_theta = mean(shiftdim(FC_mnFCBand(1,Chind_pari, Chind_pariR)),'all');
        PaR_alpha = mean(shiftdim(FC_mnFCBand(2,Chind_pari, Chind_pariR)),'all');
        % occipital
        OcR_theta = mean(shiftdim(FC_mnFCBand(1,Chind_occip, Chind_occipR)),'all');
        OcR_alpha = mean(shiftdim(FC_mnFCBand(2,Chind_occip, Chind_occipR)),'all');
        
        FCmn = [FrR_theta; FrR_alpha; CeR_theta; CeR_alpha; PaR_theta; PaR_alpha; OcR_theta; OcR_alpha];

        clear FrR_theta FrR_alpha CeR_theta CeR_alpha PaR_theta PaR_alpha OcR_theta OcR_alpha FC_mnFCBand
        clear Ch_labels Chind_front Chind_frontR Chind_cent Chind_centR Chind_pari Chind_pariR Chind_occip Chind_occipR
        

    else % no sufficient data for this condition
        FCmn = nan(Regions_N*Freqs_N,1);
    end
    % Rename vars
    F_fc = FCmn; clear FCmn



    % Create table with new data
    % Columns: IDfull, ID, session, condition (4), region (4) frequency
    % band (2), N trials, FC (dbWPLI)
    % Ntrials
    N_trials = [repmat(S_Neps, Regions_N*Freqs_N,1); repmat(T_Neps, Regions_N*Freqs_N,1); ...
        repmat(A_Neps, Regions_N*Freqs_N,1); repmat(F_Neps, Regions_N*Freqs_N,1)];
    % Values of measurements: FC
    FC = [S_fc; T_fc; A_fc; F_fc];
    % assemble table
    dataTable_currID = table(ClinicalID, IDfull, IDcode, Site, Sessions, Condition, Regions, Frequencies, N_trials, FC);

    % add current ID data to table
    cd(MOI_location)
    if exist('RCT1_FC_longformat_CxRxFB.mat','file') == 2
        load RCT1_FC_longformat_CxRxFB.mat
        NewRowInds = [(ss-1)*height(dataTable_currID)+1 ss*height(dataTable_currID)];
        RCT1_FC_longformat(NewRowInds(1,1):NewRowInds(1,2),:) = dataTable_currID;
        clear NewRowInds
    else % first ID
        RCT1_FC_longformat = dataTable_currID;
    end

    % save the table
    save('RCT1_FC_longformat_CxRxFB.mat','RCT1_FC_longformat')
    
    % clear up
    clear ClinicalID IDparts IDfull IDcode Site Sessions N_trials FC
    clear S_fc T_fc A_fc F_fc FC DATA DATA_fc
    clear CurrID dataTable_currID S_Neps T_Neps A_Neps F_Neps 
    clear RCT1_FC_longformat

end


%% Save FC in theta and alpha range in .csv file for R
cd /Users/riannehaartsen/Documents/02g_ArbaclofenTrial/RCT1_EEG/DATA/01C_PreprocFC
load RCT1_FC_longformat_CxRxFB.mat
writetable(RCT1_FC_longformat, 'RCT1_FC_longformat_CxRxFB.csv')
clear RCT1_FC_longformat
