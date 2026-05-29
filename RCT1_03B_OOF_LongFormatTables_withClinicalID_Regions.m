%% RCT_03B: 1/f features for statistical analyses

% This script takes csv files from Python that fitted one-over-f models
% onto the absolute power data from script RCT1_02B. It takes the 1/f
% offset, 1/f slope, and Rsq and error fit of the 1/f model per condition
% per region, for the narrow range (3-28Hz) and the wide range (1-120Hz).
% Only features with extracted on power spectra based on 20 or more epochs
% and with a 1/f Rsq higher than .95 are included in the tables, otherwise
% the value is set to NaN. 

% Here we are the 28 dataset with bump artefacts in the frequencies above
% 500Hz for the 1/f features across the wide range. 

% Data are saved in a .csv file for further statistical analyses in
% Rstudio.

% Rianne Haartsen, PhD.; Feb - 2025
% Birkbeck University of London

%% Extract FOOOF metrics and save into long format table
% load data and info
MOI_location = 'xxx/DATA/01B_PreprocOOF/';
cd(MOI_location)
load miniMADE_MeasuresOfInterest_report.mat
Path_FOOOFspectra = 'xxx/DATA/01B_PreprocOOF/FOOOF_powerspectra_Region/';
FOOOF_Narrow = readtable(strcat(Path_FOOOFspectra, 'Narrow_FBoI/RCT1_Region_Narrow_data.csv'),'VariableNamingRule','preserve'); 
FOOOF_Wide = readtable(strcat(Path_FOOOFspectra, 'Wide_3_120Hz_LNint/RCT1_Region_Wide_3_120Hz_LNint_data.csv'),'VariableNamingRule','preserve'); 
PSCcode_conversion = readtable('xxx/rct1_psc1_psc2.xlsx');

% Settings
Conditions_N = 4;
Regions_N = 4;
Ranges_N = 2;

Condition = [repmat({'Social'}, Regions_N*Ranges_N, 1); ...
    repmat({'Toy'}, Regions_N*Ranges_N, 1);...
    repmat({'Abstract'}, Regions_N*Ranges_N, 1);...
    repmat({'Fix'}, Regions_N*Ranges_N, 1)];
Regions = repmat([repmat({'Frontal'}, Ranges_N, 1); ...
    repmat({'Central'}, Ranges_N, 1);...
    repmat({'Parietal'}, Ranges_N, 1);...
    repmat({'Occipital'}, Ranges_N, 1)],Conditions_N,1);
Ranges = repmat({'Narrow';'Wide'},Conditions_N*Regions_N,1);

% Spectral power: correct for HF bumps in spectrum
% based on visual inspection of the power spectra (1-120Hz)
% exclude for Low vs High Freq: set High Freq to NaN
Subjects_with_HFbumps = [xxx];


% Loop through participants
for ss = 1:height(miniMADE_MOI_table)

    % info
    CurrID = miniMADE_MOI_table.ID{ss};
    HFBUMPS_present = 0;
    if ismember(ss, Subjects_with_HFbumps)
        HFBUMPS_present = 1;
    end
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

    % Social condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S_Neps = miniMADE_MOI_table.Seps(ss);
    if S_Neps >= 20
        % find 1/f results for this ID & condition
            % Frontal
            % Narrow
                Name = strcat(CurrID, '_Soc_Frontal');
                Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
                if ~isempty(Nrow) && S_Neps >= 20
                    fNr_Rsq = FOOOF_Narrow{Nrow,2};
                    fNr_FitE = FOOOF_Narrow{Nrow,3};
                    if fNr_Rsq > .95
                        fNr_offset = FOOOF_Narrow.intercept(Nrow);
                        fNr_slope = FOOOF_Narrow.slope(Nrow);
                    else
                        fNr_offset = NaN;
                        fNr_slope = NaN;
                    end
                else
                    fNr_Rsq = NaN;
                    fNr_FitE = NaN;
                    fNr_offset = NaN;
                    fNr_slope = NaN;
                end
                clear Nrow Name
            % Wide
                Name = strcat(CurrID, '_Soc_Frontal');
                Nrow = find(strcmp(FOOOF_Wide.Var1,Name)==1);
                if ~isempty(Nrow) && S_Neps >= 20 && HFBUMPS_present ~= 1
                    fWr_Rsq = FOOOF_Wide{Nrow,2};
                    fWr_FitE = FOOOF_Wide{Nrow,3};
                    if fWr_Rsq > .95
                        fWr_offset = FOOOF_Wide.intercept(Nrow);
                        fWr_slope = FOOOF_Wide.slope(Nrow);
                    else
                        fWr_offset = NaN;
                        fWr_slope = NaN;
                    end
                else
                    fWr_Rsq = NaN;
                    fWr_FitE = NaN;
                    fWr_offset = NaN;
                    fWr_slope = NaN;
                end
                clear Nrow Name

            % Central
            % Narrow
                Name = strcat(CurrID, '_Soc_Central');
                Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
                if ~isempty(Nrow) && S_Neps >= 20
                    cNr_Rsq = FOOOF_Narrow{Nrow,2};
                    cNr_FitE = FOOOF_Narrow{Nrow,3};
                    if cNr_Rsq > .95
                        cNr_offset = FOOOF_Narrow.intercept(Nrow);
                        cNr_slope = FOOOF_Narrow.slope(Nrow);
                    else
                        cNr_offset = NaN;
                        cNr_slope = NaN;
                    end
                else
                    cNr_Rsq = NaN;
                    cNr_FitE = NaN;
                    cNr_offset = NaN;
                    cNr_slope = NaN;
                end
                clear Nrow Name
            % Wide
                Name = strcat(CurrID, '_Soc_Central');
                Nrow = find(strcmp(FOOOF_Wide.Var1,Name)==1);
                if ~isempty(Nrow) && S_Neps >= 20 && HFBUMPS_present ~= 1
                    cWr_Rsq = FOOOF_Wide{Nrow,2};
                    cWr_FitE = FOOOF_Wide{Nrow,3};
                    if cWr_Rsq > .95
                        cWr_offset = FOOOF_Wide.intercept(Nrow);
                        cWr_slope = FOOOF_Wide.slope(Nrow);
                    else
                        cWr_offset = NaN;
                        cWr_slope = NaN;
                    end
                else
                    cWr_Rsq = NaN;
                    cWr_FitE = NaN;
                    cWr_offset = NaN;
                    cWr_slope = NaN;
                end
                clear Nrow Name

            % Parietal
            % Narrow
                Name = strcat(CurrID, '_Soc_Parietal');
                Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
                if ~isempty(Nrow) && S_Neps >= 20
                    pNr_Rsq = FOOOF_Narrow{Nrow,2};
                    pNr_FitE = FOOOF_Narrow{Nrow,3};
                    if pNr_Rsq > .95
                        pNr_offset = FOOOF_Narrow.intercept(Nrow);
                        pNr_slope = FOOOF_Narrow.slope(Nrow);
                    else
                        pNr_offset = NaN;
                        pNr_slope = NaN;
                    end
                else
                    pNr_Rsq = NaN;
                    pNr_FitE = NaN;
                    pNr_offset = NaN;
                    pNr_slope = NaN;
                end
                clear Nrow Name
            % Wide
                Name = strcat(CurrID, '_Soc_Parietal');
                Nrow = find(strcmp(FOOOF_Wide.Var1,Name)==1);
                if ~isempty(Nrow) && S_Neps >= 20 && HFBUMPS_present ~= 1
                    pWr_Rsq = FOOOF_Wide{Nrow,2};
                    pWr_FitE = FOOOF_Wide{Nrow,3};
                    if pWr_Rsq > .95
                        pWr_offset = FOOOF_Wide.intercept(Nrow);
                        pWr_slope = FOOOF_Wide.slope(Nrow);
                    else
                        pWr_offset = NaN;
                        pWr_slope = NaN;
                    end
                else
                    pWr_Rsq = NaN;
                    pWr_FitE = NaN;
                    pWr_offset = NaN;
                    pWr_slope = NaN;
                end
                clear Nrow Name

            % Occipital
            % Narrow
                Name = strcat(CurrID, '_Soc_Occipital');
                Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
                if ~isempty(Nrow) && S_Neps >= 20
                    oNr_Rsq = FOOOF_Narrow{Nrow,2};
                    oNr_FitE = FOOOF_Narrow{Nrow,3};
                    if oNr_Rsq > .95
                        oNr_offset = FOOOF_Narrow.intercept(Nrow);
                        oNr_slope = FOOOF_Narrow.slope(Nrow);
                    else
                        oNr_offset = NaN;
                        oNr_slope = NaN;
                    end
                else
                    oNr_Rsq = NaN;
                    oNr_FitE = NaN;
                    oNr_offset = NaN;
                    oNr_slope = NaN;
                end
                clear Nrow Name
            % Wide
                Name = strcat(CurrID, '_Soc_Occipital');
                Nrow = find(strcmp(FOOOF_Wide.Var1,Name)==1);
                if ~isempty(Nrow) && S_Neps >= 20 && HFBUMPS_present ~= 1
                    oWr_Rsq = FOOOF_Wide{Nrow,2};
                    oWr_FitE = FOOOF_Wide{Nrow,3};
                    if oWr_Rsq > .95
                        oWr_offset = FOOOF_Wide.intercept(Nrow);
                        oWr_slope = FOOOF_Wide.slope(Nrow);
                    else
                        oWr_offset = NaN;
                        oWr_slope = NaN;
                    end
                else
                    oWr_Rsq = NaN;
                    oWr_FitE = NaN;
                    oWr_offset = NaN;
                    oWr_slope = NaN;
                end
                clear Nrow Name

    else % set all to nan
            % frontal
            fNr_Rsq = NaN;
            fNr_FitE = NaN;
            fNr_offset = NaN;
            fNr_slope = NaN;
            fWr_Rsq = NaN;
            fWr_FitE = NaN;
            fWr_offset = NaN;
            fWr_slope = NaN;
            % central
            cNr_Rsq = NaN;
            cNr_FitE = NaN;
            cNr_offset = NaN;
            cNr_slope = NaN;
            cWr_Rsq = NaN;
            cWr_FitE = NaN;
            cWr_offset = NaN;
            cWr_slope = NaN;
            % parietal
            pNr_Rsq = NaN;
            pNr_FitE = NaN;
            pNr_offset = NaN;
            pNr_slope = NaN;
            pWr_Rsq = NaN;
            pWr_FitE = NaN;
            pWr_offset = NaN;
            pWr_slope = NaN;
            % occipital
            oNr_Rsq = NaN;
            oNr_FitE = NaN;
            oNr_offset = NaN;
            oNr_slope = NaN;
            oWr_Rsq = NaN;
            oWr_FitE = NaN;
            oWr_offset = NaN;
            oWr_slope = NaN;
    end

    % summary
    S_Offset = [fNr_offset; fWr_offset; cNr_offset; cWr_offset; pNr_offset; pWr_offset; oNr_offset; oWr_offset];
    S_Slope = [fNr_slope; fWr_slope; cNr_slope; cWr_slope; pNr_slope; pWr_slope; oNr_slope; oWr_slope];
    S_Rsq = [fNr_Rsq; fWr_Rsq; cNr_Rsq; cWr_Rsq; pNr_Rsq; pWr_Rsq; oNr_Rsq; oWr_Rsq];
    S_FitE = [fNr_FitE; fWr_FitE; cNr_FitE; cWr_FitE; pNr_FitE; pWr_FitE; oNr_FitE; oWr_FitE];
    % clean up
    clear fNr_offset fWr_offset cNr_offset cWr_offset pNr_offset pWr_offset oNr_offset oWr_offset
    clear fNr_slope fWr_slope cNr_slope cWr_slope pNr_slope pWr_slope oNr_slope oWr_slope
    clear fNr_Rsq fWr_Rsq cNr_Rsq cWr_Rsq pNr_Rsq pWr_Rsq oNr_Rsq oWr_Rsq
    clear fNr_FitE fWr_FitE cNr_FitE cWr_FitE pNr_FitE pWr_FitE oNr_FitE oWr_FitE


    % Toy condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T_Neps = miniMADE_MOI_table.Teps(ss);
    if T_Neps >= 20
        % find 1/f results for this ID & condition
            % Frontal
            % Narrow
                Name = strcat(CurrID, '_Toy_Frontal');
                Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
                if ~isempty(Nrow) && T_Neps >= 20
                    fNr_Rsq = FOOOF_Narrow{Nrow,2};
                    fNr_FitE = FOOOF_Narrow{Nrow,3};
                    if fNr_Rsq > .95
                        fNr_offset = FOOOF_Narrow.intercept(Nrow);
                        fNr_slope = FOOOF_Narrow.slope(Nrow);
                    else
                        fNr_offset = NaN;
                        fNr_slope = NaN;
                    end
                else
                    fNr_Rsq = NaN;
                    fNr_FitE = NaN;
                    fNr_offset = NaN;
                    fNr_slope = NaN;
                end
                clear Nrow Name
            % Wide
                Name = strcat(CurrID, '_Toy_Frontal');
                Nrow = find(strcmp(FOOOF_Wide.Var1,Name)==1);
                if ~isempty(Nrow) && T_Neps >= 20 && HFBUMPS_present ~= 1
                    fWr_Rsq = FOOOF_Wide{Nrow,2};
                    fWr_FitE = FOOOF_Wide{Nrow,3};
                    if fWr_Rsq > .95
                        fWr_offset = FOOOF_Wide.intercept(Nrow);
                        fWr_slope = FOOOF_Wide.slope(Nrow);
                    else
                        fWr_offset = NaN;
                        fWr_slope = NaN;
                    end
                else
                    fWr_Rsq = NaN;
                    fWr_FitE = NaN;
                    fWr_offset = NaN;
                    fWr_slope = NaN;
                end
                clear Nrow Name

            % Central
            % Narrow
                Name = strcat(CurrID, '_Toy_Central');
                Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
                if ~isempty(Nrow) && T_Neps >= 20
                    cNr_Rsq = FOOOF_Narrow{Nrow,2};
                    cNr_FitE = FOOOF_Narrow{Nrow,3};
                    if cNr_Rsq > .95
                        cNr_offset = FOOOF_Narrow.intercept(Nrow);
                        cNr_slope = FOOOF_Narrow.slope(Nrow);
                    else
                        cNr_offset = NaN;
                        cNr_slope = NaN;
                    end
                else
                    cNr_Rsq = NaN;
                    cNr_FitE = NaN;
                    cNr_offset = NaN;
                    cNr_slope = NaN;
                end
                clear Nrow Name
            % Wide
                Name = strcat(CurrID, '_Toy_Central');
                Nrow = find(strcmp(FOOOF_Wide.Var1,Name)==1);
                if ~isempty(Nrow) && T_Neps >= 20 && HFBUMPS_present ~= 1
                    cWr_Rsq = FOOOF_Wide{Nrow,2};
                    cWr_FitE = FOOOF_Wide{Nrow,3};
                    if cWr_Rsq > .95
                        cWr_offset = FOOOF_Wide.intercept(Nrow);
                        cWr_slope = FOOOF_Wide.slope(Nrow);
                    else
                        cWr_offset = NaN;
                        cWr_slope = NaN;
                    end
                else
                    cWr_Rsq = NaN;
                    cWr_FitE = NaN;
                    cWr_offset = NaN;
                    cWr_slope = NaN;
                end
                clear Nrow Name

            % Parietal
            % Narrow
                Name = strcat(CurrID, '_Toy_Parietal');
                Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
                if ~isempty(Nrow) && T_Neps >= 20
                    pNr_Rsq = FOOOF_Narrow{Nrow,2};
                    pNr_FitE = FOOOF_Narrow{Nrow,3};
                    if pNr_Rsq > .95
                        pNr_offset = FOOOF_Narrow.intercept(Nrow);
                        pNr_slope = FOOOF_Narrow.slope(Nrow);
                    else
                        pNr_offset = NaN;
                        pNr_slope = NaN;
                    end
                else
                    pNr_Rsq = NaN;
                    pNr_FitE = NaN;
                    pNr_offset = NaN;
                    pNr_slope = NaN;
                end
                clear Nrow Name
            % Wide
                Name = strcat(CurrID, '_Toy_Parietal');
                Nrow = find(strcmp(FOOOF_Wide.Var1,Name)==1);
                if ~isempty(Nrow) && T_Neps >= 20 && HFBUMPS_present ~= 1
                    pWr_Rsq = FOOOF_Wide{Nrow,2};
                    pWr_FitE = FOOOF_Wide{Nrow,3};
                    if pWr_Rsq > .95
                        pWr_offset = FOOOF_Wide.intercept(Nrow);
                        pWr_slope = FOOOF_Wide.slope(Nrow);
                    else
                        pWr_offset = NaN;
                        pWr_slope = NaN;
                    end
                else
                    pWr_Rsq = NaN;
                    pWr_FitE = NaN;
                    pWr_offset = NaN;
                    pWr_slope = NaN;
                end
                clear Nrow Name

            % Occipital
            % Narrow
                Name = strcat(CurrID, '_Toy_Occipital');
                Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
                if ~isempty(Nrow) && T_Neps >= 20
                    oNr_Rsq = FOOOF_Narrow{Nrow,2};
                    oNr_FitE = FOOOF_Narrow{Nrow,3};
                    if oNr_Rsq > .95
                        oNr_offset = FOOOF_Narrow.intercept(Nrow);
                        oNr_slope = FOOOF_Narrow.slope(Nrow);
                    else
                        oNr_offset = NaN;
                        oNr_slope = NaN;
                    end
                else
                    oNr_Rsq = NaN;
                    oNr_FitE = NaN;
                    oNr_offset = NaN;
                    oNr_slope = NaN;
                end
                clear Nrow Name
            % Wide
                Name = strcat(CurrID, '_Toy_Occipital');
                Nrow = find(strcmp(FOOOF_Wide.Var1,Name)==1);
                if ~isempty(Nrow) && T_Neps >= 20 && HFBUMPS_present ~= 1
                    oWr_Rsq = FOOOF_Wide{Nrow,2};
                    oWr_FitE = FOOOF_Wide{Nrow,3};
                    if oWr_Rsq > .95
                        oWr_offset = FOOOF_Wide.intercept(Nrow);
                        oWr_slope = FOOOF_Wide.slope(Nrow);
                    else
                        oWr_offset = NaN;
                        oWr_slope = NaN;
                    end
                else
                    oWr_Rsq = NaN;
                    oWr_FitE = NaN;
                    oWr_offset = NaN;
                    oWr_slope = NaN;
                end
                clear Nrow Name

    else % set all to nan
            % frontal
            fNr_Rsq = NaN;
            fNr_FitE = NaN;
            fNr_offset = NaN;
            fNr_slope = NaN;
            fWr_Rsq = NaN;
            fWr_FitE = NaN;
            fWr_offset = NaN;
            fWr_slope = NaN;
            % central
            cNr_Rsq = NaN;
            cNr_FitE = NaN;
            cNr_offset = NaN;
            cNr_slope = NaN;
            cWr_Rsq = NaN;
            cWr_FitE = NaN;
            cWr_offset = NaN;
            cWr_slope = NaN;
            % parietal
            pNr_Rsq = NaN;
            pNr_FitE = NaN;
            pNr_offset = NaN;
            pNr_slope = NaN;
            pWr_Rsq = NaN;
            pWr_FitE = NaN;
            pWr_offset = NaN;
            pWr_slope = NaN;
            % occipital
            oNr_Rsq = NaN;
            oNr_FitE = NaN;
            oNr_offset = NaN;
            oNr_slope = NaN;
            oWr_Rsq = NaN;
            oWr_FitE = NaN;
            oWr_offset = NaN;
            oWr_slope = NaN;
    end

    % summary
    T_Offset = [fNr_offset; fWr_offset; cNr_offset; cWr_offset; pNr_offset; pWr_offset; oNr_offset; oWr_offset];
    T_Slope = [fNr_slope; fWr_slope; cNr_slope; cWr_slope; pNr_slope; pWr_slope; oNr_slope; oWr_slope];
    T_Rsq = [fNr_Rsq; fWr_Rsq; cNr_Rsq; cWr_Rsq; pNr_Rsq; pWr_Rsq; oNr_Rsq; oWr_Rsq];
    T_FitE = [fNr_FitE; fWr_FitE; cNr_FitE; cWr_FitE; pNr_FitE; pWr_FitE; oNr_FitE; oWr_FitE];
    % clean up
    clear fNr_offset fWr_offset cNr_offset cWr_offset pNr_offset pWr_offset oNr_offset oWr_offset
    clear fNr_slope fWr_slope cNr_slope cWr_slope pNr_slope pWr_slope oNr_slope oWr_slope
    clear fNr_Rsq fWr_Rsq cNr_Rsq cWr_Rsq pNr_Rsq pWr_Rsq oNr_Rsq oWr_Rsq
    clear fNr_FitE fWr_FitE cNr_FitE cWr_FitE pNr_FitE pWr_FitE oNr_FitE oWr_FitE


    % Abstract condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A_Neps = miniMADE_MOI_table.Aeps(ss);
    if A_Neps >= 20
        % find 1/f results for this ID & condition
            % Frontal
            % Narrow
                Name = strcat(CurrID, '_Abs_Frontal');
                Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
                if ~isempty(Nrow) && A_Neps >= 20
                    fNr_Rsq = FOOOF_Narrow{Nrow,2};
                    fNr_FitE = FOOOF_Narrow{Nrow,3};
                    if fNr_Rsq > .95
                        fNr_offset = FOOOF_Narrow.intercept(Nrow);
                        fNr_slope = FOOOF_Narrow.slope(Nrow);
                    else
                        fNr_offset = NaN;
                        fNr_slope = NaN;
                    end
                else
                    fNr_Rsq = NaN;
                    fNr_FitE = NaN;
                    fNr_offset = NaN;
                    fNr_slope = NaN;
                end
                clear Nrow Name
            % Wide
                Name = strcat(CurrID, '_Abs_Frontal');
                Nrow = find(strcmp(FOOOF_Wide.Var1,Name)==1);
                if ~isempty(Nrow) && A_Neps >= 20 && HFBUMPS_present ~= 1
                    fWr_Rsq = FOOOF_Wide{Nrow,2};
                    fWr_FitE = FOOOF_Wide{Nrow,3};
                    if fWr_Rsq > .95
                        fWr_offset = FOOOF_Wide.intercept(Nrow);
                        fWr_slope = FOOOF_Wide.slope(Nrow);
                    else
                        fWr_offset = NaN;
                        fWr_slope = NaN;
                    end
                else
                    fWr_Rsq = NaN;
                    fWr_FitE = NaN;
                    fWr_offset = NaN;
                    fWr_slope = NaN;
                end
                clear Nrow Name

            % Central
            % Narrow
                Name = strcat(CurrID, '_Abs_Central');
                Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
                if ~isempty(Nrow) && A_Neps >= 20
                    cNr_Rsq = FOOOF_Narrow{Nrow,2};
                    cNr_FitE = FOOOF_Narrow{Nrow,3};
                    if cNr_Rsq > .95
                        cNr_offset = FOOOF_Narrow.intercept(Nrow);
                        cNr_slope = FOOOF_Narrow.slope(Nrow);
                    else
                        cNr_offset = NaN;
                        cNr_slope = NaN;
                    end
                else
                    cNr_Rsq = NaN;
                    cNr_FitE = NaN;
                    cNr_offset = NaN;
                    cNr_slope = NaN;
                end
                clear Nrow Name
            % Wide
                Name = strcat(CurrID, '_Abs_Central');
                Nrow = find(strcmp(FOOOF_Wide.Var1,Name)==1);
                if ~isempty(Nrow) && A_Neps >= 20 && HFBUMPS_present ~= 1
                    cWr_Rsq = FOOOF_Wide{Nrow,2};
                    cWr_FitE = FOOOF_Wide{Nrow,3};
                    if cWr_Rsq > .95
                        cWr_offset = FOOOF_Wide.intercept(Nrow);
                        cWr_slope = FOOOF_Wide.slope(Nrow);
                    else
                        cWr_offset = NaN;
                        cWr_slope = NaN;
                    end
                else
                    cWr_Rsq = NaN;
                    cWr_FitE = NaN;
                    cWr_offset = NaN;
                    cWr_slope = NaN;
                end
                clear Nrow Name

            % Parietal
            % Narrow
                Name = strcat(CurrID, '_Abs_Parietal');
                Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
                if ~isempty(Nrow) && A_Neps >= 20
                    pNr_Rsq = FOOOF_Narrow{Nrow,2};
                    pNr_FitE = FOOOF_Narrow{Nrow,3};
                    if pNr_Rsq > .95
                        pNr_offset = FOOOF_Narrow.intercept(Nrow);
                        pNr_slope = FOOOF_Narrow.slope(Nrow);
                    else
                        pNr_offset = NaN;
                        pNr_slope = NaN;
                    end
                else
                    pNr_Rsq = NaN;
                    pNr_FitE = NaN;
                    pNr_offset = NaN;
                    pNr_slope = NaN;
                end
                clear Nrow Name
            % Wide
                Name = strcat(CurrID, '_Abs_Parietal');
                Nrow = find(strcmp(FOOOF_Wide.Var1,Name)==1);
                if ~isempty(Nrow) && A_Neps >= 20 && HFBUMPS_present ~= 1
                    pWr_Rsq = FOOOF_Wide{Nrow,2};
                    pWr_FitE = FOOOF_Wide{Nrow,3};
                    if pWr_Rsq > .95
                        pWr_offset = FOOOF_Wide.intercept(Nrow);
                        pWr_slope = FOOOF_Wide.slope(Nrow);
                    else
                        pWr_offset = NaN;
                        pWr_slope = NaN;
                    end
                else
                    pWr_Rsq = NaN;
                    pWr_FitE = NaN;
                    pWr_offset = NaN;
                    pWr_slope = NaN;
                end
                clear Nrow Name

            % Occipital
            % Narrow
                Name = strcat(CurrID, '_Abs_Occipital');
                Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
                if ~isempty(Nrow) && A_Neps >= 20
                    oNr_Rsq = FOOOF_Narrow{Nrow,2};
                    oNr_FitE = FOOOF_Narrow{Nrow,3};
                    if oNr_Rsq > .95
                        oNr_offset = FOOOF_Narrow.intercept(Nrow);
                        oNr_slope = FOOOF_Narrow.slope(Nrow);
                    else
                        oNr_offset = NaN;
                        oNr_slope = NaN;
                    end
                else
                    oNr_Rsq = NaN;
                    oNr_FitE = NaN;
                    oNr_offset = NaN;
                    oNr_slope = NaN;
                end
                clear Nrow Name
            % Wide
                Name = strcat(CurrID, '_Abs_Occipital');
                Nrow = find(strcmp(FOOOF_Wide.Var1,Name)==1);
                if ~isempty(Nrow) && A_Neps >= 20 && HFBUMPS_present ~= 1
                    oWr_Rsq = FOOOF_Wide{Nrow,2};
                    oWr_FitE = FOOOF_Wide{Nrow,3};
                    if oWr_Rsq > .95
                        oWr_offset = FOOOF_Wide.intercept(Nrow);
                        oWr_slope = FOOOF_Wide.slope(Nrow);
                    else
                        oWr_offset = NaN;
                        oWr_slope = NaN;
                    end
                else
                    oWr_Rsq = NaN;
                    oWr_FitE = NaN;
                    oWr_offset = NaN;
                    oWr_slope = NaN;
                end
                clear Nrow Name

    else % set all to nan
            % frontal
            fNr_Rsq = NaN;
            fNr_FitE = NaN;
            fNr_offset = NaN;
            fNr_slope = NaN;
            fWr_Rsq = NaN;
            fWr_FitE = NaN;
            fWr_offset = NaN;
            fWr_slope = NaN;
            % central
            cNr_Rsq = NaN;
            cNr_FitE = NaN;
            cNr_offset = NaN;
            cNr_slope = NaN;
            cWr_Rsq = NaN;
            cWr_FitE = NaN;
            cWr_offset = NaN;
            cWr_slope = NaN;
            % parietal
            pNr_Rsq = NaN;
            pNr_FitE = NaN;
            pNr_offset = NaN;
            pNr_slope = NaN;
            pWr_Rsq = NaN;
            pWr_FitE = NaN;
            pWr_offset = NaN;
            pWr_slope = NaN;
            % occipital
            oNr_Rsq = NaN;
            oNr_FitE = NaN;
            oNr_offset = NaN;
            oNr_slope = NaN;
            oWr_Rsq = NaN;
            oWr_FitE = NaN;
            oWr_offset = NaN;
            oWr_slope = NaN;
    end

    % summary
    A_Offset = [fNr_offset; fWr_offset; cNr_offset; cWr_offset; pNr_offset; pWr_offset; oNr_offset; oWr_offset];
    A_Slope = [fNr_slope; fWr_slope; cNr_slope; cWr_slope; pNr_slope; pWr_slope; oNr_slope; oWr_slope];
    A_Rsq = [fNr_Rsq; fWr_Rsq; cNr_Rsq; cWr_Rsq; pNr_Rsq; pWr_Rsq; oNr_Rsq; oWr_Rsq];
    A_FitE = [fNr_FitE; fWr_FitE; cNr_FitE; cWr_FitE; pNr_FitE; pWr_FitE; oNr_FitE; oWr_FitE];
    % clean up
    clear  fNr_offset fWr_offset cNr_offset cWr_offset pNr_offset pWr_offset oNr_offset oWr_offset
    clear fNr_slope fWr_slope cNr_slope cWr_slope pNr_slope pWr_slope oNr_slope oWr_slope
    clear fNr_Rsq fWr_Rsq cNr_Rsq cWr_Rsq pNr_Rsq pWr_Rsq oNr_Rsq oWr_Rsq
    clear fNr_FitE fWr_FitE cNr_FitE cWr_FitE pNr_FitE pWr_FitE oNr_FitE oWr_FitE


    % Fixation condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F_Neps = miniMADE_MOI_table.Feps(ss);
    if F_Neps >= 20
        % find 1/f results for this ID & condition
            % Frontal
            % Narrow
                Name = strcat(CurrID, '_Fix_Frontal');
                Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
                if ~isempty(Nrow) && F_Neps >= 20
                    fNr_Rsq = FOOOF_Narrow{Nrow,2};
                    fNr_FitE = FOOOF_Narrow{Nrow,3};
                    if fNr_Rsq > .95
                        fNr_offset = FOOOF_Narrow.intercept(Nrow);
                        fNr_slope = FOOOF_Narrow.slope(Nrow);
                    else
                        fNr_offset = NaN;
                        fNr_slope = NaN;
                    end
                else
                    fNr_Rsq = NaN;
                    fNr_FitE = NaN;
                    fNr_offset = NaN;
                    fNr_slope = NaN;
                end
                clear Nrow Name
            % Wide
                Name = strcat(CurrID, '_Fix_Frontal');
                Nrow = find(strcmp(FOOOF_Wide.Var1,Name)==1);
                if ~isempty(Nrow) && F_Neps >= 20 && HFBUMPS_present ~= 1
                    fWr_Rsq = FOOOF_Wide{Nrow,2};
                    fWr_FitE = FOOOF_Wide{Nrow,3};
                    if fWr_Rsq > .95
                        fWr_offset = FOOOF_Wide.intercept(Nrow);
                        fWr_slope = FOOOF_Wide.slope(Nrow);
                    else
                        fWr_offset = NaN;
                        fWr_slope = NaN;
                    end
                else
                    fWr_Rsq = NaN;
                    fWr_FitE = NaN;
                    fWr_offset = NaN;
                    fWr_slope = NaN;
                end
                clear Nrow Name

            % Central
            % Narrow
                Name = strcat(CurrID, '_Fix_Central');
                Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
                if ~isempty(Nrow) && F_Neps >= 20
                    cNr_Rsq = FOOOF_Narrow{Nrow,2};
                    cNr_FitE = FOOOF_Narrow{Nrow,3};
                    if cNr_Rsq > .95
                        cNr_offset = FOOOF_Narrow.intercept(Nrow);
                        cNr_slope = FOOOF_Narrow.slope(Nrow);
                    else
                        cNr_offset = NaN;
                        cNr_slope = NaN;
                    end
                else
                    cNr_Rsq = NaN;
                    cNr_FitE = NaN;
                    cNr_offset = NaN;
                    cNr_slope = NaN;
                end
                clear Nrow Name
            % Wide
                Name = strcat(CurrID, '_Fix_Central');
                Nrow = find(strcmp(FOOOF_Wide.Var1,Name)==1);
                if ~isempty(Nrow) && F_Neps >= 20 && HFBUMPS_present ~= 1
                    cWr_Rsq = FOOOF_Wide{Nrow,2};
                    cWr_FitE = FOOOF_Wide{Nrow,3};
                    if cWr_Rsq > .95
                        cWr_offset = FOOOF_Wide.intercept(Nrow);
                        cWr_slope = FOOOF_Wide.slope(Nrow);
                    else
                        cWr_offset = NaN;
                        cWr_slope = NaN;
                    end
                else
                    cWr_Rsq = NaN;
                    cWr_FitE = NaN;
                    cWr_offset = NaN;
                    cWr_slope = NaN;
                end
                clear Nrow Name

            % Parietal
            % Narrow
                Name = strcat(CurrID, '_Fix_Parietal');
                Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
                if ~isempty(Nrow) && F_Neps >= 20
                    pNr_Rsq = FOOOF_Narrow{Nrow,2};
                    pNr_FitE = FOOOF_Narrow{Nrow,3};
                    if pNr_Rsq > .95
                        pNr_offset = FOOOF_Narrow.intercept(Nrow);
                        pNr_slope = FOOOF_Narrow.slope(Nrow);
                    else
                        pNr_offset = NaN;
                        pNr_slope = NaN;
                    end
                else
                    pNr_Rsq = NaN;
                    pNr_FitE = NaN;
                    pNr_offset = NaN;
                    pNr_slope = NaN;
                end
                clear Nrow Name
            % Wide
                Name = strcat(CurrID, '_Fix_Parietal');
                Nrow = find(strcmp(FOOOF_Wide.Var1,Name)==1);
                if ~isempty(Nrow) && F_Neps >= 20 && HFBUMPS_present ~= 1
                    pWr_Rsq = FOOOF_Wide{Nrow,2};
                    pWr_FitE = FOOOF_Wide{Nrow,3};
                    if pWr_Rsq > .95
                        pWr_offset = FOOOF_Wide.intercept(Nrow);
                        pWr_slope = FOOOF_Wide.slope(Nrow);
                    else
                        pWr_offset = NaN;
                        pWr_slope = NaN;
                    end
                else
                    pWr_Rsq = NaN;
                    pWr_FitE = NaN;
                    pWr_offset = NaN;
                    pWr_slope = NaN;
                end
                clear Nrow Name

            % Occipital
            % Narrow
                Name = strcat(CurrID, '_Fix_Occipital');
                Nrow = find(strcmp(FOOOF_Narrow.Var1,Name)==1);
                if ~isempty(Nrow) && F_Neps >= 20
                    oNr_Rsq = FOOOF_Narrow{Nrow,2};
                    oNr_FitE = FOOOF_Narrow{Nrow,3};
                    if oNr_Rsq > .95
                        oNr_offset = FOOOF_Narrow.intercept(Nrow);
                        oNr_slope = FOOOF_Narrow.slope(Nrow);
                    else
                        oNr_offset = NaN;
                        oNr_slope = NaN;
                    end
                else
                    oNr_Rsq = NaN;
                    oNr_FitE = NaN;
                    oNr_offset = NaN;
                    oNr_slope = NaN;
                end
                clear Nrow Name
            % Wide
                Name = strcat(CurrID, '_Fix_Occipital');
                Nrow = find(strcmp(FOOOF_Wide.Var1,Name)==1);
                if ~isempty(Nrow) && F_Neps >= 20 && HFBUMPS_present ~= 1
                    oWr_Rsq = FOOOF_Wide{Nrow,2};
                    oWr_FitE = FOOOF_Wide{Nrow,3};
                    if oWr_Rsq > .95
                        oWr_offset = FOOOF_Wide.intercept(Nrow);
                        oWr_slope = FOOOF_Wide.slope(Nrow);
                    else
                        oWr_offset = NaN;
                        oWr_slope = NaN;
                    end
                else
                    oWr_Rsq = NaN;
                    oWr_FitE = NaN;
                    oWr_offset = NaN;
                    oWr_slope = NaN;
                end
                clear Nrow Name

    else % set all to nan
            % frontal
            fNr_Rsq = NaN;
            fNr_FitE = NaN;
            fNr_offset = NaN;
            fNr_slope = NaN;
            fWr_Rsq = NaN;
            fWr_FitE = NaN;
            fWr_offset = NaN;
            fWr_slope = NaN;
            % central
            cNr_Rsq = NaN;
            cNr_FitE = NaN;
            cNr_offset = NaN;
            cNr_slope = NaN;
            cWr_Rsq = NaN;
            cWr_FitE = NaN;
            cWr_offset = NaN;
            cWr_slope = NaN;
            % parietal
            pNr_Rsq = NaN;
            pNr_FitE = NaN;
            pNr_offset = NaN;
            pNr_slope = NaN;
            pWr_Rsq = NaN;
            pWr_FitE = NaN;
            pWr_offset = NaN;
            pWr_slope = NaN;
            % occipital
            oNr_Rsq = NaN;
            oNr_FitE = NaN;
            oNr_offset = NaN;
            oNr_slope = NaN;
            oWr_Rsq = NaN;
            oWr_FitE = NaN;
            oWr_offset = NaN;
            oWr_slope = NaN;
    end

    % summary
    F_Offset = [fNr_offset; fWr_offset; cNr_offset; cWr_offset; pNr_offset; pWr_offset; oNr_offset; oWr_offset];
    F_Slope = [fNr_slope; fWr_slope; cNr_slope; cWr_slope; pNr_slope; pWr_slope; oNr_slope; oWr_slope];
    F_Rsq = [fNr_Rsq; fWr_Rsq; cNr_Rsq; cWr_Rsq; pNr_Rsq; pWr_Rsq; oNr_Rsq; oWr_Rsq];
    F_FitE = [fNr_FitE; fWr_FitE; cNr_FitE; cWr_FitE; pNr_FitE; pWr_FitE; oNr_FitE; oWr_FitE];
    % clean up
    clear  fNr_offset fWr_offset cNr_offset cWr_offset pNr_offset pWr_offset oNr_offset oWr_offset
    clear fNr_slope fWr_slope cNr_slope cWr_slope pNr_slope pWr_slope oNr_slope oWr_slope
    clear fNr_Rsq fWr_Rsq cNr_Rsq cWr_Rsq pNr_Rsq pWr_Rsq oNr_Rsq oWr_Rsq
    clear fNr_FitE fWr_FitE cNr_FitE cWr_FitE pNr_FitE pWr_FitE oNr_FitE oWr_FitE



    % Create table with new data
    % Columns: ID, condition (4), N trials, frequency range (2), 1/f offset,
    % 1/f slope, Rsq, Fit error
    N_trials = [repmat(S_Neps, Regions_N*Ranges_N,1); repmat(T_Neps, Regions_N*Ranges_N,1); ...
        repmat(A_Neps, Regions_N*Ranges_N,1); repmat(F_Neps, Regions_N*Ranges_N,1)];
    % Values of measurements
    Offset = [S_Offset; T_Offset; A_Offset; F_Offset];
    Slope = [S_Slope; T_Slope; A_Slope; F_Slope];
    Rsq = [S_Rsq; T_Rsq; A_Rsq; F_Rsq];
    FitError = [S_FitE; T_FitE; A_FitE; F_FitE];
    % assemble table
    dataTable_currID = table(ClinicalID, IDfull, IDcode, Site, Sessions, Condition, Regions, Ranges, N_trials, Offset, Slope, Rsq, FitError);
    % add current ID data to table
    cd(MOI_location)
    if exist('RCT1_1oF_longformat_CxRxFB.mat','file') == 2
        load RCT1_1oF_longformat_CxRxFB.mat
        NewRowInds = [(ss-1)*height(dataTable_currID)+1 ss*height(dataTable_currID)];
        RCT1_1oF_longformat(NewRowInds(1,1):NewRowInds(1,2),:) = dataTable_currID;
        clear NewRowInds
    else % first ID
        RCT1_1oF_longformat = dataTable_currID;
    end

    % save the table
    save('RCT1_1oF_longformat_CxRxFB.mat','RCT1_1oF_longformat')
    
    % clear up
    clear ClinicalID ClinicalID1 IDfull IDcode Sessions N_trials Offset Slope Rsq FitError IDparts HFBUMPS_present
    clear S_Offset T_Offset A_Offset F_Offset S_Slope T_Slope A_Slope F_Slope
    clear S_Rsq T_Rsq A_Rsq F_Rsq S_FitE T_FitE A_FitE F_FitE
    clear CurrID dataTable_currID S_Neps T_Neps A_Neps F_Neps 
    clear RCT1_1oF_longformat

end

%% Save table into csv format for Rstudio
cd(MOI_location)
load RCT1_1oF_longformat_CxRxFB.mat
% FOOOF in .csv file for R
cd xxx/DATA/01B_PreprocOOF/
writetable(RCT1_1oF_longformat, 'RCT1_1oF_longformat_CxRxFB.csv')
clear RCT1_1oF_longformat

