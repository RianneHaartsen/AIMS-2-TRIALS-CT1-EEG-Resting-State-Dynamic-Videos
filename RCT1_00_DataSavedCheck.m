%% Rbaclofen RCT1 00: initial data check for saved data on laptop vs data tracker from Teresa Del Bianco & Eleonora Broggi

% 1) Script reads in excel file with data and all files in the DATA folder,
% then checks which files are missing from the folder.
% 2) Converting the fieldtrip data to EEGlab for the pre-processing
% pipelines. 

% Created by dr. Rianne Haartsen, 2024

%% 1) Load in data and file names %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd xxx/DATA
TrackerData = readtable('arba-tracker-30May2024.xlsx','Sheet','Available-IDS-only_OWEY+Castor');

FilesSaved = dir("/xxx/fieldtrip_*.mat");
Conversion = readtable('/xxx/rct1_psc1_psc2.xlsx');

TrackerRH = TrackerData(:,[1 2 3]);
TrackerRH.Properties.VariableNames = {'id_eeg','tp','site'};

%% Check if fieldtrip file is present - psc1 or psc2
Nsaved = 0;
Nconverted = 0;
for ii = 1:height(TrackerRH)
    % for tracker ID
        if strcmp(TrackerRH.tp{ii},'1')
            CurNameA = strcat('fieldtrip_', TrackerRH.id_eeg{ii});
            CurSessionA = ('_test.mat');
        elseif strcmp(TrackerRH.tp{ii},'2') || strcmp(TrackerRH.tp{ii},'7')
            CurNameA = strcat('fieldtrip_', TrackerRH.id_eeg{ii});
            CurSessionA = ('_retest.mat');
        else
            warning('Unknown time point')
        end
    % check if ID is in tracker
    Saved = 0;
    for cc = 1:height(FilesSaved)
        if contains(FilesSaved(cc).name, CurNameA) && contains(FilesSaved(cc).name, CurSessionA)
            Saved = 1;
            FilesSaved(cc).InTracker = 'in tracker';
            NameEEGfile = FilesSaved(cc).name;
        end
    end
    clear cc

    if Saved == 0
        % if id was not converted from psc1 to psc2
        % find psc1 for psc2
            for tt = 1:height(Conversion)
                if strcmp(TrackerRH.id_eeg{ii}, num2str(Conversion.PSC2(tt)))
                    if strcmp(TrackerRH.tp{ii},'1')
                        CurNameB = strcat('fieldtrip_', Conversion.PSC1{tt});
                        CurSessionB = ('_test.mat');
                    elseif strcmp(TrackerRH.tp{ii},'2') || strcmp(TrackerRH.tp{ii},'7')
                        CurNameB = strcat('fieldtrip_', Conversion.PSC1{tt});
                        CurSessionB = ('_retest.mat');
                    else
                        warning('Unknown time point')
                    end
                end
            end
        % check if data in tracker
            for cc = 1:height(FilesSaved)
                if contains(FilesSaved(cc).name, CurNameB) && contains(FilesSaved(cc).name, CurSessionB)
                    Saved = 1;
                    Nconverted = Nconverted + 1;
                    FilesSaved(cc).InTracker = 'in tracker';
                    NameEEGfile = FilesSaved(cc).name;
                end
            end


    end

    % update tracker if saved or not
    if Saved == 1
        TrackerRH.Raw_ft_file{ii} = NameEEGfile;
        Nsaved = Nsaved+1;
    else
        TrackerRH.Raw_ft_file{ii} = [];
    end
    clear Saved cc
end

disp(Nsaved)
disp(Nconverted)

%% 
save("TrackerRH4.mat","TrackerRH")

save("Ft_folder_FilesSaved4.mat","FilesSaved")



%% 2) Convert fieldtrip saved data to EEGlab format %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('xxx/fieldtrip-20210629');  % add fieldtrip path
ft_defaults
addpath(genpath('xxx/eeglab2024.0'))
FilesSaved = dir("xxx/fieldtrip_*.mat");

for ii = 1:height(FilesSaved)

        % check if eeglab file exists
        cd xxx/00b_V1_05_export_EEGlab
        newname1 = strcat(extractAfter(FilesSaved(ii).name, 'fieldtrip_'));
        newname = strcat(extractBefore(newname1, 'test.mat'),'test_rawEEGlab.set');
        if exist(newname, 'file') == 2
            continue
        end

        % load fieldtrip data if eeglab file does not exist yet
        load(strcat(FilesSaved(ii).folder, '/',FilesSaved(ii).name))
        % rename header variables so recognised by ft/eeglab
        if exist('ft_data', 'var')
            ft = ft_data;
        end
        ft.enobio_hdr.Fs = ft.enobio_hdr.fs;
        ft.enobio_hdr.nChans = ft.enobio_hdr.numChans;
        ft.enobio_hdr.nSamplesPre = 0;
        ft.enobio_hdr.nTrials = 1;
        ft.enobio_hdr.label = ft.label;
        % add offset and duration to events field
        for tt = 1:height(ft.events)
            ft.events(tt).offset = 0;
            ft.events(tt).duration = 1;
        end
        clear tt
        % convert ft to eeglab
        eeglab_data = fieldtrip2eeglab(ft.enobio_hdr,cat(3,ft.trial{:}), ft.events);
        % save the eeglab file
        cd xxx/00b_V1_05_export_EEGlab
        newname1 = strcat(extractAfter(FilesSaved(ii).name, 'fieldtrip_'));
        newname = strcat(extractBefore(newname1, 'test.mat'),'test_rawEEGlab.set');
        pop_saveset(eeglab_data, 'filename', newname);

        % clean up
        clear ft newname eeglab_data
end

save("Ft_folder_FilesSaved4.mat","FilesSaved")


%% 3) Check for missing events %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd xxx/DATA
load TrackerRH4.mat

for ii = 1:height(TrackerRH)
        % load fieldtrip data
        if ~isempty(TrackerRH.Raw_ft_file{ii})
            load(strcat('xxx/',TrackerRH.Raw_ft_file{ii}))
            % rename header variables so recognised by ft/eeglab
            if exist('ft_data', 'var')
                ft = ft_data;
            end
            % check N events 
            if ~isempty(ft.events)
                N_EEGevents = size(ft.events,1);
            else
                N_EEGevents = 0;
            end
        else
            N_EEGevents = NaN;
        end
        % add to tracker
        TrackerRH.N_EEGevents{ii} = N_EEGevents;
        
        % clean up
        clear ft N_EEGevents ft_data
end


%% Save the tracker file for further preprocessing

save("TrackerRH4.mat","TrackerRH")
