%% RCT_02Ca: Visual inspection of functional connectivity data

% This script takes preprocessing report from the pipeline and creates a
% figure of the data quantity after each round of artefact rejection, and
% another one on clean data epochs per dataset. 
% The other half of the script visualises the global dbWPLI spectra
% (1-32Hz) and global dbWPLI within the theta and alpha frequency band for
% the different conditions. 

% Dr. Rianne Haartsen, Sept 2024
% Birkbeck University of London

%% 
clear all
clc

output_location = 'xxx/DATA/01C_PreprocFC';
cd(output_location)
load miniMADE_preprocessing_report.mat

%% Ns drop after each AR round
Ndrops = figure('Position', [378, 400, 1041, 396]);
Ns_all = [miniMADE_report_table.total_epochs_before_artifact_rejection miniMADE_report_table.Neps_postAR1 miniMADE_report_table.total_epochs_after_artifact_rejection];
MaxN = max(Ns_all,[],'all');
YUpLim = MaxN + 6;
subplot 511
    bar(Ns_all)
    xlim([0 50.5])
    ylabel('N epochs'); xlabel('Datasets'); ylim([0 YUpLim])
    legend({'Present','After AR1 - blinks', 'After AR2 - other'})
subplot 512
    bar(Ns_all)
    xlim([50.5 100.5])
    ylabel('N epochs'); xlabel('Datasets'); ylim([0 YUpLim])
    legend({'Present','After AR1 - blinks', 'After AR2 - other'})
subplot 513
    bar(Ns_all)
    xlim([100.5 150.5])
    ylabel('N epochs'); xlabel('Datasets'); ylim([0 YUpLim])
    legend({'Present','After AR1 - blinks', 'After AR2 - other'})
subplot 514
    bar(Ns_all)
    xlim([150.5 200.5])
    ylabel('N epochs'); xlabel('Datasets'); ylim([0 YUpLim])
    legend({'Present','After AR1 - blinks', 'After AR2 - other'})
subplot 515
    bar(Ns_all)
    xlim([200.5 250.5])
    ylabel('N epochs'); xlabel('Datasets'); ylim([0 YUpLim])
    legend({'Present','After AR1 - blinks', 'After AR2 - other'})


sgtitle('Ns after AR1 & 2')

saveas(Ndrops, 'Ntrials_dropout.png');
close(Ndrops)


%% Data quantity
cd(output_location)
load miniMADE_MeasuresOfInterestFC_report.mat

DataQuantityFC = figure;

subplot(2,1,1) % N epochs per dataset
y = [miniMADE_MOIfc_table.Seps miniMADE_MOIfc_table.Teps miniMADE_MOIfc_table.Aeps miniMADE_MOIfc_table.Feps];
bar(y)
ylabel('N epochs'); xlabel('Dataset')
xticks(1:10:length(miniMADE_MOIfc_table.ID))
xticklabels(1:10:length(miniMADE_MOIfc_table.ID))
xtickangle(45)
ylim([0 850])
ax = gca; ax.XAxis.FontSize = 12; ax.YAxis.FontSize = 12;
yline(90,'LineStyle','-')
legend({'Soc', 'Toy','Abstract','Fixation','Threshold FC'},'Location','southoutside','Orientation','horizontal', 'FontSize', 12)
title('Data quantity per dataset','FontSize',12)

subplot(2,1,2) % inclusion rates overall
y2 = [sum(miniMADE_MOIfc_table.Seps >= 90) sum(miniMADE_MOIfc_table.Teps >= 90) ...
    sum(miniMADE_MOIfc_table.Aeps >= 90) sum(miniMADE_MOIfc_table.Feps >= 90)];
bar_handle = bar(y2);
bar_handle.FaceColor = 'flat'; % Use 'flat' to allow setting colors individually
bar_handle.CData = [0 0.4470 0.7410; 
                    0.8500 0.3250 0.0980; 
                    0.9290 0.6940 0.1250; 
                    0.4940 0.1840 0.5560]; 
for i = 1:length(y2)
    text(i, y2(i) + 1, sprintf('%.1f%%', y2(i)/length(miniMADE_MOIfc_table.ID)*100), ... % Add 1 to position above the bar
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 12);
end
ylabel('N datasets'); xlabel('Condition')
xticks(1:1:4)
xticklabels({'Soc', 'Toy','Abstract','Fixation'})
ax = gca; ax.XAxis.FontSize = 12; ax.YAxis.FontSize = 12;
ylim([0 length(miniMADE_MOIfc_table.ID)])
grid on; ax.XGrid = 'off'; ax.YGrid = 'on';
title('Overall inclusion rates across datasets','FontSize',12)

sgtitle('Data quantity for functional connectivity (global dbWPLI)','FontSize',14)


saveas(DataQuantityFC, 'DataQuantity_percondition.png');
close(DataQuantityFC)

%% Global dbWPLI values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1) FC spectra
RCT1 = figure;
Colours = colormap(parula(height(miniMADE_MOIfc_table)));
fig_linewidth = .5;
freqs = miniMADE_MOIfc_table.Sfreqs{1};

spS = subplot(2,2,1); % social 
plot([0,32],[0,0]);
hold on
for ss = 1:height(miniMADE_MOIfc_table.Seps)
    if miniMADE_MOIfc_table.Seps(ss) >= 90
        plot(freqs, miniMADE_MOIfc_table.SGlobdbWPLI{ss},'Color', Colours(ss,:),'LineStyle','-','LineWidth', fig_linewidth)
    end
end
xlim([1 32]); xlabel('Frequency (Hz)'); ylabel('Global dbWPLI')
title(strcat('Social: N=', num2str(sum((miniMADE_MOIfc_table.Seps >= 90),1))))

spT = subplot(2,2,2); % toy 
plot([0,32],[0,0]);
hold on
for ss = 1:height(miniMADE_MOIfc_table.Teps)
    if miniMADE_MOIfc_table.Teps(ss) >= 90
        plot(freqs, miniMADE_MOIfc_table.TGlobdbWPLI{ss},'Color', Colours(ss,:),'LineStyle','-','LineWidth', fig_linewidth)
    end
end
xlim([1 32]); xlabel('Frequency (Hz)'); ylabel('Global dbWPLI')
title(strcat('Toy: N=', num2str(sum((miniMADE_MOIfc_table.Teps >= 90),1))))

spA = subplot(2,2,3); % abstract 
plot([0,32],[0,0]);
hold on
for ss = 1:height(miniMADE_MOIfc_table.Aeps)
    if miniMADE_MOIfc_table.Aeps(ss) >= 90
        plot(freqs, miniMADE_MOIfc_table.AGlobdbWPLI{ss},'Color', Colours(ss,:),'LineStyle','-','LineWidth', fig_linewidth)
    end
end
xlim([1 32]); xlabel('Frequency (Hz)'); ylabel('Global dbWPLI')
title(strcat('Abstract: N=', num2str(sum((miniMADE_MOIfc_table.Aeps >= 90),1))))

spF = subplot(2,2,4); % fixation 
plot([0,32],[0,0]);
hold on
for ss = 1:height(miniMADE_MOIfc_table.Feps)
    if miniMADE_MOIfc_table.Feps(ss) >= 90
        plot(freqs, miniMADE_MOIfc_table.FGlobdbWPLI{ss},'Color', Colours(ss,:),'LineStyle','-','LineWidth', fig_linewidth)
    end
end
xlim([1 32]); xlabel('Frequency (Hz)'); ylabel('Global dbWPLI')
title(strcat('Fixation: N=', num2str(sum((miniMADE_MOIfc_table.Feps >= 90),1))))

% general
linkaxes([spS, spT, spA, spF], 'xy')
sgtitle('Functional connectivity spectra for individual datasets')


%% 2) Global dbWPLI values in FoI
% theta band: 6-7Hz, alpha band: 8-12Hz
addpath xxx/RainCloudPlots-master/tutorial_matlab

Freqs = miniMADE_MOIfc_table.Sfreqs{1};
Freq_boundaries = [6 7; 8 12];
Freqs_inds = [find(Freqs == Freq_boundaries(1,1)), find(Freqs == Freq_boundaries(1,2)); ...
    find(Freqs == Freq_boundaries(2,1)), find(Freqs == Freq_boundaries(2,2))];

RaincloudPlots = figure;
Colours = colormap(jet(size(Freqs_inds,1)));

spS = subplot(2,2,1); % social
    for ff = 1:size(Freqs_inds,1)
        FC_vals_all = cell2mat(miniMADE_MOIfc_table.SGlobdbWPLI(miniMADE_MOIfc_table.Seps >= 90)');
        FC_vals_all_subjxfreq = FC_vals_all';
        FC_vals = mean(FC_vals_all_subjxfreq(:,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
        raincloud_plot(FC_vals, 'box_on', 1, 'color', Colours(ff,:), 'alpha', 0.2,...
         'box_dodge', 1, 'box_dodge_amount', (0.15+(.2*(ff-1))), 'dot_dodge_amount', (0.15+(.2*(ff-1))),...
         'box_col_match', 0,'line_width',1);
        clear FC_vals FC_vals_all FC_vals_all_subjxfreq
    end
    % set(gca, 'YLim', [-20 30]); set(gca, 'XLim', [-.1 .2]);
    box off; view([-90 90]);
    title(strcat('Social: N=', num2str(sum((miniMADE_MOIfc_table.Seps >= 90),1))))
    xlabel('Global dbWPLI')

spT = subplot(2,2,2); % toy
    for ff = 1:size(Freqs_inds,1)
        FC_vals_all = cell2mat(miniMADE_MOIfc_table.TGlobdbWPLI(miniMADE_MOIfc_table.Teps >= 90)');
        FC_vals_all_subjxfreq = FC_vals_all';
        FC_vals = mean(FC_vals_all_subjxfreq(:,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
        raincloud_plot(FC_vals, 'box_on', 1, 'color', Colours(ff,:), 'alpha', 0.2,...
         'box_dodge', 1, 'box_dodge_amount', (0.15+(.2*(ff-1))), 'dot_dodge_amount', (0.15+(.2*(ff-1))),...
         'box_col_match', 0,'line_width',1);
        clear FC_vals FC_vals_all FC_vals_all_subjxfreq
    end
    % set(gca, 'YLim', [-10 20]); set(gca, 'XLim', [-.1 .2]);
    box off; view([-90 90]);
    title(strcat('Toy: N=', num2str(sum((miniMADE_MOIfc_table.Teps >= 90),1))))
    xlabel('Global dbWPLI')

spA = subplot(2,2,3); % abstract
    for ff = 1:size(Freqs_inds,1)
        FC_vals_all = cell2mat(miniMADE_MOIfc_table.AGlobdbWPLI(miniMADE_MOIfc_table.Aeps >= 90)');
        FC_vals_all_subjxfreq = FC_vals_all';
        FC_vals = mean(FC_vals_all_subjxfreq(:,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
        raincloud_plot(FC_vals, 'box_on', 1, 'color', Colours(ff,:), 'alpha', 0.2,...
         'box_dodge', 1, 'box_dodge_amount', (0.15+(.2*(ff-1))), 'dot_dodge_amount', (0.15+(.2*(ff-1))),...
         'box_col_match', 0,'line_width',1);
        clear FC_vals FC_vals_all FC_vals_all_subjxfreq
    end
    % set(gca, 'YLim', [-12 23]); set(gca, 'XLim', [-.1 .2]);
    box off; view([-90 90]);
    title(strcat('Abstract: N=', num2str(sum((miniMADE_MOIfc_table.Aeps >= 90),1))))
    xlabel('Global dbWPLI')

spF = subplot(2,2,4); % fixation
    for ff = 1:size(Freqs_inds,1)
        FC_vals_all = cell2mat(miniMADE_MOIfc_table.AGlobdbWPLI(miniMADE_MOIfc_table.Feps >= 90)');
        FC_vals_all_subjxfreq = FC_vals_all';
        FC_vals = mean(FC_vals_all_subjxfreq(:,Freqs_inds(ff,1):Freqs_inds(ff,2)),2);
        raincloud_plot(FC_vals, 'box_on', 1, 'color', Colours(ff,:), 'alpha', 0.2,...
         'box_dodge', 1, 'box_dodge_amount', (0.15+(.2*(ff-1))), 'dot_dodge_amount', (0.15+(.2*(ff-1))),...
         'box_col_match', 0,'line_width',1);
        clear FC_vals FC_vals_all FC_vals_all_subjxfreq
    end
    % set(gca, 'YLim', [-15 30]); set(gca, 'XLim', [-.1 .2]);
    box off; view([-90 90]);
    title(strcat('Fixation: N=', num2str(sum((miniMADE_MOIfc_table.Feps >= 90),1))))
    xlabel('Global dbWPLI')

% general
linkaxes([spS, spT, spA, spF], 'xy')
sgtitle('Global dbWPLI within theta & alpha bands for individual datasets')
