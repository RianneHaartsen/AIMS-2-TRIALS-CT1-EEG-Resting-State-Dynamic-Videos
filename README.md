Analyses scripts for AIMS-2-TRIALS RCT1 EEG resting state and social/ non-social videos

This folder contains the EEG pre-processing scripts used to extract spectral power, 1/f slope and offset, and connectivity for EEG data in the randomised controlled trial 1 examining the effects of the GABA-B agonist arbaclofen as part of the AIMS-2-TRIALS. The protocol of the study has been published in (Parellada et al., 2021). Hypotheses for the EEG metrics and plans for the pre-processed were pre-registered: (Del Bianco et al., 2023). The pre-processing scripts are based on the miniMADE pipeline with appropriate adjustments for the calculation of our metrics of interest (Troller-Renfree et al., 2021) and artefact rejection steps from (Mason et al., 2022)
	This folder also contains the R-studio script used for the statistical analyses in the study. 
	Details on the specific EEG preprocessing steps and statistical approach can be found in the main manuscript. The following provides a brief overview of the scripts. 

<img width="1008" height="958" alt="Overview of flow of scripts" src="https://github.com/user-attachments/assets/840d6071-f5fa-4ebd-847c-ab2dc2470a16" />

 
Scripts included:

Data check (all metrics): 
RCT1_00_DataSavedCheck.m
This script checks what EEG datasets have been collected and if this matches the log overview files. 

Initial preprocessing (all metrics):
RCT1_01_PreprocFilterContinuous.m
This script reads in the EEG datasets and ensures it is in the right EEGlab format for further preprocessing. It also applies a filter to the whole dataset. 

After this script, the pipeline forks into 3 different streams for the 3 different metrics:

A) Spectral power (log, absolute, and relative power)
RCT1_01A_PreprocPower_TimesSeries.m
This script segments the data during the tasks of interest (abstract video, fixation cross, social video, and toy video) into 1-second epochs with 50% overlap. It then cleans the data resulting in cleaned time series for each epoch for each channel. 
RCT1_02A_SpectralPower.m
This script takes the clean time series from the RCT1_01A script and calculates spectral power. 
RCT1_02Aa_VisualInspection.m
This script contains code to visualise the spectral power values from the datasets. 
RCT1_03A_Power_LongformatTables_withClinicalID_Regions_upto45Hz.m
This script takes the output from the RCT1_02A script and extracts the values of interest from power for the statistical analyses (per frequency band, per region, per condition for log, absolute, and relative power). Data are saved into long format tables for further analyses in RStudio. 

B) 1/f features and aperiodic-adjusted power
RCT1_01B_PreprocOOF_TimesSeries.m
This script segments the data during the tasks of interest (abstract video, fixation cross, social video, and toy video) into 2-second epochs with 50% overlap. It then cleans the data resulting in cleaned time series for each epoch for each channel. 
RCT1_02B_SpectralPower4FOOOF.m
This script takes the clean time series from the RCT1_01B script and calculates spectral power. It then saves the power spectrum in a .csv file for the FOOOF analyses in the next script. 
RCT1_02Ba_FOOOF_Region_NarrowFrequencyRange.py
This script contains code to extract 1/f features using Fitting of One Over F (FOOOF). 
RCT1_03B_OOF_LongformatTables_withClinicalID_Regions.m
This script takes the output from the RCT1_02B and RCT1_02Ba scripts and extracts the values of interest for the statistical analyses of 1/f (1/f slope and offset, R2 and error of fit of the FOOOF model per region, per condition). Data are saved into long format tables for further analyses in RStudio. 
RCT1_03Ba_ApAdjPow_LongformatTables_withClinicalID_Regions.m
This script takes the output from the RCT1_02B and RCT1_02Ba scripts and calculates aperiodic-adjusted power for the statistical analyses (based on 1/f slope and offset and absolute power per frequency band, per region, per condition). Data are saved into long format tables for further analyses in RStudio. 

C) Functional connectivity
RCT1_01C_PreprocFC_TimesSeries.m
This script segments the data during the tasks of interest (abstract video, fixation cross, social video, and toy video) into 1-second epochs with 50% overlap. It then cleans the data resulting in cleaned time series for each epoch for each channel. 
RCT1_02C_FunctionalConnectivity.m
This script takes the clean time series from the RCT1_01C script and estimates the debiased weighted phase lag index within each condition. This script calls to the function PLIbasic_and_directional2.m which extracts the dbWPLI. 
RCT1_02Ca_VisualInspection.m
This script contains code to visualise the connectivity values from the datasets. 
RCT1_03C_FC_LongformatTables_withClinicalID_Regions.m
This script takes the output from the RCT1_02C script and extracts the values of interest from the connectivity matrices for the statistical analyses (per frequency band, per region, per condition). Data are saved into long format tables for further analyses in RStudio. 

After running the scripts for the 3 pipeline streams, the data are then analysed in R studio using the script RCT1_04_Statistical Analyses. 

Other packages
The scripts call onto functions from TaskEngine (Jones et al., 2019), EEGlab (Delorme & Makeig, 2004), Fieldtrip(Oostenveld et al., 2011), and FOOOF (Donoghue et al., 2020), among others. See scripts for more details and the citations for the code. 




 
References
Del Bianco, T., Haartsen, R., & Jones, E. J. H. (2023). A Phase II Randomised, Double Blind, Placebo-Controlled Study of Effects of Arbaclofen on Periodic, Aperiodic and Event-Related Measures of EEG in Autistic Children and Adolescents: Pre-Registration of Processing and Analysis Plan. Open Science Framework.
Delorme, A., & Makeig, S. (2004). EEGLAB: an open source toolbox for analysis of single-trial EEG dynamics including independent component analysis. Journal of Neuroscience Methods, 134, 9–21. https://doi.org/10.1007/3-540-35375-5
Donoghue, T., Haller, M., Peterson, E. J., Varma, P., Sebastian, P., Gao, R., Noto, T., Lara, A. H., Wallis, J. D., Knight, R. T., Shestyuk, A., & Voytek, B. (2020). Parameterizing neural power spectra into periodic and aperiodic components. Nature Neuroscience, 23(December), 1655–1665. https://doi.org/10.1038/s41593-020-00744-x
Jones, E. J. H., Mason, L., Begum Ali, J., van den Boomen, C., Braukmann, R., Cauvet, E., Demurie, E., Hessels, R. S., Ward, E. K., Hunnius, S., Bolte, S., Tomalski, P., Kemner, C., Warreyn, P., Roeyers, H., Buitelaar, J., Falck-Ytter, T., Charman, T., & Johnson, M. H. (2019). Eurosibs: Towards robust measurement of infant neurocognitive predictors of autism across Europe. Infant Behavior and Development, 57(August 2018), 101316. https://doi.org/10.1016/j.infbeh.2019.03.007
Mason, L., Moessnang, C., Chatham, C., Ham, L., Tillmann, J., Dumas, G., Ellis, C., Leblond, C. S., Cliquet, F., Bourgeron, T., Beckmann, C., Charman, T., Oakley, B., Banaschewski, T., Meyer-Lindenberg, A., Baron-Cohen, S., Bölte, S., Buitelaar, J. K., Durston, S., … Jones, E. J. H. (2022). Stratifying the autistic phenotype using electrophysiological indices of social perception. Science Translational Medicine, 14(658). https://doi.org/10.1126/scitranslmed.abf8987
Oostenveld, R., Fries, P., Maris, E., & Schoffelen, J.-M. (2011). FieldTrip: Open Source Software for Advanced Analysis of MEG, EEG, and Invasive Electrophysiological Data. Computational Intelligence and Neuroscience, 2011, 1–9. https://doi.org/10.1155/2011/156869
Parellada, M., San José Cáceres, A., Palmer, M., Delorme, R., Jones, E. J. H., Parr, J. R., Anagnostou, E., Murphy, D. G. M., Loth, E., Wang, P. P., Charman, T., Strydom, A., & Arango, C. (2021). A Phase II Randomized, Double-Blind, Placebo-Controlled Study of the Efficacy, Safety, and Tolerability of Arbaclofen Administered for the Treatment of Social Function in Children and Adolescents With Autism Spectrum Disorders: Study Protocol for AIMS-2-TRIALS-CT1. Frontiers in Psychiatry, 12(August), 1–13. https://doi.org/10.3389/fpsyt.2021.701729
Troller-Renfree, S. V., Morales, S., Leach, S. C., Bowers, M. E., Debnath, R., Fifer, W. P., Fox, N. A., & Noble, K. G. (2021). Feasibility of assessing brain activity using mobile, in-home collection of electroencephalography: methods and analysis. Developmental Psychobiology, 63(6). https://doi.org/10.1002/dev.22128
 
