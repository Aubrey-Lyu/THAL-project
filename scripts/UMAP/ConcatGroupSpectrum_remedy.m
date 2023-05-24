clear
home_dir = '/data/dian';
dir_base = getenv('DATA');
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/Cohort_Organization'))
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/CECS_Pipeline/CCEP_pipeline/auto_pip'))
addpath(fullfile(home_dir, 'Dropbox/scripts/my_functions'))

%% sort directories
project_name = 'CCEP';
center = 'Stanford';

data_result_folder = ...
    fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore4_spectrum');
if exist(data_result_folder, 'dir') == 0; mkdir(data_result_folder); end

metaT_folder =  ...
    fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore3_new2Sbjs');

%% load metaT
metaT = readtable(fullfile(metaT_folder, 'table_CCEPnewpipOutput_wholebrain_anatomical_info_activationRedone.csv'));
sublist = {'S21_172_KS', 'S21_171_MM', 'S21_169_BH'};%unique(metaT.subject);

for is = 1:length(sublist)

    subject = sublist{is};
    disp(['Processing for ' subject '...'])
    subIDX = find(strcmp(metaT.subject, subject));
    subT = metaT(subIDX,:);
    % save data for each subject (because of limited working memory)
    CCEP_power = nan(size(subT,1), 501, 59); % trial-averaged power, dim: elecpair * time * freq,
    CCEP_phase = nan(size(subT,1), 501, 59); % trial-averaged phase
    CCEP_itpc  = nan(size(subT,1), 501, 59); % inter-trial phase coherence
    % get other parameters
    sdat = getCCEPspectr(subT(1,:), dir_base);
    time         = sdat.time;
    fsample      = sdat.fsample;
    wavelet_span = sdat.wavelet_span;
    spectrum     = sdat.spectrum;
    freqs        = sdat.freqs;

    %% main loop
    tic
    idx_in_metaT = zeros(size(subT,1),1);
    for i = 1:size(subT,1)
        sdat = getCCEPspectr(subT(i,:), dir_base);
        CCEP_power(i,:,:) = squeeze(mean(sdat.wave, 2, 'omitnan'))';
        CCEP_phase(i,:,:) = squeeze(mean(sdat.phase, 2, 'omitnan'))';
        CCEP_itpc(i,:,:)  = calcITPC(sdat.phase)';
        idx_in_metaT(i)   = subIDX(i);
    end
    toc
    %% save data
   % if ~exist(fullfile(data_result_folder, sprintf('spectCCEP_power_%s.mat', subject)), 'file')
        save(fullfile(data_result_folder, sprintf('spectCCEP_power_%s.mat', subject)),...
            'CCEP_power', 'idx_in_metaT', ...
            'time', 'fsample', 'wavelet_span', 'spectrum', 'freqs','-v7.3')
   % end

   % if ~exist(fullfile(data_result_folder, sprintf('spectCCEP_phase_%s.mat', subject)), 'file')
        save(fullfile(data_result_folder, sprintf('spectCCEP_phase_%s.mat', subject)),...
            'CCEP_phase', 'idx_in_metaT', ...
            'time', 'fsample', 'wavelet_span', 'spectrum', 'freqs','-v7.3')
   % end

   % if ~exist(fullfile(data_result_folder, sprintf('spectCCEP_itpc_%s.mat', subject)), 'file')
        save(fullfile(data_result_folder, sprintf('spectCCEP_itpc_%s.mat', subject)),...
            'CCEP_itpc', 'idx_in_metaT', ...
            'time', 'fsample', 'wavelet_span', 'spectrum', 'freqs', '-v7.3')
   % end
end