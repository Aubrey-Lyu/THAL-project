clear
home_dir = '/data/dian';
dir_base = getenv('DATA');
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/Cohort_Organization'))
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/CECS_Pipeline/CCEP_pipeline/auto_pip'))
addpath(fullfile(home_dir, 'Dropbox/scripts/my_functions'))
%% data analysis stesp
gather_data = 0;
plot_data = 1;
%% sort directories
project_name = 'CCEP';
center = 'Stanford';

data_result_folder = ...
    fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore5_locked/spectCCEP');
if exist(data_result_folder, 'dir') == 0; mkdir(data_result_folder); end

metaT_folder =  ...
    fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore5_locked');

%% load metaT
metaT = readtable(fullfile(metaT_folder, 'table_CCEPnewpipOutput_wholebrain_anatomical_info2.csv'));
load(fullfile(metaT_folder, 'CCEP_all_flat_meanTr_cleaned.mat'), 'badChan')
goodChan = ones(size(metaT,1), 1);
goodChan(badChan) = 0;

sublist = unique(metaT.subject);
cd(metaT_folder)

%% gather data
if gather_data == 1
    for is = 1:length(sublist)

        subject = sublist{is};
        disp(['Processing for ' subject '...'])
        subIDX = find(strcmp(metaT.subject, subject));
        subT = metaT(subIDX,:);
        % save data for each subject (because of limited working memory)
        %     CCEP_power = nan(size(subT,1), 501, 59); % trial-averaged power, dim: elecpair * time * freq,
        %     CCEP_phase = nan(size(subT,1), 501, 59); % trial-averaged phase
        %     CCEP_itpc  = nan(size(subT,1), 501, 59); % inter-trial phase coherence
        CCEP_spectra = nan(size(subT,1), 59); % overall spectrum
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
        parfor i = 1:size(subT,1)
            sdat = getCCEPspectr(subT(i,:), dir_base);
            %         CCEP_power(i,:,:) = squeeze(mean(sdat.wave, 2, 'omitnan'))';
            %         CCEP_phase(i,:,:) = squeeze(mean(sdat.phase, 2, 'omitnan'))';
            %         CCEP_itpc(i,:,:)  = calcITPC(sdat.phase)';
            CCEP_spectra(i,:) = sdat.spectrum;
            idx_in_metaT(i)   = subIDX(i);
        end
        toc
        %% save data
        %     if ~exist(fullfile(data_result_folder, sprintf('spectCCEP_power_%s.mat', subject)), 'file')
        %         save(fullfile(data_result_folder, sprintf('spectCCEP_power_%s.mat', subject)),...
        %             'CCEP_power', 'idx_in_metaT', ...
        %             'time', 'fsample', 'wavelet_span', 'spectrum', 'freqs','-v7.3')
        %     end
        %
        %     if ~exist(fullfile(data_result_folder, sprintf('spectCCEP_phase_%s.mat', subject)), 'file')
        %         save(fullfile(data_result_folder, sprintf('spectCCEP_phase_%s.mat', subject)),...
        %             'CCEP_phase', 'idx_in_metaT', ...
        %             'time', 'fsample', 'wavelet_span', 'spectrum', 'freqs','-v7.3')
        %     end
        %
        %     if ~exist(fullfile(data_result_folder, sprintf('spectCCEP_itpc_%s.mat', subject)), 'file')
        %         save(fullfile(data_result_folder, sprintf('spectCCEP_itpc_%s.mat', subject)),...
        %             'CCEP_itpc', 'idx_in_metaT', ...
        %             'time', 'fsample', 'wavelet_span', 'spectrum', 'freqs', '-v7.3')
        %     end
        if ~exist(fullfile(data_result_folder, sprintf('spectCCEP_spectra_%s.mat', subject)), 'file')
            save(fullfile(data_result_folder, sprintf('spectCCEP_spectra_%s.mat', subject)),...
                'CCEP_spectra', 'idx_in_metaT', ...
                'time', 'fsample', 'wavelet_span',  'freqs','-v7.3')
        end
    end
end

%% plot data
if plot_data == 1
    all_spect = [];
    idx = [];
    for is = 1:length(sublist)
        subject = sublist{is};
        load(fullfile(data_result_folder, sprintf('spectCCEP_spectra_%s.mat', subject)));
        all_spect = [all_spect; CCEP_spectra];
        idx = [idx; idx_in_metaT];
    end
    prefilterIdx =  find(goodChan & ...
        metaT.eudDist>5 & ...
        ~ismember(metaT.JP_label_in, {'', 'empty', 'NAN', 'NA'}) & ...
        ~ismember(metaT.JP_label_out, {'', 'empty', 'NAN', 'NA'}) & ...
        metaT.sCrossBorder == 0 & metaT.rCrossBorder == 0)  ;

    filterIdx = intersect(idx, prefilterIdx, 'stable');

end

spectra = normalize(log10(all_spect'),1);
x = normalize((1:size(spectra,1))');
Res = spectra;
parfor i = 1:size(spectra,2)
    y = spectra(:,i);
    b = x\y;
    res = y - b.*x ;
    Res(:,i) = res;
end

fromThal = find(ismember(metaT.JP_label_out, {'antTH', 'midTH', 'pstTH'}));
fromCor = find(~ismember(metaT.JP_label_out, {'antTH', 'midTH', 'pstTH', 'CLT','BG','AMY'}));

Res1 = Res(:, intersect(fromThal, filterIdx));
Res2 = Res(:, intersect(fromCor, filterIdx));

close all
figure; plot(Res1); xticks(1:5:59);xticklabels(freqs(1:5:59));
figure; plot(Res2); xticks(1:5:59);xticklabels(freqs(1:5:59));