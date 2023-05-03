% test bipolar significance
clear
home_dir = getuserdir;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/lbcn_preproc/data_convertion'));
addpath(genpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/lbcn_preproc/cleaning')));
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/lbcn_preproc/vizualization'));
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/lbcn_preproc/data_processing'));
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/lbcn_preproc/lbcn_personal/personal'));
addpath(genpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY')));
addpath(genpath(fullfile(home_dir, 'Dropbox/scripts/external/eeglab')))
addpath(fullfile(home_dir, 'Dropbox/scripts/external/spm12'))
% initialize EEGlab
eeglab;
% for using python external packages:)
setenv('PYTHONPATH', '/home/dian/.conda/envs/mne_old/bin/:/data/dian/Dropbox/scripts/Stanford/CECS_Pipeline_COPY/SSD_code/SSD_code/pedro_zip/Python/Packages/')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[server_root, comp_root, code_root] = AddPaths('exo-lbcn', 'CCEP');
%comp_root = '/data/dian/Working_projects/data';
result_folder = fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore2_redoPreProc_dataSort_Decomposition');
plot_folder =  fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/PLOTS');

metaT = readtable(fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/COHORT/table_CCEPnewpipOutput_wholebrain_anatomical_info.csv'));
vnames = metaT.Properties.VariableNames;
load(fullfile(result_folder, 'CCEP_all_flat_meanTr_cleaned.mat'));
nearIdx  = metaT.eudDist > 5 & metaT.eudDist < 10;
farIdx = metaT.eudDist > 100;
actIdx = metaT.activated == 0;
block_loaded = 'dummy';
% near pairs
filterIdx = find(nearIdx & actIdx);
filterIdx(ismember(filterIdx, badChan)) = [];


%% random select 5 times
for n_rand = 1:5
    randSelect = randsample(length(filterIdx),1);
    % meanTr CCEP
    zccep_far = zccep_clean(filterIdx(randSelect),:);
    %figure; %plot(zccep_far')

    %% single trial CCEP
    block = metaT.block_name{filterIdx(randSelect)};
    stimpair = metaT.stim_chan{filterIdx(randSelect)};

    sbj_id = metaT.subject{filterIdx(randSelect)};

    rawdata_fn = fullfile(comp_root, 'neuralData/originalData', sbj_id, block, [block '.edf']);

    %% process original ccep without bipolar reference
    if strcmp(block, block_loaded) == 0 % do not rewrite this part
        dirs = InitializeDirs('CCEP', sbj_id, comp_root, server_root, code_root);
        if ismember(sbj_id, ...
                {'S22_178_AF', 'S21_172_KS','S21_166_TM', 'S20_152_HT', 'S21_167_MQ', 'S21_169_BH', 'S21_171_MM',...
    'S19_137_AF', 'S21_165_WN', 'S22_176_LB', 'S22_177_JM', 'S22_182_DH', 'S22_183_CR'}) % self cohort
            rawdata_fn = strrep(rawdata_fn,'/data/dian/Working_projects/', '/media/dian/wasabi-dl577/');
        end
        block_loaded = block;
        % onset detect
        cfg.dirs                = dirs;
        cfg.rsample             = 1000;
        cfg.sbj_name            = sbj_id;
        cfg.block_name          = block;
        cfg.validation_test     = 0;
        cfg.manual_select_elec  = 0;
        cfg.use_Trigger         = 1;
        %-------------
        if ismember(cfg.sbj_name, {'S19_137_AF'}) % may be true for other early pts
            cfg.interval        = 2000;
        else
            cfg.interval        = 2070;
        end
        %-------------
        % if autodetect failed, change manual_select_elec to 1 (written in the script)
        autoDetect_failedCases
        [eeg_bi, chanNames, Triggers, dat, original_chanNames] = BipolarNihonKohden_trigger(rawdata_fn,1); % second input to do %filter & notch
        thr = 4.5;
        [onset, cfg] = onsetDetect0(eeg_bi/1000, chanNames, stimpair, [], thr,Triggers, cfg);
        zccep0 = epochCCEP(onset, dat'/1000, original_chanNames, cfg);
    end
    %% Compare bipolar referenced zccep and original ccep
    % select rc1 and rc2
    rc1 = lower(metaT.rc1{filterIdx(randSelect)});
    rc2 = lower(metaT.rc2{filterIdx(randSelect)});
    wave = zccep0.wave(:, strcmpi(zccep0.recordingchannel, rc1)| ...
        strcmpi(zccep0.recordingchannel, rc2),:);

    % exclude bad trials
    reject_trials = table2array(metaT(filterIdx(randSelect),contains(vnames, 'reject_trials')))';
    reject_trials(isnan(reject_trials)) = [];

    wv_meanTr  = mean(abs(wave),[2 3], 'omitnan');
    thr_data   = mean(wv_meanTr, 'omitnan') + 4.*std(wv_meanTr, 'omitnan');
    reject_trials  = [reject_trials; find(wv_meanTr >= thr_data)'];
    % [rb, cb]       = find(wave>100); % containing any extreme values bigger than 100 z-scores
    % reject_trials  = [reject_trials; unique(cb)];
    reject_trials  = unique(reject_trials);

    % detected from the pipeline
    trial_No = 1:size(wave, 1);
    trial_No(ismember(trial_No, reject_trials)) = [];

    wave_clean = wave(trial_No,:,:);
    wave0 = squeeze(mean(wave_clean, 1, 'omitnan'));

    %% plot to compare
    close; figure('Position', [675   266   785   685])
    sgtitle({'Comparing original zCCEP vs bipolar-ref zCCEP (not activated)',...
        ['Eucidean Distance = ' num2str(metaT.eudDist(filterIdx(randSelect))) 'mm']})
    subplot(2,1,1)
    title('Original zCCEP (cleaned)')
    colors = [47 144 185; 147 181 207]/256;
    hold on
    for i = 1:size(wave0,1)
        plot(wave0(i,:)', 'Color', colors(i,:), 'LineWidth',1.5);
    end
    yline(0, '--k')
    hold off
    legend({['rc1:' rc1 ' of ' strrep(sbj_id, '_', '\_')], ['rc2:' rc2]})
    %------------------
    subplot(2,1,2)
    title('Bipolar-referenced zCCEP (cleaned)')
    colors =  [192 44 56; 249 215 112]/256;
    zccep_referenced = zccep_far;
    zccep0_contrast = (wave0(1,:)-wave0(2,:));
    d = [zccep_referenced; zccep0_contrast];
    hold on
    for i = 1:size(d,1)
        plot(d(i,:)', 'Color', colors(i,:), 'LineWidth',1.5);
    end
    yline(0, '--k')
    hold off
    legend({['bipoar-ref zccep:' metaT.record_chan{filterIdx(randSelect)}],...
        'origial zccep contrast (rc1-rc2)'})
    saveas(gcf, fullfile(plot_folder, ['bipolar-ref_effects_No' num2str(filterIdx(randSelect)), '.png']));
end