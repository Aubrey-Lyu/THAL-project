% THIS SCRIPT WILL PRODUCE THE ICA TO ALL CCEP DATA OF THE THAL COHORT
% The input data are partly imported from SELF-project, and partly from
% _CCEP_Behzad_Neda
clear
%home_dir = '/data/dian';
home_dir = getuserdir;

addpath(fullfile(home_dir, 'Dropbox/scripts/external/ColorBrewer'));
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/lbcn_preproc/vizualization'));
addpath(genpath(fullfile(home_dir, 'Dropbox/scripts/external/pca_ica')));
addpath(fullfile(home_dir, 'Dropbox/scripts/my_functions'));
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/personal_validation'))
addpath(genpath(fullfile(home_dir,'Dropbox/scripts/external/bpc_scripts/bpc_scripts/')))

%% load data
comp_root = '/data/dian/Working_projects/data';
result_folder = fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore1_dataSort_Denoise');
prelim_plot_folder = fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/CCEP/Plots/explore1');

d1 = load(fullfile(result_folder, 'CCEP_all_flat_meanTr.mat'));
d2 = load(fullfile(result_folder, 'CCEP_all_flat_medianTr.mat'));
metaT = readtable(fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/COHORT/table_CCEPnewpipOutput_wholebrain_anatomical_info_importedData.csv'));

%% denoise CCEP
[badChan, zccep_Tr]   = detectCCEP_badChan(d1.CCEP_flat, 1); clear d1
[badChan2, zccep_Tr2] = detectCCEP_badChan(d2.CCEP_flat2, 1); clear d2

%% perform ICA to subregions of the brain
% use only ipsilateral data
ipsiPairIdx = metaT.MNIout_coord_1 .* metaT.MNIin_coord_1 > 0;
% remove stim-rec pairs closer than 5 mm
distIdx = metaT.eudDist > 5;
% goodChan (not naned overall)
zccep_Tr = smoothdata(zccep_Tr', 'gaussian',15)';
goodChan = ~(isnan(mean(zccep_Tr,2, 'omitnan')));
% group the preselect index
preIdx = ipsiPairIdx & distIdx & goodChan;

%% separate between subcor-subcor, subcor-cor, cor-subcor, cor-cor
figure('Position', [  1254         705        1271         621]);
n = 1;
N = 4;
%{'subcor-subcor', 'subcor-cor', 'cor-subcor', 'cor-cor'}
SUBCOR_GROUP = {'THALAMUS ANT','THALAMUS MID' ,'THALAMUS POST' ,'BG', 'AMY', 'CLAUSTRUM' };
%{'thal-cor', 'thal-thal', 'cor-thal', 'cor-cor'}
THAL_GROUP = {'THALAMUS ANT','THALAMUS MID' ,'THALAMUS POST'};
COR1_GROUP = {'FG','MCC','OCC','SOMATOMOTOR','SOMATOSENSORY','STG'};
COR2_GROUP = {'HPC','HPC ANT','HPC MID','HPC POST','INSULA','INSULA/LFC',...
    'IPL','ITG','LFC','MFC','MTG','OFC','PMC','SPL','TP'};
for type = {'THAL-COR2', 'COR2-THAL', 'THAL-THAL', 'COR2-COR2'}
    ty = type{1};
    switch ty
        case 'THAL-COR2'
            filterIdx = preIdx & ( ...
                ismember(metaT.JP_label_in1, COR2_GROUP) & ... % record
                ismember(metaT.JP_label_out1, THAL_GROUP)); % stim
            filterIdx1 = filterIdx;
        case 'COR2-THAL'
            filterIdx = preIdx &( ...
                ismember(metaT.JP_label_in1, THAL_GROUP) & ...
                ismember(metaT.JP_label_out1, COR2_GROUP));
            fiterIdx2 = filterIdx;
        case 'THAL-THAL'
            filterIdx = preIdx &( ...
                ismember(metaT.JP_label_in1, THAL_GROUP) & ...
                ismember(metaT.JP_label_out1, THAL_GROUP));
            filterIdx3 = filterIdx;
        case 'COR2-COR2'
            filterIdx = preIdx &( ...
                ismember(metaT.JP_label_in1, COR2_GROUP) & ...
                ismember(metaT.JP_label_out1, COR2_GROUP));
            filterIdx4 = filterIdx;
    end
    SampleSize = sum(filterIdx);
    disp(['Calculating ICA for ' ty ' pairs...'])
    ccep = zccep_Tr(filterIdx,:);
    %perform rica
    Mdl = rica(ccep, N, 'IterationLimit', 10000);
    ICs = Mdl.TransformWeights;
    % --- compare algorithms ---
    %{
    N = 5;
    cm = cbrewer2('Set2',N);
%----------------------------------
% in case of out of memory
%ccep = ccep(:,200:800);
%perform rica 
    Mdl = rica(ccep, N, 'IterationLimit', 10000);
    ICs1 = Mdl.TransformWeights;
    %-----------------------------
    % perform fast ICA
    Zfica = fastICA(ccep, N);
    %-----------------------------
    % Perform max-kurtosis ICA
    Zkica = kICA(ccep, N);
    %-----------------------------
    % Perform PCA
    Zpca = PCA(ccep, N);
    %-----------------------------
    figure;
    % Fast ICA
    subplot(4,1,1);
    for i = 1:N
        plot(Zfica(i,:),'-','Color',cm(i,:)); hold on;
    end
    title('Independent components [Fast ICA]');
    axis tight;

    % Max-kurtosis
    subplot(4,1,2);
    for i = 1:N
        plot(Zkica(i,:),'-','Color',cm(i,:)); hold on;
    end
    title('Independent components [max-kurtosis]');
    axis tight;

    % PCA
    subplot(4,1,3);
    for i = 1:N
        plot(Zpca(i,:),'-','Color',cm(i,:)); hold on;
    end
    title('Principal components');
    axis tight;
    % rica
    subplot(4,1,4);
    for i = 1:N
        plot(ICs(:,i),'-','Color',cm(i,:)); hold on;
    end
    title('rica');
    axis tight;
    %}
    % % plot
    subplot(2,2,n); n=n+1;
    plot(ICs(1:1500,:)); title(sprintf('%s (N = %d)', ty, SampleSize));
    legend(cellstr([repmat('IC',N,1) num2str((1:N)')]))
    save(fullfile(result_folder, sprintf('IC_CCEPmeanTr_%s.mat',ty)), 'ICs')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BPC to each categrories
%{'subcor-subcor', 'subcor-cor', 'cor-subcor', 'cor-cor'}
SUBCOR_GROUP = {'THALAMUS ANT','THALAMUS MID' ,'THALAMUS POST' ,'BG', 'AMY', 'CLAUSTRUM' };
%{'thal-cor', 'thal-thal', 'cor-thal', 'cor-cor'}
THAL_GROUP = {'THALAMUS ANT','THALAMUS MID' ,'THALAMUS POST'};
COR1_GROUP = {'FG','MCC','OCC','SOMATOMOTOR','SOMATOSENSORY','STG'};
COR2_GROUP = {'HPC','HPC ANT','HPC MID','HPC POST','INSULA','INSULA/LFC',...
    'IPL','ITG','LFC','MFC','MTG','OFC','PMC','SPL','TP'};
type = {'THAL-COR', 'COR-THAL', 'THAL-THAL', 'COR-COR'};
for t = 4:length(type)
    ty = type{t};
    switch ty
        case 'THAL-COR'
            filterIdx = preIdx & ( ...
                ~ismember(metaT.JP_label_in1, SUBCOR_GROUP) & ... % record
                ismember(metaT.JP_label_out1, THAL_GROUP)); % stim
            filterIdx1 = filterIdx;
        case 'COR-THAL'
            filterIdx = preIdx &( ...
                ismember(metaT.JP_label_in1, THAL_GROUP) & ...
                ~ismember(metaT.JP_label_out1, SUBCOR_GROUP));
            fiterIdx2 = filterIdx;
        case 'THAL-THAL'
            filterIdx = preIdx &( ...
                ismember(metaT.JP_label_in1, THAL_GROUP) & ...
                ismember(metaT.JP_label_out1, THAL_GROUP));
            filterIdx3 = filterIdx;
        case 'COR-COR'
            filterIdx = preIdx &( ...
                ~ismember(metaT.JP_label_in1, SUBCOR_GROUP) & ...
                ~ismember(metaT.JP_label_out1, SUBCOR_GROUP));
            filterIdx4 = filterIdx;
    end
    % Trial-averaged level
    
    sc = strcat(metaT.subject(filterIdx),'.', metaT.stim_chan(filterIdx));
    rc = strcat(metaT.subject(filterIdx),'.', metaT.record_chan(filterIdx));

    Blocks = metaT.block_name(filterIdx);
    vnames = metaT.Properties.VariableNames;
    ALL_reject_trials = table2array(metaT(filterIdx,contains(vnames, 'reject_trials')));
    %------- prepare input for BPC function %-------
    Cset = [sc; rc];
    % ic is a unique number/ID for each elec
    [C, ~, ic] = unique(Cset); % mind the order, reshape(C(ic),[],2) is in the same order as metaT(filterIdx,:), or reshape(Cset,[],2), or P.pair
    pair = reshape(ic, [] ,2);
    pair_types_ = table();
    pair_types_.pair = pair;
    block_loaded = 'dummy';
    NTr = nan(height(pair_types_),1);
    V_ = nan(2201,55,length(Blocks));% all timeseries: Timepoints*trial(50)*block
    tic
    % organize trial-level data
    parfor ib = 1:length(Blocks) % same order as P

        block = Blocks{ib};

        sbj = strsplit(sc{ib}, '.'); sbj = sbj{1};
        ccep_dir = fullfile(comp_root, 'computed_data/CCEP', sbj);
        rec_Chan = strsplit(rc{ib}, '.'); rec_Chan = rec_Chan{2};
        % load data
        % if ~strcmp(block_loaded, block)
        d = load(fullfile(ccep_dir, [block '_CCEP.mat'])); %d.CCEP
        block_loaded = block;
        % end

        rec_idx = strcmp(d.CCEP.recordingchannel, rec_Chan);
        wave = squeeze(d.CCEP.wave(:,rec_idx,:))';

        %% exclude bad trials
        reject_trials = ALL_reject_trials(ib,:)';
        reject_trials(isnan(reject_trials)) = [];

        wv_meanT  = mean(abs(wave),1, 'omitnan');
        thr_data   = mean(wv_meanT, 'omitnan') + 4.*std(wv_meanT, 'omitnan');
        reject_trials  = [reject_trials;...
            find(wv_meanT >= thr_data)';...
            find(isnan(wv_meanT))'];
        % [rb, cb]       = find(wave>100); % containing any extreme values bigger than 100 z-scores
        % reject_trials  = [reject_trials; unique(cb)];
        reject_trials  = unique(reject_trials);
        reject_trials(reject_trials>size(wave,2))=[];
        % exclude bad trials
        wave(:,reject_trials) = nan;
        %% put data together
        nTr = size(wave,2)-length(reject_trials);
        NTr(ib) = nTr;
        wave_ = nan(2201,55); wave_(:,1:size(wave,2)) = wave(1:2201,:);
        V_(:,:,ib) = wave_;
    end
    toc % ~3min

    % ----> > > reorganize the data matrix >>>
    V = reshape(V_, 2201, []);
    V(:,isnan(mean(V,1,'omitnan')))=[];
    %----------------------------------
    indices = cell(length(NTr),1);
    parfor ib = 1:length(NTr)
        NTr0 = [0; NTr];
        nTr = NTr(ib);
        from = sum(NTr0(1:ib))+1;
        to = sum(NTr0(1:ib))+nTr;
        indices{ib} = from:to;
    end

    pair_types_.indices = indices;
    pair_types = table2struct(pair_types_);

    % save data
    subfolder = 'BPC';
    outputDir = fullfile(result_folder, subfolder);
    if ~exist(outputDir, 'dir'); mkdir(outputDir); end
    %ty = 'THAL-COR2';
    save(fullfile(outputDir, ['BPC_input_' ty '.mat']), 'pair_types', 'V', 'NTr', 'filterIdx', '-v7.3');
    %{
%% run BPC for group
[B,excluded_pairs] = bpc_identify(V, pair_types);
% save results
save(fullfile(outputDir, sprintf('BPC_output_%s_group.mat', ty), 'B','excluded_pairs'));
    %}
    %% run BPC for each subject
    tic
    %load(fullfile(outputDir, ['BPC_input_' ty '.mat']));

    sbjs = unique(metaT.subject(filterIdx));
    for ss = 1:length(sbjs)
        sbj = sbjs{ss};
        sidx = strcmp(metaT.subject(filterIdx), sbj);
        pair_types_ = pair_types(sidx);
        NTr_ = NTr(sidx);
        v = [];
        indices_new = cell(length(pair_types_),1);
        ntr = 0;
        for n = 1:length(pair_types_)
            indices = pair_types_(n).indices;
            ntr = length(indices);
            indices_new{n} = reshape((size(v,2)+1) : (size(v,2)+ntr),[],1);
            v = [v, V(:, indices)];

        end
        pair_types_sbj = struct2table(pair_types_);
        pair_types_sbj.indices = indices_new;
        pair_types_sbj = table2struct(pair_types_sbj);
        if length(pair_types_sbj) < 2
            continue;
        end
        %run bpc
        [B,excluded_pairs] = bpc_identify(v, pair_types_sbj);
        % save results
        if ~exist(fullfile(outputDir, ty), 'dir'); mkdir(fullfile(outputDir, ty)); end
        save(fullfile(outputDir, ty, sprintf('BPC_output_%s_%s.mat', sbj, ty)), 'B','excluded_pairs');
    end
end
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ICA to whole brain
N = 5;
%---------------------meanTr-------------------------
Mdl = rica(d1.CCEP_flat, N, 'IterationLimit', 100000);
ICs = Mdl.TransformWeights;
save(fullfile(result_folder, 'IC_CCEPmeanTr_flat.mat'), 'ICs') % IC(:,3) is noise
close; figure; plot(ICs); legend
%--------------------medianTr-------------------------
Mdl = rica(d2.CCEP_flat2, N, 'IterationLimit', 100000);
ICs = Mdl.TransformWeights;
save(fullfile(result_folder, 'IC_CCEPmedianTr_flat.mat'), 'ICs') % IC(:,2) is noise
close; figure; plot(ICs); legend

% define noise IC
%load(fullfile(result_folder, 'IC_CCEPmedianTr_flat.mat'))
%noiseIC = ICs(:,2);
%save(fullfile(result_folder, 'noiseIC.mat'), 'noiseIC')
load (fullfile(result_folder, 'noiseIC.mat'))

%% regress out the noise IC component from the IC
% for meanTrial CCEP
CCEP_flat = d1.CCEP_flat;
parfor i = 1:size(CCEP_flat,1)
    ccep1 = CCEP_flat(i,:);
    gl = fitlm(noiseIC, ccep1);
    yhat = predict(gl, noiseIC);
    res = ccep1 - yhat';
    CCEP_flat(i,:) = res;
end
save(fullfile(result_folder, 'CCEP_all_flat_meanTr_scrubbed.mat'), 'CCEP_flat');

% for medianTrial CCEP
CCEP_flat2 = d2.CCEP_flat2;
parfor i = 1:size(CCEP_flat2,1)
    ccep1 = CCEP_flat2(i,:);
    gl = fitlm(noiseIC, ccep1);
    yhat = predict(gl, noiseIC);
    res = ccep1 - yhat';
    CCEP_flat2(i,:) = res;
end
save(fullfile(result_folder, 'CCEP_all_flat_medianTr_scrubbed.mat'), 'CCEP_flat2');

%% apply smoothing to scrubbed CCEP
load(fullfile(result_folder, 'CCEP_all_flat_meanTr_scrubbed.mat'))
load(fullfile(result_folder, 'CCEP_all_flat_medianTr_scrubbed.mat'))
CCEP_flat  = smoothdata(CCEP_flat', 'gaussian',25)';
CCEP_flat2 = smoothdata(CCEP_flat2', 'gaussian',35)';


%% plot for quality control
% comparing four options:
filterIdx = find( metaT.eudDist>5 & ...
    metaT.min_pk_time > 15 &... % those that have been detected active
    ~ismember(metaT.JP_label_in1, {'', 'empty', 'NAN', 'EXCLUDE', 'NA'}) & ...
    ~ismember(metaT.JP_label_out1, {'', 'empty', 'NAN', 'EXCLUDE', 'NA'}));

ylims = [-50 30];
%---- plot -----
Irand = randi([1,length(filterIdx)],1,200);% randomly exhibit 200 CCEP to show
Idx = filterIdx(Irand);
close; figure;

subplot(2,2,1)
dat = d1.CCEP_flat(Idx,:)';
% baseline = dat(1:190,:); unit = std(baseline, 1, 'omitnan') ;
% dat = dat./unit; % normalize accoridng to baseline
plot(dat)
%ylim(ylims)
title('CCEP (mean)')

subplot(2,2,2)
dat = CCEP_flat(Idx,:)';
% baseline = dat(1:190,:); unit = std(baseline, 1, 'omitnan') ;
% dat = dat./unit; % normalize accoridng to baseline
plot(dat)
%ylim(ylims)
title({'CCEP scrubbed (mean)', '+ smoothed (25 pts)'})

subplot(2,2,3)
dat = d2.CCEP_flat2(Idx,:)';
% baseline = dat(1:190,:); unit = std(baseline, 1, 'omitnan') ;
% dat = dat./unit; % normalize accoridng to baseline
plot(dat)
%ylim(ylims)
title('CCEP (median)')

subplot(2,2,4)
dat = CCEP_flat2(Idx,:)';
% baseline = dat(1:190,:); unit = std(baseline, 1, 'omitnan') ;
% dat = dat./unit; % normalize accoridng to baseline
plot(dat)
%ylim(ylims)
title({'CCEP scrubbed (median)','+ smoothed (35 pts)'})


%% plot PC components

%-------
time_onset  = -50;
time_offset = 450;
time = time_onset+1:time_offset;
time_loc = time_onset-(-200)+1:time_offset+200;
%-------

close
figure;
colors = brewermap(N,'Set1');

for i_c = 1:N
    [pk, loc] = findpeaks(ICs(time_loc,i_c), 'MinPeakHeight', 0.04, 'MinPeakDistance', 100);

    hold on
    if i_c == 2
        plot(time, ICs(time_loc, i_c), 'LineWidth', 1.5, 'Color', colors(i_c,:),'LineStyle','--', 'DisplayName', ['IC' num2str(i_c)])
    else
        plot(time, ICs(time_loc, i_c), 'LineWidth', 1.5, 'Color', colors(i_c,:), 'DisplayName', ['IC' num2str(i_c)])
    end
    %plot([time_onset time_offset],[0 0], 'LineWidth',1, 'Color','k')
    %plot([time(loc) time(loc)],[0 pk])
    scatter(time(loc), pk+0.003, 55, 'v', 'filled', 'MarkerFaceColor', colors(i_c,:))
    text(time(loc)-10, pk+0.01, sprintf('%.f ms',time(loc)))
end
% y_min = -0.06;
% y_max = 0.17;
% ylim([y_min y_max])
xlabel('Time (ms)'); ylabel('Learned transformation weights')
%legend({'IC1','IC2','IC3','IC4'})
% add areas
seg = [20 78; 78 135; 135 375];
colors = colors([1 4 3],:);

y_points = [y_min y_max y_max y_min];
for i_fill = 1:3
    x_points = [seg(i_fill, 1) seg(i_fill, 1) seg(i_fill, 2) seg(i_fill, 2)];
    color = colors(i_fill,:);
    a = fill(x_points, y_points, color, 'EdgeColor', 'none');
    a.FaceAlpha = 0.22;
end
% plot onset
plot([0 0],[y_min y_max], 'LineWidth',0.8, 'Color',[0.3 0.3 0.3])
