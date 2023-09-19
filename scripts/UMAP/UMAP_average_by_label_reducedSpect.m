home_dir = '/data/dian'; %'/home/lvdian'; %
% use reduced spect due to working memory limitation
addpath(genpath(fullfile(home_dir, 'Dropbox/scripts/external/cbrewer2/cbrewer2')));
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/personal_validation'))
addpath('~/scripts/my_functions')
addpath(genpath('~/scripts/Stanford/ThalamocoricalLoop-project'))
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/CECS_Pipeline/CCEP_pipeline/HFB plot/subfunctions'));

data_folder = '~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore5_locked/UMAP_learn/resample3';
result_folder =  '~/Dropbox/Stanford_Matters/data/THAL/UMAP/ALLDATA_semisupervise';
outfigDir = '~/Dropbox/Stanford_Matters/data/THAL/PLOTS/UMAP/supervisedLearn_Group';

%% load labels
T = readtable(fullfile(result_folder, 'brainInfo.csv'));
features = readcell(fullfile(result_folder, 'cleanFeatures.txt'), Delimiter = "");
isnoise = textread(fullfile(result_folder, 'isnoise.txt'), '%d', 'delimiter', '\n');
T0 = T(isnoise==0,:);
T1 = T(isnoise==1,:);
features0 = features(isnoise==0);

% check noise ratio for each subject
sblist = unique(T.aSubID);
%{
ratio = ones(length(sblist),1);
for sn = 1:length(sblist)
    sb = sblist{sn};
    nall = length(find(strcmp(T.aSubID, sb)));
    ngood = length(find(strcmp(T0.aSubID, sb)));
    ratio(sn) = ngood/nall;
end

close;figure;bar(categorical(strrep(sblist, '_', '\_')), ratio, 'FaceColor','none'), box off
%}

%% load data based on each clean feature categrories
Ftrs = unique(features0);
vectCCEP_featured = struct();

for i = 1:length(Ftrs)
    ftr = Ftrs{i};
    vectCCEP_featured(i).feature = ftr;

    fIdx = strcmp(features0, ftr);
    t0 = T0(fIdx,:);

    VRPW = nan(length(fIdx), 1380);
    VRPC  = nan(length(fIdx), 1200);
    filterIdx = [];
    tic
    for sn = 1:length(sblist)

        sIdx = find(strcmp(t0.aSubID, sblist{sn}));

        if isempty(sIdx)
            continue;
        end

        load(fullfile(data_folder, ['SpecReduceCollapse_' t0.subject{sIdx(1)}]));

        [filter_idx,~, row2select] = intersect(t0.index(sIdx)+1, filteridx_metaT,  'stable'); % +1 account for python emuerating starting from 0

        spect1 = Vrpw(row2select,:,:);
        spect2 = Vrpc(row2select,:,:);

        VRPW(sIdx,:) = spect1;
        VRPC(sIdx,:)  = spect2;

        filterIdx = [filterIdx; filter_idx];

    end
    VRPW(isnan(mean(VRPW, 2, 'omitnan')),:) = [];
    VRPC(isnan(mean(VRPC, 2, 'omitnan')),:) = [];
    vectCCEP_featured(i).Vrpw = VRPW;
    vectCCEP_featured(i).Vrpc = VRPC;
    vectCCEP_featured(i).idx_in_metaT = filterIdx;

end
toc

%% plot feature spectCCEP
% get time/frequency information
info = load(sprintf('/data/dian/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore5_locked/spectCCEP/spectCCEP_power_%s.mat', ...
    t0.subject{1}), 'time', 'freqs'); % example sbj
[~, rfreq, rtime] = reduceSpectrum([], info.freqs,info.time*1000,...
    [4, 4, 4, 4, 4, 3], [5 25 20 10],0);

freq = cell(length(rfreq),1);
for i = 1:length(rfreq)
    freq{i} = sprintf('%.1f', rfreq(i));
end

time = rtime/1000;


close all

figure('Position', [995         691        1151         620]);
for i = 1:4%length(vectCCEP_featured)
    subplot(2,2,i)
    feature = vectCCEP_featured(i).feature;
    meanspect1 = reshape(mean(vectCCEP_featured(i).Vrpw, 1, 'omitnan'), [],60);

    cmap = myColors('inferno',[],[],[],0);
    cmap = cmap.cm;
    labelstring = 'zscore-log(power)'; % 'zscore-log(power)', 'zscore-sqrt(itpc)'
    subplotCCEPspectrum(meanspect1, feature,time,freq,cmap,labelstring,5,3)
    caxis([-0.8, 1])
end
sgtitle('Power')
% save
figname = sprintf('umapfeature_Power2.png');
exportgraphics(gcf, fullfile(outfigDir, figname),'Resolution',300)

%----------------------------------------------------
figure('Position', [995         691        1151         620]);
for i = 1:length(vectCCEP_featured)
    subplot(2,2,i)
    feature = vectCCEP_featured(i).feature;
    meanspect2 = reshape(mean(vectCCEP_featured(i).Vrpc, 1, 'omitnan'), [], 60);

    cmap = myColors('viridis',[],[],[],0);
    cmap = cmap.cm;
    labelstring = 'zscore-sqrt(itpc)'; % 'zscore-log(power)', 'zscore-sqrt(itpc)'
    subplotCCEPspectrum(meanspect2, feature,time,freq(1:end-1),cmap,labelstring,5,3)
    caxis([-0.8, 1])
end
sgtitle('Inter-trial phase coherence')
% save
figname = sprintf('umapfeature_ITPC2.png');
%exportgraphics(gcf, fullfile(outfigDir, figname),'Resolution',300)

%% save data
%save(fullfile(result_folder, 'vectCCEP_featured.mat'), 'vectCCEP_featured', '-v7.3');

%% plot the extra cluster
%{
feature_plus = 'THAL-contr2';
clstLab = textread(fullfile(result_folder, 'manualClst.txt'), '%d', 'delimiter', '\n');
Tc = T(clstLab==1,:);

vectCCEP_featured(5).feature = feature_plus;

fIdx = strcmp(features0, ftr);
t0 = Tc;

VRPW = nan(length(fIdx), 1380);
VRPC  = nan(length(fIdx), 1200);
filterIdx = [];

for sn = 1:length(sblist)

    sIdx = find(strcmp(t0.aSubID, sblist{sn}));

    if isempty(sIdx)
        continue;
    end

    load(fullfile(data_folder, ['SpecReduceCollapse_' t0.subject{sIdx(1)}]));

    [filter_idx,~, row2select] = intersect(t0.index(sIdx)+1, filteridx_metaT,  'stable'); % +1 account for python emuerating starting from 0

    spect1 = Vrpw(row2select,:,:);
    spect2 = Vrpc(row2select,:,:);

    VRPW(sIdx,:) = spect1;
    VRPC(sIdx,:)  = spect2;

    filterIdx = [filterIdx; filter_idx];

end
vectCCEP_featured(5).Vrpw = VRPW;
vectCCEP_featured(5).Vrpc = VRPC;
vectCCEP_featured(5).idx_in_metaT = filterIdx;
%--------------- plot extra cluster-------------------------------------
figure;

    subplot(2,1,1)
    feature = vectCCEP_featured(5).feature;
    meanspect1 = reshape(mean(vectCCEP_featured(5).Vrpw, 1, 'omitnan'), [],60);

    cmap = myColors('inferno',[],[],[],0);
    cmap = cmap.cm;
    labelstring = 'zscore-log(power)'; % 'zscore-log(power)', 'zscore-sqrt(itpc)'
    subplotCCEPspectrum(meanspect1, feature,time,freq,cmap,labelstring,5,3)


    subplot(2,1,2)
    feature = vectCCEP_featured(5).feature;
    meanspect2 = reshape(mean(vectCCEP_featured(5).Vrpc, 1, 'omitnan'), [], 60);

    cmap = myColors('viridis',[],[],[],0);
    cmap = cmap.cm;
    labelstring = 'zscore-sqrt(itpc)'; % 'zscore-log(power)', 'zscore-sqrt(itpc)'
    subplotCCEPspectrum(meanspect2, feature,time,freq(1:end-1),cmap,labelstring,5,3)

% save
figname = sprintf('umapfeature_extraClst.png');
%exportgraphics(gcf, fullfile(outfigDir, figname),'Resolution',300)
%}




