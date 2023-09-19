home_dir = '/data/dian'; %'/home/lvdian'; %
% use reduced spect due to working memory limitation
addpath(genpath(fullfile(home_dir, 'Dropbox/scripts/external/cbrewer2/cbrewer2')));
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/personal_validation'))
addpath('~/scripts/my_functions')
addpath(genpath('~/scripts/Stanford/ThalamocoricalLoop-project'))
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/CECS_Pipeline/CCEP_pipeline/HFB plot/subfunctions'));

data_folder = '~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore5_locked/UMAP_learn/resample3';
metaT_folder = '~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore5_locked';
result_folder =  '~/Dropbox/Stanford_Matters/data/THAL/UMAP/ALLDATA_semisupervise';
umapact_folder = '~/Dropbox/Stanford_Matters/data/THAL/UMAP/WITHINSBJ_semisupervise';

metaT = readtable(fullfile(metaT_folder, 'table_CCEPnewpipOutput_wholebrain_anatomical_info_activationRedone2.csv'));

umapAct = nan(size(metaT,1),1);

%% exclusion label

sbj_list = unique(T.subject);

for ss = 1:length(sbj_list)
sbj = sbj_list{ss};
isact = textread(fullfile(umapact_folder, sprintf('actLearn_2clusters_anot_%s.txt', sbj)), '%f', 'delimiter', '\n');
isact(isact==0) = nan; % original codes: 0=noise, 1=inactive, 2=active
isact = isact-1;
% load index
I = load(fullfile(data_folder, sprintf('SpecReduceCollapse_%s.mat', sbj)), 'filteridx_metaT');
umapAct(I.filteridx_metaT) = isact;
end


%% activation label - refined by supervised learning
T = readtable(fullfile(result_folder, 'brainInfo.csv'));
features = readcell(fullfile(result_folder, 'cleanFeatures.txt'), Delimiter = "");
isnoise = textread(fullfile(result_folder, 'isnoise.txt'), '%d', 'delimiter', '\n');
T1 = T(isnoise==1,:);

idx_noise = T1.index + 1; % account for python calculator starting from 0
umapAct(idx_noise) = NaN; % new codes: nan=noise, 0 = inactive, 1=acive

save(fullfile(metaT_folder, 'umapAct.mat'), 'umapAct');

%% activation stratege: democratic voting
cd(metaT_folder)
diary UmapActivationEval.log

activated = metaT.activated_default;
activated_redo = metaT.activated_SimRule;
CECs_activation = metaT.CECS_activation;

act_vote = 0+(sum([activated_redo, activated, CECs_activation], 2) >= 2) ;

% default-act
evaluation = makeConfusionMat(activated, umapAct);
fprintf(...
    '\ndefault-act has error rate = %.2f, accuracy = %.2f, \nsensitivity = %.2f, specificity = %.2f, \nprecision = %.2f, false positive rate = %.2f, \ncorrelation with the observation is %.2f\n', ...
    evaluation.ERR, evaluation.ACC, evaluation.SN, evaluation.SP, ...
    evaluation.PREC, evaluation.FPR, evaluation.MCC)
% simpRule-act
evaluation = makeConfusionMat(activated_redo, umapAct);
fprintf(...
    '\nsimpRule-act has error rate = %.2f, accuracy = %.2f, \nsensitivity = %.2f, specificity = %.2f, \nprecision = %.2f, false positive rate = %.2f, \ncorrelation with the observation is %.2f\n', ...
    evaluation.ERR, evaluation.ACC, evaluation.SN, evaluation.SP, ...
    evaluation.PREC, evaluation.FPR, evaluation.MCC)
% CEC-act
evaluation = makeConfusionMat( CECs_activation, umapAct);
fprintf(...
    '\nCEC-act has error rate = %.2f, accuracy = %.2f, \nsensitivity = %.2f, specificity = %.2f, \nprecision = %.2f, false positive rate = %.2f, \ncorrelation with the observation is %.2f\n', ...
    evaluation.ERR, evaluation.ACC, evaluation.SN, evaluation.SP, ...
    evaluation.PREC, evaluation.FPR, evaluation.MCC)
% act_vote
evaluation = makeConfusionMat(act_vote, umapAct);
fprintf(...
    '\nact_consensus has error rate = %.2f, accuracy = %.2f, \nsensitivity = %.2f, specificity = %.2f, \nprecision = %.2f, false positive rate = %.2f, \ncorrelation with the observation is %.2f\n', ...
    evaluation.ERR, evaluation.ACC, evaluation.SN, evaluation.SP, ...
    evaluation.PREC, evaluation.FPR, evaluation.MCC)
% save report
diary off
