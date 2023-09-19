% this script will separate the subject-level spectCCEP file into single stim-pair file
% for speeding up loading files in parallel

home_dir = '/data/dian';
addpath(genpath(fullfile(home_dir, 'Dropbox/scripts/external/cbrewer2/cbrewer2')));
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/personal_validation'))
addpath('~/scripts/my_functions')
addpath(genpath('~/scripts/Stanford/ThalamocoricalLoop-project'))
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/CECS_Pipeline/CCEP_pipeline/HFB plot/subfunctions'));

data_folder = '~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore5_locked';
result_folder =  fullfile(data_folder, 'spectCCEP');
if ~exist(result_folder,'dir'); mkdir(result_folder);end
% plot_folder = '~/Dropbox/Stanford_Matters/data/THAL/CCEP/Plots/explore5/render2-reduced';
% if ~exist(plot_folder, 'dir'); mkdir(plot_folder);end
 if ~exist(fullfile(result_folder, 'SingleTXT'), 'dir')
        mkdir(fullfile(result_folder, 'SingleTXT')); end
% load metaTable
metaT = readtable('~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore5_locked/table_CCEPnewpipOutput_wholebrain_anatomical_info_activationRedone2.csv');
%
% load data
load('~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore5_locked/CCEP_all_flat_meanTr_cleaned.mat');%, 'CCEP_flat' each row corresponding to the metaT
%
sblist = unique(metaT.subject);
vars = metaT.Properties.VariableNames;
ttpIdx = find(contains(vars, 'pks_time'));
peakTimeMat = table2array(metaT(:,ttpIdx));

%filterIdx = load('~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore3_new2Sbjs/CCEP_all_flat_meanTr_cleaned.mat', 'badChan');
%badChan = badChan_new;
goodChan = ones(size(metaT,1),1);
goodChan(badChan) = 0;
cd(fullfile(data_folder, 'spectCCEP'));

tic
for is = 1:length(sblist)

    sb = sblist{is};

    %% step 1: get data and filter out bad channels

    power = load(sprintf('spectCCEP_power_%s.mat', sb));
    phase = load(sprintf('spectCCEP_phase_%s.mat', sb));
    itpc = load(sprintf('spectCCEP_itpc_%s.mat', sb));

    prefilterIdx =  find(goodChan & ...
        ~ismember(metaT.JP_label_in, {'', 'empty', 'NAN', 'NA'}) & ...
        ~ismember(metaT.JP_label_out, {'', 'empty', 'NAN', 'NA'}))  ;

    [filterIdx, preIdxsub, ~] = intersect(phase.idx_in_metaT, prefilterIdx, 'stable');

    subT = metaT(filterIdx,:);

    powerDat = power.CCEP_power(preIdxsub,:,:);
    phaseDat = phase.CCEP_phase(preIdxsub,:,:);
    itpcDat  = itpc.CCEP_itpc(preIdxsub,:,:);
    filteridx_metaT = filterIdx;

    %% step 2: transform data to get near-normal distribution
    toflip=[];
       
       % if ~exist(['flipCCEP_' sb '.mat'], 'file')
            [phaseflipped, toflip] = ...
            flipPhase(phaseDat, zccep_clean(filterIdx,:), peakTimeMat(filterIdx,:), toflip);
            %======================================
            save(['flipCCEP_' sb '.mat'], 'toflip');
            %======================================
       % else
          %  load(['flipCCEP_' sb '.mat'], 'toflip')
           % phaseflipped = ...
          %  flipPhase(phaseDat, [], [], toflip);
        %end

    % log transform to power
    pw = log10(abs(powerDat));
    % use sign flipped phase
    ph = phaseflipped;
    % sqrt transform to itpc
    pc = sqrt(itpcDat);

    
    %% save data in other formats
    col_idx = ismember(vars,...
        {'aSubID','JP_label_out','JP_label_in', 'Yeo7_out1', 'Yeo7_in1','eudDist','activated'});
    T = metaT(filteridx_metaT, col_idx);
    T.Yeo7_out1 = YeoNetShort(T.Yeo7_out1);
    T.Yeo7_in1  = YeoNetShort(T.Yeo7_in1);
    T.JP_label_in(strcmpi(T.JP_label_in, 'INSULA/LFC')) = {'INSULAorLFC'};
    T.JP_label_out(strcmpi(T.JP_label_out, 'INSULA/LFC')) = {'INSULAorLFC'};

    parfor i = 1:length(filteridx_metaT)
        labelname = sprintf('from%s-%s_to%s-%s_%dmm_act%d_%s_IDX%d.txt',...
            T.JP_label_out{i}, T.Yeo7_out1{i}, ...
            T.JP_label_in{i},  T.Yeo7_in1{i}, ...
            round(T.eudDist(i)), T.activated(i), ...
            T.aSubID{i}(1:3), filteridx_metaT(i)); %#ok<PFBNS>
        % write text file
        writematrix(Vrpw(i,:), fullfile(result_folder, 'SingleTXT' , ['pw_' labelname]),'Delimiter','tab')
        writematrix(Vrph(i,:), fullfile(result_folder, 'SingleTXT' , ['ph_' labelname]),'Delimiter','tab')
        writematrix(Vrpc(i,:), fullfile(result_folder, 'SingleTXT' , ['pc_' labelname]),'Delimiter','tab')
    end

end % ~ 70s per round
toc