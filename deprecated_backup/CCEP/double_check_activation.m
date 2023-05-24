% double check activation detection algorithm
clear
home_dir = '/data/dian';%getuserdir; %'C:\Users\lvdia';

addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/personal_validation'))
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/Cohort_Organization'))
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/ThalamocoricalLoop-project/CCEP'))
addpath(fullfile(home_dir, 'Dropbox/scripts/my_functions'))

dir_base = fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL');

result_folder = fullfile(dir_base, 'CCEP' , 'results', 'explore2_redoPreProc_dataSort_Decomposition');
cd(result_folder)
plot_folder = fullfile(dir_base,'PLOTS');
if ~exist(plot_folder, 'dir'); mkdir(plot_folder); end

% load data
metaT = readtable(fullfile(dir_base, 'COHORT', 'table_CCEPnewpipOutput_wholebrain_anatomical_info.csv'));
load CCEP_all_flat_meanTr_cleaned

% use multiple detection rules
[activated_redo, act_loc] = redetect_activation(metaT, zccep_clean); %may take a few minutes
save(fullfile(result_folder, 'activation_redo.mat'), 'activated_redo', 'act_loc');
% read activated columne from the meta table
activated = metaT.activated;
CECs_activation = metaT.CECS_activation;
idx_mismatch = activated~=activated_redo;
idx_mismatch2 = CECs_activation~=activated_redo;
idx_mismatch3 = CECs_activation~=activated;
fprintf('\ndefault-act and simpRule-act has %.2f%% mismatch cases.', 100*length(find(idx_mismatch))/size(metaT,1));
fprintf('\nCEC-act and simpRule-act has %.2f%% mismatch cases.', 100*length(find(idx_mismatch2))/size(metaT,1));
fprintf('\nCEC-act and default-act has %.2f%% mismatch cases.', 100*length(find(idx_mismatch3))/size(metaT,1));
disp(' ')
%% more measures
idx_11_20 = activated == 1 & activated_redo == 0;
idx_10_21 = activated == 0 & activated_redo == 1;
idx_11_30 = activated == 1 & CECs_activation == 0;
idx_10_31 = activated == 0 & CECs_activation == 1;
idx_21_30 = activated_redo == 0 & CECs_activation == 1;
idx_20_31 = activated_redo == 1 & CECs_activation == 0;
fprintf('\ndefault-act=1 and simpRule-act=0 has %.2f%% cases.', 100*length(find(idx_11_20))/size(metaT,1));
fprintf('\ndefault-act=0 and simpRule-act=1 has %.2f%% cases.', 100*length(find(idx_10_21))/size(metaT,1));
fprintf('\ndefault-act=1 and CEC-act=0 has %.2f%% cases.', 100*length(find(idx_11_30))/size(metaT,1));
fprintf('\ndefault-act=0 and CEC-act=1 has %.2f%% cases.', 100*length(find(idx_10_31))/size(metaT,1));
fprintf('\nsimpRule-act=1 and CEC-act=0 has %.2f%% cases.', 100*length(find(idx_21_30))/size(metaT,1));
fprintf('\nsimpRule-act=0 and CEC-act=1 has %.2f%% cases.', 100*length(find(idx_20_31))/size(metaT,1));
disp(' ')
%% more measures more
idx_10_20_30 = activated == 0 & activated_redo == 0 & CECs_activation == 0;
idx_10_20_31 = activated == 0 & activated_redo == 0 & CECs_activation == 1;
idx_11_20_30 = activated == 1 & activated_redo == 0 & CECs_activation == 0;
idx_11_20_31 = activated == 1 & activated_redo == 0 & CECs_activation == 1;
fprintf('\nWhen simpRule-act=0, default-act=0 and CEC-act=0 has %.2f%% cases.', 100*length(find(idx_10_20_30))/size(metaT,1));
fprintf('\nWhen simpRule-act=0, default-act=0 and CEC-act=1 has %.2f%% cases.', 100*length(find(idx_10_20_31))/size(metaT,1));
fprintf('\nWhen simpRule-act=0, default-act=1 and CEC-act=0 has %.2f%% cases.', 100*length(find(idx_11_20_30))/size(metaT,1));
fprintf('\nWhen simpRule-act=0, default-act=1 and CEC-act=1 has %.2f%% cases.', 100*length(find(idx_11_20_31))/size(metaT,1));
idx_10_21_30 = activated == 0 & activated_redo == 1 & CECs_activation == 0;
idx_10_21_31 = activated == 0 & activated_redo == 1 & CECs_activation == 1;
idx_11_21_30 = activated == 1 & activated_redo == 1 & CECs_activation == 0;
idx_11_21_31 = activated == 1 & activated_redo == 1 & CECs_activation == 1;
fprintf('\nWhen simpRule-act=1, default-act=0 and CEC-act=0 has %.2f%% cases.', 100*length(find(idx_10_21_30))/size(metaT,1));
fprintf('\nWhen simpRule-act=1, default-act=0 and CEC-act=1 has %.2f%% cases.', 100*length(find(idx_10_21_31))/size(metaT,1));
fprintf('\nWhen simpRule-act=1, default-act=1 and CEC-act=0 has %.2f%% cases.', 100*length(find(idx_11_21_30))/size(metaT,1));
fprintf('\nWhen simpRule-act=1, default-act=1 and CEC-act=1 has %.2f%% cases.', 100*length(find(idx_11_21_31))/size(metaT,1));
disp(' ')

%% visualise mismatch cases
% preset exclude index
exclude = metaT.eudDist<=5 | strcmp(metaT.JP_label_in, 'EXCLUDE')|...
strcmp(metaT.JP_label_out, 'EXCLUDE');
close all
%{
for counter = 1:20

figure('Position', [ 233          57        1555         894]); 
sgtitle('ramdom select 10 zccep')
%>>> when activated_simRule = 1 and activated_default = 0;
filterIdx = find(activated == 0 & activated_redo == 1 & ~exclude);
n = 6; % random select 5 examples
reind = randsample(length(filterIdx), n);
ccep1 = zccep_clean(filterIdx(reind),:);
labels = cell(n,1);
for i = 1:n
    rn = reind(i);
    labels{i} = sprintf('%s (%s) --> %s (%s)_{%d}, dist=%dmm, CEC-act=%d',...
        metaT.stim_chan{filterIdx(rn)}, metaT.JP_label_out{filterIdx(rn)},...
        metaT.record_chan{filterIdx(rn)}, metaT.JP_label_in{filterIdx(rn)},...
        filterIdx(rn),...
        round(metaT.eudDist(filterIdx(rn))), metaT.CECS_activation(filterIdx(rn))...
        );
end
%------------------------------------------------
subplot(4,1,1)
plot(ccep1'); legend(labels, 'Location', 'northeastoutside')
ylim([-10 10])
title({'activated-simpRule = 1 & activated-default = 0'})
%fn = sprintf('doubleCheckActivation_%s,png', ...
%    strjoin(strrep(cellstr(num2str(filterIdx(reind)))',' ',''),'_'));
%saveas(gcf, fullfile(plot_folder, fn));
%% >>> when activated_simRule = 0 and activated_default = 1
filterIdx = find(activated == 1 & activated_redo == 0 & ~exclude);
reind = randsample(length(filterIdx), n);
ccep2 = zccep_clean(filterIdx(reind),:);
labels = cell(n,1);
for i = 1:n
    rn = reind(i);
    labels{i} = sprintf('%s (%s) --> %s (%s)_{%d}, dist=%dmm, CEC-act=%d',...
        metaT.stim_chan{filterIdx(rn)}, metaT.JP_label_out{filterIdx(rn)},...
        metaT.record_chan{filterIdx(rn)}, metaT.JP_label_in{filterIdx(rn)},...
        filterIdx(rn),...
        round(metaT.eudDist(filterIdx(rn))), metaT.CECS_activation(filterIdx(rn))...
        );
end
%------------------------------------------------
subplot(4,1,2)

plot(ccep2'); legend(labels, 'Location', 'northeastoutside')
ylim([-10 10])
title({'activated-simpRule = 0 & activated-default = 1'})
%fn = sprintf('doubleCheckActivation_%s,png', ...
%    strjoin(strrep(cellstr(num2str(filterIdx(reind)))',' ',''),'_'));
%saveas(gcf, fullfile(plot_folder, fn));

%% >>> when CECs_activation = 0 and activated_default = 1
filterIdx = find(CECs_activation == 0 & activated == 1 & ~exclude);
reind = randsample(length(filterIdx), n);
ccep3 = zccep_clean(filterIdx(reind),:);
labels = cell(n,1);
for i = 1:n
    rn = reind(i);
    labels{i} = sprintf('%s (%s) --> %s (%s)_{%d}, dist=%dmm, simpRule-act=%d',...
        metaT.stim_chan{filterIdx(rn)}, metaT.JP_label_out{filterIdx(rn)},...
        metaT.record_chan{filterIdx(rn)}, metaT.JP_label_in{filterIdx(rn)},...
        filterIdx(rn),...
        round(metaT.eudDist(filterIdx(rn))), activated_redo(filterIdx(rn))...
        );
end
%------------------------------------------------
subplot(4,1,3) 

plot(ccep3'); legend(labels, 'Location', 'northeastoutside')
ylim([-10 10])
title({'CECs-activation = 0 & activated-default = 1'})
%fn = sprintf('doubleCheckActivation_%s,png', ...
%    strjoin(strrep(cellstr(num2str(filterIdx(reind)))',' ',''),'_'));
%saveas(gcf, fullfile(plot_folder, fn));

%% >>> when CECs_activation = 1 and activated_default = 0
filterIdx = find(CECs_activation == 1 & activated == 0 & ~exclude);
reind = randsample(length(filterIdx), n);
ccep4 = zccep_clean(filterIdx(reind),:);
labels = cell(n,1);
for i = 1:n
    rn = reind(i);
    labels{i} = sprintf('%s (%s) --> %s (%s)_{%d}, dist=%dmm, simpRule-act=%d',...
        metaT.stim_chan{filterIdx(rn)}, metaT.JP_label_out{filterIdx(rn)},...
        metaT.record_chan{filterIdx(rn)}, metaT.JP_label_in{filterIdx(rn)},...
        filterIdx(rn),...
        round(metaT.eudDist(filterIdx(rn))), activated_redo(filterIdx(rn))...
        );
end
%------------------------------------------------
subplot(4,1,4)

plot(ccep4'); legend(labels, 'Location', 'northeastoutside')
ylim([-10 10])
title({'CECs-activation = 1 & activated-default = 0'})
%fn = sprintf('doubleCheckActivation_%s,png', ...
%    strjoin(strrep(cellstr(num2str(filterIdx(reind)))',' ',''),'_'));
fn = sprintf('doubleCheckActivation_example-%d.png', counter);
saveas(gcf, fullfile(plot_folder, fn));

end
%}

%% activation stratege: democratic voting
act_vote = 0+(sum([activated_redo, activated, CECs_activation], 2) >= 2) ;
% pass badChan filter
act_vote(ismember(find(act_vote==1), badChan)) = 0;
% add manual inspection
%===========================================================
act_vote = manualcheck_activation(act_vote, [], 1);%======== 0:interactive
%===========================================================

% default-act
evaluation = makeConfusionMat(activated, act_vote);
fprintf(...
    '\ndefault-act has error rate = %.2f, accuracy = %.2f, \nsensitivity = %.2f, specificity = %.2f, \nprecision = %.2f, false positive rate = %.2f, \ncorrelation with the consensus is %.2f\n', ...
    evaluation.ERR, evaluation.ACC, evaluation.SN, evaluation.SP, ...
    evaluation.PREC, evaluation.FPR, evaluation.MCC)
% simpRule-act
evaluation = makeConfusionMat(activated_redo, act_vote);
fprintf(...
    '\nsimpRule-act has error rate = %.2f, accuracy = %.2f, \nsensitivity = %.2f, specificity = %.2f, \nprecision = %.2f, false positive rate = %.2f, \ncorrelation with the consensus is %.2f\n', ...
    evaluation.ERR, evaluation.ACC, evaluation.SN, evaluation.SP, ...
    evaluation.PREC, evaluation.FPR, evaluation.MCC)
% CEC-act
evaluation = makeConfusionMat( CECs_activation, act_vote);
fprintf(...
    '\nCEC-act has error rate = %.2f, accuracy = %.2f, \nsensitivity = %.2f, specificity = %.2f, \nprecision = %.2f, false positive rate = %.2f, \ncorrelation with the consensus is %.2f\n', ...
    evaluation.ERR, evaluation.ACC, evaluation.SN, evaluation.SP, ...
    evaluation.PREC, evaluation.FPR, evaluation.MCC)

%% add the act_vote to meta table, and reorganize

varnames = metaT.Properties.VariableNames;
pkIdx = contains(varnames, 'pks_time_');
pktimes = table2array(metaT(:,pkIdx));
pks_time = cell(length(act_loc),1);
min_pk_time = nan(length(act_loc),1);

for i = 1:length(act_loc)
    if act_vote(i) == 1
        if activated_redo(i) == 1 && activated(i) == 1
            pks_time{i} = pktimes(i,:);
        elseif activated_redo(i) == 1 && activated(i) == 0
            pks_time{i} = act_loc{i};
        elseif activated_redo(i) == 0 && activated(i) == 1
            pks_time{i} = pktimes(i,:);
        else
            pks_time{i} = nan;
        end
        min_pk_time(i) = min(pks_time{i},[], 'omitnan');
    else
        pks_time{i} = nan;
    end
    
end

metaT.activated_default = activated;
metaT.activated_simpRule = activated_redo;
metaT.activated = act_vote;
metaT(:, pkIdx) = [];
metaT.pks_time = pks_time;
metaT.min_pk_time = min_pk_time;
writetable(metaT, fullfile(result_folder , 'table_CCEPnewpipOutput_wholebrain_anatomical_info_activationRedone.csv'));

