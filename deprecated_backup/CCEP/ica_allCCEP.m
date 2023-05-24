% THIS SCRIPT WILL PRODUCE THE ICA TO ALL CCEP DATA OF THE THAL COHORT

home_dir = '/data/dian';

addpath(fullfile(home_dir, 'Dropbox/scripts/external/ColorBrewer'));
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/lbcn_preproc/vizualization'));
addpath(genpath(fullfile(home_dir, 'Dropbox/scripts/external/pca_ica')));
addpath(fullfile(home_dir, 'Dropbox/scripts/my_functions'));

%% load data
comp_root = '/data/dian/Working_projects/data';
result_folder = fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore1');
prelim_plot_folder = fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/CCEP/Plots/explore1');

d1 = load(fullfile(result_folder, 'CCEP_all_flat_meanTr.mat'));
d2 = load(fullfile(result_folder, 'CCEP_all_flat_medianTr.mat'));
metaT = readtable(fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/COHORT/table_CCEPnewpipOutput_wholebrain_anatomical_info.csv'));

%% perform ICA
% separate between subcor-subcor, subcor-cor, cor-subcor, cor-cor
figure; n = 1;
N = 4;
for type = {'subcor-subcor', 'subcor-cor', 'cor-subcor', 'cor-cor'}
    ty = type{1};
    switch ty
        case 'subcor-subcor'
            filterIdx = find( ...
                ismember(metaT.JP_label_in1, {'THALAMUS ANT','THALAMUS MID' ,'THALAMUS MID' ,'BG' }) & ... % record
                ismember(metaT.JP_label_out1, {'THALAMUS ANT','THALAMUS MID' ,'THALAMUS MID' ,'BG' })); % stim
        case 'subcor-cor'
            filterIdx = find( ...
                ~ismember(metaT.JP_label_in1, {'THALAMUS ANT','THALAMUS MID' ,'THALAMUS MID' ,'BG' }) & ...
                ismember(metaT.JP_label_out1, {'THALAMUS ANT','THALAMUS MID' ,'THALAMUS MID' ,'BG' }));
        case 'cor-subcor'
            filterIdx = find( ...
                ismember(metaT.JP_label_in1, {'THALAMUS ANT','THALAMUS MID' ,'THALAMUS MID' ,'BG' }) & ...
                ~ismember(metaT.JP_label_out1, {'THALAMUS ANT','THALAMUS MID' ,'THALAMUS MID' ,'BG' }));
        case 'cor-cor'
            filterIdx = find( ...
                ~ismember(metaT.JP_label_in1, {'THALAMUS ANT','THALAMUS MID' ,'THALAMUS MID' ,'BG' }) & ...
                ~ismember(metaT.JP_label_out1, {'THALAMUS ANT','THALAMUS MID' ,'THALAMUS MID' ,'BG' }));
    end
    disp(['Calculating ICA for ' ty ' pairs...'])
    ccep = d1.CCEP_flat(filterIdx,:);
    %perform rica 
    Mdl = rica(ccep, N, 'IterationLimit', 10000);
    ICs = Mdl.TransformWeights;
     %% --- compare algorithms ---
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
   %% % plot
    subplot(2,2,n); n=n+1;
    plot(ICs); title(type); legend(cellstr([repmat('IC',N,1) num2str((1:N)')]))
    save(fullfile(result_folder, sprintf('IC_CCEPmeanTr_%s.mat',ty)), 'ICs')
end
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

%% quality check the entropy of each CCEP,
% get rid of the invalid stimulation paris which may show up as active
% use second-level derirative
CCEP_flat_diffdiff = diff(diff(zscore(CCEP_flat', [],1)));
CCEP_flat2_diffdiff = diff(diff(zscore(CCEP_flat2', [], 1)));
dmax = max(abs(CCEP_flat_diffdiff),[],1, 'omitnan');

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
