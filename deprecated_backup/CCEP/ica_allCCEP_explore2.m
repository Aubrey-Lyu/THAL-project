% THIS SCRIPT WILL PRODUCE THE ICA TO ALL CCEP DATA OF THE THAL COHORT
% The input data are partly imported from SELF-project, and partly from
% _CCEP_Behzad_Neda

%home_dir = '/data/dian';
home_dir = getuserdir;

addpath(fullfile(home_dir, 'Dropbox/scripts/external/ColorBrewer'));
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/lbcn_preproc/vizualization'));
addpath(genpath(fullfile(home_dir, 'Dropbox/scripts/external/pca_ica')));
addpath(fullfile(home_dir, 'Dropbox/scripts/my_functions'));
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/personal_validation'))
addpath(genpath('/media/lvdian/Dropbox/scripts/external/bpc_scripts/bpc_scripts/'))

%% load data
comp_root = '/data/dian/Working_projects/data';
result_folder = fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore2_dataSort_Denoise_TAadded');
prelim_plot_folder = fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/CCEP/Plots/explore2');

d1 = load(fullfile(result_folder, 'CCEP_all_flat_meanTr.mat'));
d2 = load(fullfile(result_folder, 'CCEP_all_flat_medianTr.mat'));
metaT = readtable(fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/COHORT/table_CCEPnewpipOutput_wholebrain_anatomical_info.csv'));

%% denoise CCEP
[badChan, zccep_Tr]   = detectCCEP_badChan(d1.CCEP_flat, 1);
[badChan2, zccep_Tr2] = detectCCEP_badChan(d2.CCEP_flat2, 1);

%% perform ICA
% use only ipsilateral data
ipsiPairIdx = metaT.MNIout_coord_1 .* metaT.MNIin_coord_1 > 0;
% remove stim-rec pairs closer than 5 mm
closeIdx = metaT.eudDist > 5;
% goodChan (not naned overall)
zccep_Tr = smoothdata(zccep_Tr', 'gaussian',15)';
goodChan = ~(isnan(mean(zccep_Tr,2, 'omitnan')));
% group the preselect index
preIdx = ipsiPairIdx & closeIdx & goodChan;

%% separate between subcor-subcor, subcor-cor, cor-subcor, cor-cor
close; figure('Position', [  1254         705        1271         621]);
n = 1;
N = 4;
for type = {'subcor-subcor', 'subcor-cor', 'cor-subcor', 'cor-cor'}
    ty = type{1};
    switch ty
        case 'subcor-subcor'
            filterIdx = preIdx & ( ...
                ismember(metaT.JP_label_in1, {'THALAMUS ANT','THALAMUS MID' ,'THALAMUS MID' ,'BG' }) & ... % record
                ismember(metaT.JP_label_out1, {'THALAMUS ANT','THALAMUS MID' ,'THALAMUS MID' ,'BG' })); % stim
            filterIdx1 = filterIdx;
        case 'subcor-cor'
           filterIdx = preIdx &( ...
                ~ismember(metaT.JP_label_in1, {'THALAMUS ANT','THALAMUS MID' ,'THALAMUS MID' ,'BG' }) & ...
                ismember(metaT.JP_label_out1, {'THALAMUS ANT','THALAMUS MID' ,'THALAMUS MID' ,'BG' }));
             fiterIdx2 = filterIdx;
        case 'cor-subcor'
            filterIdx = preIdx &( ...
                ismember(metaT.JP_label_in1, {'THALAMUS ANT','THALAMUS MID' ,'THALAMUS MID' ,'BG' }) & ...
                ~ismember(metaT.JP_label_out1, {'THALAMUS ANT','THALAMUS MID' ,'THALAMUS MID' ,'BG' }));
             filterIdx3 = filterIdx;
        case 'cor-cor'
            filterIdx = preIdx &( ...
                ~ismember(metaT.JP_label_in1, {'THALAMUS ANT','THALAMUS MID' ,'THALAMUS MID' ,'BG' }) & ...
                ~ismember(metaT.JP_label_out1, {'THALAMUS ANT','THALAMUS MID' ,'THALAMUS MID' ,'BG' }));
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

%% BPC to each categrories
% Trial-averaged level
sc1 = strcat(metaT.subject(filterIdx1),'.', metaT.sc1(filterIdx1));
rc1 = strcat(metaT.subject(filterIdx1),'.', metaT.rc1(filterIdx1));
Blocks = metaT.block_name(filterIdx1);
% load data
for ib = 1:length(Blocks)
    block = Blocks{ib};
    sbj = strsplit(sc1{ib}, '.');
    sbj = sbj{1};
    d = load(fullfile())
end

ccep = zccep_Tr(filterIdx1,:);
Cset = [sc1; rc1];
[C, ~, ic] = unique(Cset);
pair = reshape(ic, [] ,2);
P = table();
P.pair = pair;
P.indices = (1:length(pair))';

V = ccep(:,200:800)';

pair_types = table2struct(P);

[B,excluded_pairs]=bpc_identify(V,pair_types);

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
