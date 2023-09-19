% This code is inherited from prelim_CCEP_explore5_2.m of the SELF-project
%dl577@stanford.edu; Mar 7, 2023
% use newly preprocessed data
clear; close all
home_dir = getuserdir; %'/home/dian';%

server_root = '/Volumes/neurology_jparvizi$';

addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/Cohort_Organization'))
addpath(fullfile(home_dir, 'Dropbox/scripts/my_functions'));
addpath(fullfile(home_dir, 'Dropbox/scripts/external'));
addpath(genpath(fullfile(home_dir, 'Dropbox/scripts/external/cbrewer2/cbrewer2')));
addpath(fullfile(home_dir, 'Dropbox/scripts/external/ColorBrewer'));
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/lbcn_preproc/vizualization'));
addpath(fullfile(home_dir, 'Dropbox/scripts/external/pca_ica'));
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/CECS_Pipeline/CCEP_pipeline/tools'));
addpath(fullfile(home_dir, 'Dropbox/scripts/JAKE_2/lbcn_preproc/freesurfer'))
addpath(fullfile(home_dir, 'Dropbox/scripts/external/tight_subplot'))
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/personal_validation'))
addpath(fullfile(home_dir, 'Dropbox/scripts/PlottingTools/pipeline'))
addpath(fullfile(home_dir, 'Dropbox/scripts/PlottingTools/Brain_3D'))
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/ThalamocoricalLoop-project/CCEP'))

comp_root = '/data/dian/Working_projects/data';
result_folder = fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore2_redoPreProc_dataSort_Decomposition');
prelim_plot_folder = fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/CCEP/Plots/explore2');

%% load stimulation metatable
metaT = readtable(fullfile(result_folder, 'table_CCEPnewpipOutput_wholebrain_anatomical_info_activationRedone.csv'));
metaT.JP_label_out = ListSortAnatLabel_THAL(metaT.JP_label_out, 1);
metaT.JP_label_in = ListSortAnatLabel_THAL(metaT.JP_label_in, 1);
%{
sublist = {'S22_178_AF', 'S21_172_KS','S21_166_TM', 'S20_152_HT', 'S21_167_MQ', 'S21_169_BH', 'S21_171_MM',...
    'S19_137_AF', 'S21_165_WN', 'S22_176_LB', 'S22_177_JM', 'S22_182_DH', 'S22_183_CR', ...
    'S22_185_TW', 'S22_188_CB', 'S22_189_LMA', 'S22_190_AS', 'S22_191_KM', 'S22_192_LG', 'S22_193_AM', 'S23_194_PS', 'S23_195_MZ'};%, , 'S22_181_CB'};
%}
sublist = unique(metaT.aSubID)';

%%
time = -199:500;
%%
if ~exist(result_folder, 'dir'); mkdir(result_folder); end
if ~exist(prelim_plot_folder, 'dir'); mkdir(prelim_plot_folder); end
cd(result_folder)
% Make sure your are connected to CISCO and logged in the server

%% what analyses

to_plot_regional_CCEP  = 0;
plot_brain             = 1;

%% Load metatable
% assume metaT has been cleaned up (JP_label sorted)
% exclude cross border stim/rec pairs

vnames = metaT.Properties.VariableNames;
ttpIdx = find(contains(vnames, 'pks_time'));

%% load cleaned ccep data
load(fullfile(result_folder, 'CCEP_all_flat_meanTr_cleaned.mat'));%, 'CCEP_flat' each row corresponding to the metaT

%% assume flatten_CCEP step has been done
%%   ==================================== PLOT REGION AVERATED CCEP ====================================
%% loop over self hot and self cold elecgtrodes
if to_plot_regional_CCEP == 1
    %[badChan, zccep_clean]   = detectCCEP_badChan(CCEP_flat, 1);
    plot_regional_CCEP(zccep_clean, metaT, {'antTH','midTH','pstTH'}, [], prelim_plot_folder);
    
end

%% plot brain
% plot inflow pathway
if plot_brain == 1
    tic
    %    parpool(4)
    %-----------------------------------------
    output_folder = fullfile(prelim_plot_folder, 'render3');
    
    if ~exist(output_folder, 'dir'); mkdir(output_folder); end
    
    artPeak_time = 10;
    art_time = [-10:artPeak_time, 800:1999];
    onset = 200;
    %% ==================PREPARE DATA==================
    zccep_Tr = zccep_clean; clear zccep_clean%smoothdata(zccep_clean, 2, 'movmean', 10, 'omitnan');
    
    % coordinates
    MNIout = [metaT.MNIout_coord_1 metaT.MNIout_coord_2 metaT.MNIout_coord_3];
    % filter data
    goodChan = ~(isnan(mean(zccep_Tr,2, 'omitnan')));
    
    prefilterIdx =  goodChan &...
        metaT.MNIout_coord_1 .* metaT.MNIin_coord_1 > 0 &... % ipsilateral
        metaT.eudDist>5 & ...
        (sum(table2array(metaT(:, ttpIdx)), 2, 'omitnan') ) > artPeak_time   &... % exclude case when the only peak happens ealier than 10 ms
        ~ismember(metaT.JP_label_in, {'', 'empty', 'NAN', 'NA'}) & ...
        ~ismember(metaT.JP_label_out, {'', 'empty', 'NAN', 'NA'}) & ...
        metaT.sCrossBorder == 0  ;
    
    toantTH = strcmp(metaT.JP_label_in , 'antTH') & prefilterIdx;
    tomidTH = strcmp(metaT.JP_label_in , 'midTH') & prefilterIdx;
    topstTH = strcmp(metaT.JP_label_in , 'pstTH') & prefilterIdx;
    
    %------------------------------------------------------------
    toROIs = {toantTH, tomidTH, topstTH};
    AllStim = cellstr(join([string(metaT.aSubID), string(metaT.stim_chan)]));
    %------
    timepoints = -5:600;
    weights_cell = cell(length(toROIs),1);
    coordinates_cell = cell(length(toROIs),1);
    jplabel_cell = cell(length(toROIs),1);
    %------
    for t = 1:length(toROIs)
        to = find(toROIs{t});
        [stimlist, it] = unique(AllStim(to));
        idx = to(it);
        zccep_roi = zeros(length(stimlist),size(zccep_Tr,2));
        for i = 1:length(stimlist)
            keys = strsplit(stimlist{i}, ' ');
            asub = keys{1};
            stim = keys{2};
            % find the pairs with same stim ROI and rec channels
            commonROI = intersect(to, find(strcmp(metaT.stim_chan, stim) &...
                strcmp(metaT.aSubID, asub)));
            zccep_abs = abs(zccep_Tr(commonROI,:));
            % threshold
            zccep_abs(zccep_abs<2.2) = 0;
            zccep_abs(:, art_time+onset) = 0;
            % use mean of the absolute value
            zccep_roi(i,:) = mean(zccep_abs, 1, 'omitnan');
        end
        weights_cell{t} = zccep_roi(:, timepoints+onset) ;
        jplabel_cell{t} = metaT.JP_label_out(idx);
        coordinates_cell{t} = MNIout(idx,:);
    end
    anot = {'To antTH', 'To midTH', 'To pstTH'};
    % load brain surface data
    % brain surface data
    freesurferhome = getenv('FREESURFER_HOME');
    [cmcortex.vert, cmcortex.tri] = read_surf(fullfile(freesurferhome, 'subjects/fsaverage' , 'surf', 'rh.pial'));
    pers_rule = readtable( fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/COHORT/ELECTRODE_VISUALIZATION_RULE_thal.csv'));
    %clear up workspace to save memory
    clear zccep_Tr metaT zccep_clean
    
    %% plot surface weights
    clear cfg
    
    cfg.brainAlpha        = 0.78;
    cfg.suppress_colorbar = 1;
    cfg.gsp               = 70; % gaussian speading parameter
    cfg.figPosition       = [201   225   582   610];
    cfg.maxActRadius      = 350;
    cfg.anotFontSize      = 15;
    cfg.titleFontSize     = 14;
    cfg.titletext         = 'Inflow';
    cfg.recElecDotSize    = 25;
    cfg.covgAdj           = 1;
    cfg.climsPerc         = repmat([0 0.96], 3,1);
    cfg.climMin           = 0;
    
    switch cfg.titletext
        case 'Outflow'
            cfg.MarkerColor = repmat({[166,27,41]/256},3,1);
            cfg.cmaps = cbrewer2('seq','RdPu',64, 'spline');
            %litup_colors = [[166,27,41]/256; [255, 203, 1]/256;  [35, 118, 183]/256];% dark red[160, 28, 52],yellow, blue,
        case 'Inflow'
            cfg.MarkerColor = repmat({[207, 72, 19]/256},3,1);
            cfg.cmaps = cbrewer2('seq','OrRd',64, 'spline');
            %litup_colors = [[207, 72, 19]/256; [255, 203, 1]/256;  [35, 118, 183]/256];% dark red[160, 28, 52],yellow, blue,
    end
    %}
    % decide clims
    cfg.clims = repmat([0 8], 3, 1);
    %{
% the 0.97% quantile of the vals are about 7.2
    for iv = 1:length(weights_cell)

        vals = weights_cell{iv};

        if max(vals,[], 'all', 'omitnan') > 0 && min(vals,[],'all', 'omitnan') < 0 % divergent cmap
            % clim will consider all timepoints all rec chans
            val_max = max(quantile(abs(vals(vals<0)), cfg.climsPerc(iv,2), 'all'),...
                quantile(vals(vals>=0), cfg.climsPerc(iv,2), 'all'));
            cfg.clims(iv,:) = [-1.*val_max, val_max];
        else
            cfg.clims(iv,:) = [quantile(vals(vals>=0), cfg.climsPerc(iv,1), 'all'), ...
                quantile(vals(vals>=0), cfg.climsPerc(iv,2), 'all')];
        end

    end
    %}
    %--------------------------------------------------------------------------------------------------
  %  output_folder = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plotBrainSurfaceWeight(cmcortex , coordinates_cell, weights_cell, timepoints, jplabel_cell, output_folder, anot, pers_rule, cfg);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %{
    % if use parfor loop (not recommended)
for iv = 1:length(weights_cell)

    vals = weights_cell{iv};

        if max(vals,[], 'all', 'omitnan') > 0 && min(vals,[],'all', 'omitnan') < 0 % divergent cmap
            % clim will consider all timepoints all rec chans
            val_max = max(quantile(abs(vals(vals<0)), cfg.climsPerc(iv,2), 'all'),...
                quantile(vals(vals>=0), cfg.climsPerc(iv,2), 'all'));
            cfg.clims(iv,:) = [-1.*val_max, val_max];
        else
            cfg.clims(iv,:) = [quantile(vals(vals>=0), cfg.climsPerc(iv,1), 'all'), ...
                quantile(vals(vals>=0), cfg.climsPerc(iv,2), 'all')];
        end
   
end
    %--------------------------------------------------------------------------------------------------
    antccep = smoothdata(weights_cell{1}, 2, 'movmean', 10, 'omitnan');
    midccep = smoothdata(weights_cell{2}, 2, 'movmean', 10, 'omitnan');
    pstccep = smoothdata(weights_cell{3}, 2, 'movmean', 10, 'omitnan');
    clear weights_cell
    parfor it = 1:length(timepoints)
        weights_t = cell(3,1);
        timepoint = timepoints(it);
        weights_t{1} = antccep(:, it);
        weights_t{2} = midccep(:, it);
        weights_t{3} = pstccep(:, it);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        plotBrainSurfaceWeight(cmcortex , coordinates_cell, weights_t, timepoint, jplabel_cell, output_folder, anot, pers_rule, cfg);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    toc
    %}
end
