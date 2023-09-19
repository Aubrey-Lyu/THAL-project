% This code is inherited from prelim_CCEP_explore5_2.m of the SELF-project
%dl577@stanford.edu; Mar 7, 2023
% use newly preprocessed data
clear; close all
home_dir = getuserdir; %'/home/dian';

server_root = '/mnt/neurology_jparvizi$';

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
addpath(fullfile(home_dir, 'Dropbox/scripts/my_functions'))

%  brain surface manipulation tools
toolboxDir = '~/Dropbox/scripts/external/surfAnalysis-master';
addpath(genpath(toolboxDir))
addpath(genpath('~/Dropbox/scripts/external/matlab_GIfTI-master'))
addpath('~/Dropbox/scripts/external')
addpath(genpath('/Users/dl577/Dropbox/scripts/external/brainSurfer-main'))

comp_root = '/data/dian/Working_projects/data';
result_folder = fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore3_new2Sbjs');
prelim_plot_folder = fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/CCEP/Plots/explore3');

fsDir = getenv('SUBJECTS_DIR');
brainDir = '~/Dropbox/scripts/external/brainSurfer-main/brains';
wbDir = fullfile('/Volumes/Lacie', 'workbench');

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

to_plot_regional_CCEP   = 0;
plot_brain_peakCluster  = 1;

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
if plot_brain_peakCluster == 1
    
    tic
    %    parpool(4)
    %-----------------------------------------
    output_folder = fullfile(prelim_plot_folder);
    
    if ~exist(output_folder, 'dir'); mkdir(output_folder); end
    
    artPeak_time = 10;
    art_time = [-10:artPeak_time, 800:1999];
    onset = 200;
    paths = {'outflow', 'inflow'};
    clusters = {'early', 'late1', 'late2'};
    
    %% sort FS_LR standard coordinates
    FSout_coord = [metaT.LEPTOout_coord_1, metaT.LEPTOout_coord_2, metaT.LEPTOout_coord_3];
    FSin_coord  = [metaT.LEPTOin_coord_1, metaT.LEPTOin_coord_2, metaT.LEPTOin_coord_3];
    m = load(fullfile(wbDir,'affineMsurf2FSLR_group.mat'));
    
    for ss = 1:length(sublist)
        
        asub = sublist(ss);
        sbj = metaT.subject(strcmp(metaT.aSubID, asub));
        sbj = sbj{1}; keys = strsplit(sbj, '_');
        sbj_short = strjoin(keys(1:2),'_');
        
        lepto_cord_out = [metaT.LEPTOout_coord_1(strcmp(metaT.aSubID, asub),:),...
            metaT.LEPTOout_coord_2(strcmp(metaT.aSubID, asub),:),...
            metaT.LEPTOout_coord_3(strcmp(metaT.aSubID, asub),:)];
        
        lepto_cord_in = [metaT.LEPTOin_coord_1(strcmp(metaT.aSubID, asub),:),...
            metaT.LEPTOin_coord_2(strcmp(metaT.aSubID, asub),:),...
            metaT.LEPTOin_coord_3(strcmp(metaT.aSubID, asub),:)];
        % find correct affine matrix to use
        im = find(contains(m.sblist, sbj_short));
        Msurf2space = m.affineMsurf2FSLR_group{im};
        % project to the right hemisphere
        [~,~,sfs_cordsOUT_2] = convertLepto2FS_LR0(sbj_short, lepto_cord_out, Msurf2space,  wbDir, brainDir,999); % radius = 999 --> use all surface projectsions
        [~,~,sfs_cordsIN_2] = convertLepto2FS_LR0(sbj_short, lepto_cord_in, Msurf2space,  wbDir, brainDir,999);
        
        FSout_coord(strcmp(metaT.aSubID, asub),:) = sfs_cordsOUT_2;
        FSin_coord(strcmp(metaT.aSubID, asub),:) = sfs_cordsIN_2;
    end
    
    % substitude some coordinates with MNI;
    MNIin  = [metaT.MNIin_coord_1 metaT.MNIin_coord_2 metaT.MNIin_coord_3];
    MNIout = [metaT.MNIout_coord_1 metaT.MNIout_coord_2 metaT.MNIout_coord_3];
    
    FSout_coord_inter = FSout_coord;
    FSout_coord_inter(isnan(FSout_coord(:,1))&~strcmpi(metaT.JP_label_out, 'EXCLUDE')) = MNIout(isnan(FSout_coord(:,1))&~strcmpi(metaT.JP_label_out, 'EXCLUDE'),:);
    
    FSin_coord_inter = FSin_coord;
    FSin_coord_inter(isnan(FSin_coord(:,1))&~strcmpi(metaT.JP_label_in, 'EXCLUDE')) = MNIin(isnan(FSin_coord(:,1))&~strcmpi(metaT.JP_label_in, 'EXCLUDE'),:);
    
    metaT.FSout_coord = FSout_coord;
    metaT.FSin_coord = FSin_coord;
    
    %% ==================PREPARE DATA==================
    zccep_Tr = smoothdata(zccep_clean, 2, 'movmean', 10, 'omitnan');
    
    for p = 1:length(paths)
        path = paths{p};
        % coordinates
        if strcmp(path, 'outflow')
            cords = metaT.FSin_coord;
            %cords = [metaT.MNIin_coord_1 metaT.MNIin_coord_2 metaT.MNIin_coord_3];
        else
            cords = metaT.FSout_coord;
            %cords = [metaT.MNIout_coord_1 metaT.MNIout_coord_2 metaT.MNIout_coord_3];
        end
        % filter data
        goodChan = ~(isnan(mean(zccep_Tr,2, 'omitnan')));
        
        prefilterIdx =  goodChan &...
            metaT.MNIout_coord_1 .* metaT.MNIin_coord_1 > 0 &... % ipsilateral
            metaT.eudDist>5 & ...
            ~ismember(metaT.JP_label_in, {'', 'empty', 'NAN', 'NA'}) & ...
            ~ismember(metaT.JP_label_out, {'', 'empty', 'NAN', 'NA'}) & ...
            metaT.sCrossBorder == 0 & metaT.rCrossBorder == 0  ;
        % (...% exclude case when the only peak happens ealier than 10 ms
           % sum(table2array(metaT(:, ttpIdx)), 2, 'omitnan')  > artPeak_time ||...
           % sum(table2array(metaT(:, ttpIdx)), 2, 'omitnan')  ==0 ||...
           % isnan(sum(table2array(metaT(:, ttpIdx)), 2, 'omitnan'))...
           % )& ... 
        if strcmp(path, 'outflow')
            throughantTH = strcmp(metaT.JP_label_out , 'antTH') & prefilterIdx;
            throughmidTH = strcmp(metaT.JP_label_out , 'midTH') & prefilterIdx;
            throughpstTH = strcmp(metaT.JP_label_out , 'pstTH') & prefilterIdx;
        else
            throughantTH = strcmp(metaT.JP_label_in , 'antTH') & prefilterIdx;
            throughmidTH = strcmp(metaT.JP_label_in , 'midTH') & prefilterIdx;
            throughpstTH = strcmp(metaT.JP_label_in , 'pstTH') & prefilterIdx;
        end
        
        %------------------------------------------------------------
        throughROIs = {throughantTH, throughmidTH, throughpstTH};
        if strcmp(path, 'outflow')
            AllElec = cellstr(join([string(metaT.aSubID), string(metaT.record_chan)]));
        else
            AllElec = cellstr(join([string(metaT.aSubID), string(metaT.stim_chan)]));
        end
        %------
        timepoints = -10:600;
        weights_cell = cell(length(throughROIs),1);
        coordinates_cell = cell(length(throughROIs),1);
        jplabel_cell = cell(length(throughROIs),1);
        %------
        for r = 1:length(throughROIs)
            through = find(throughROIs{r});
            [eleclist, it] = unique(AllElec(through));
            idx = through(it);
            zccep_roi = zeros(length(eleclist),size(zccep_Tr,2));
            for i = 1:length(eleclist)
                keys = strsplit(eleclist{i}, ' ');
                asub = keys{1};
                el = keys{2};
                % find the pairs with same stim ROI and rec channels
                if strcmp(path, 'outflow')
                    commonROI = intersect(through, find(strcmp(metaT.record_chan, el) &...
                        strcmp(metaT.aSubID, asub)));
                else
                    commonROI = intersect(through, find(strcmp(metaT.stim_chan, el) &...
                        strcmp(metaT.aSubID, asub)));
                end
                zccep_abs = abs(zccep_Tr(commonROI,:));
                % threshold
                % threshold
                zccep_abs(zccep_abs<2.2) = 0;
                zccep_abs(:, art_time+onset) = 0;
                % use mean of the absolute value
                zccep_roi(i,:) = mean(zccep_abs, 1, 'omitnan');
            end
            
            weights_cell{r} = zccep_roi;
            
            if strcmp(path, 'outflow')
                jplabel_cell{r} = metaT.JP_label_in(idx);
            else
                jplabel_cell{r} = metaT.JP_label_out(idx);
            end
            coordinates_cell{r} = cords(idx,:);
        end
        
        %% loop through different clusters
        for c = 1:3
            
            cluster = clusters{c};
            switch cluster
                case 'early'
                    trange = [9 117]; % in ms
                case 'late1'
                    trange = [118 282];
                case 'late2'
                    trange = [284 592];
            end
            
            weights_cluster = cell(size(weights_cell));
            for r = 1:length(weights_cell)
                zccep_tmp = weights_cell{r};
                %zccep_tmp(zccep_tmp==0) = nan;
                % % take median across time points
                weights_cluster{r} = median(zccep_tmp(:, trange(1)+onset : trange(2)+onset), 2,'omitnan');
                % % take mean across time points
                % weights_cluster{r} = mean(zccep_tmp(:, trange(1)+onset : trange(2)+onset), 2,'omitnan');
                % take max across time points
                % weights_cluster{r} = max(zccep_tmp(:, trange(1)+onset : trange(2)+onset),[], 2,'omitnan');
            end
            
            if strcmp(path, 'outflow')
                anot = {'From antTH', 'From midTH', 'From pstTH'};
            else
                anot = {'To antTH', 'To midTH', 'To pstTH'};
            end
            % load brain surface data
            %             % brain surface data
            %             freesurferhome = getenv('FREESURFER_HOME');
            %             [cmcortex.vert, cmcortex.tri] = read_surf(fullfile(freesurferhome, 'subjects/fsaverage' , 'surf', 'rh.pial'));
                         pers_rule = readtable( fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/COHORT/ELECTRODE_VISUALIZATION_RULE_thal.csv'));
            %
            fn =  fullfile(brainDir, 'S900.R.midthickness_MSMAll.32k_fs_LR.surf.gii');
            this = gifti(fn);
            cortexFS0_2.vert = this.vertices;
            cortexFS0_2.tri = this.faces;
            
            %% plot surface weights
            clear cfg
            
            cfg.brainAlpha        = 0.75;
            cfg.suppress_colorbar = 1;
            cfg.gsp               = 70; % gaussian speading parameter
            cfg.figPosition       = [201   225   582   610];
            cfg.maxActRadius      = 350;
            cfg.anotFontSize      = 15;
            cfg.titleFontSize     = 14;
            cfg.titletext         = sprintf('%s %s peaks (%d-%d ms)', path, cluster, trange(1), trange(2));
            cfg.recElecDotSize    = 25;
            cfg.covgAdj           = 0.7;
            cfg.climsPerc         = repmat([0 0.96], 3,1);
            cfg.climMin           = 0;
            
            cfg.implicit_FaceAlpha = 0.05;
            cfg.implicit_EdgeAlpha = 0.1;
            
            switch path
                case 'outflow'
                    cfg.MarkerColor = repmat({[166,27,41]/256},3,1);
                    cfg.cmaps = cbrewer2('seq','RdPu',64, 'spline');
                    %litup_colors = [[166,27,41]/256; [255, 203, 1]/256;  [35, 118, 183]/256];% dark red[160, 28, 52],yellow, blue,
                case 'inflow'
                    cfg.MarkerColor = repmat({[207, 72, 19]/256},3,1);
                    cfg.cmaps = cbrewer2('seq','OrRd',64, 'spline');
                    %litup_colors = [[207, 72, 19]/256; [255, 203, 1]/256;  [35, 118, 183]/256];% dark red[160, 28, 52],yellow, blue,
            end
            % decide clims
            cfg.clims = repmat([0 8], 3, 1);
            
            %--------------------------------------------------------------------------------------------------
            figname = sprintf('BrainHeatMap CCEP %s cluster%d %d-%dms_median3.jpg', path, c, trange(1), trange(2));
            %output_folder = [];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            plotBrainSurfaceWeight(cortexFS0_2 , coordinates_cell, weights_cluster, 1, jplabel_cell, output_folder, anot, pers_rule, cfg,figname);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end % for cluster
    end % for path
    toc
    
end
