%% sort the metat table which records the CCEP newpipeline output/results and the anatomical information
% Before running it, to use the sortSubjVar_toLocal.m script to sort out the
% ScannerNativeRAS_coord and JP_labels.
%---------------------------------------------------------------------
% the output table will be used to concatenate the CCEP epochs, the order of the
% table should be corresponded to the CCEP matrix.
% dl577@stanford.edu, Apr 27 2023
clear

addpath('~/Dropbox/scripts/external')
home_dir = '/home/dian'; % getuserdir;
cmp_node = 'exo-lbcn'; % 'Dian_MacBook'; %'LBCN_iMac_Pro'; %'lvdian-Precision-3660'; 
figpath  = fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/PLOTS');

addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/SELF-project/CCEP'))
addpath(fullfile(home_dir, 'Dropbox/scripts/external/spm12'));
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/lbcn_preproc/data_convertion'));
addpath(genpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/lbcn_preproc/cleaning')));
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/lbcn_preproc/vizualization'));
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/lbcn_preproc/data_processing'));
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/lbcn_preproc/lbcn_personal/personal'));
addpath(genpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY')));
addpath(fullfile(home_dir, 'Dropbox/scripts/external/ColorBrewer'))
addpath(fullfile(home_dir, 'Dropbox/scripts/my_functions'))


%% sort path

project_name = 'CCEP';
center = 'Stanford';

data_result_folder = ...
    fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore5_locked');
if exist(data_result_folder, 'dir') == 0; mkdir(data_result_folder); end

cohort_result_folder = data_result_folder; %fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/COHORT');
if exist(cohort_result_folder, 'dir') == 0; mkdir(cohort_result_folder); end

[server_root, comp_root, code_root] = AddPaths(cmp_node, project_name);

%parpool(12) % initialize number of cores

%% what analyses?
sort_metaT      = 1;
explore_metaT   = 1;
flatten_CCEP    = 1;
redo_activation = 1;

%% Retrieve subject information

subject_names = {'S21_166_TM','S21_167_MQ', 'S21_169_BH', 'S21_170_JL','S21_171_MM',...
    'S21_172_KS', 'S22_176_LB','S22_177_JM','S22_178_AF', 'S22_181_CB', 'S22_182_DH',...
    'S22_183_CR','S22_185_TW', 'S22_188_CB',...
    'S22_189_LMA','S22_190_AS','S22_192_LG','S22_193_AM', 'S23_194_PS','S23_195_MZ', 'S23_196_HL'...
    'S23_197_TA', 'S23_198_JP', ...% 'S22_186_ML', 'S22_175_CT' (noCCEP)
    'S23_199_GB', 'S23_201_JG', 'S23_202_KC'}; % new 3 subj added

subject_names = sort(subject_names)';

Cvec = brewermap(length(subject_names), 'Set3');
sbj_name_example = 'S22_178_AF'; % short, no ITPC
dirs = InitializeDirs(project_name, sbj_name_example, comp_root, server_root, code_root);
loopin = 1:length(subject_names);

%% sort metaT
if sort_metaT == 1
tic
    d = load(fullfile(dirs.data_root, 'originalData', sbj_name_example, ['global_CCEP_' sbj_name_example '_metadata.mat']));
    vars = d.ccep_table.Properties.VariableNames'; clear d
    subjVarPool = fullfile(home_dir, 'Dropbox/Stanford_Matters/subjVars');

    %% loop through subjects
    metaT0 = cell2table(cell(0,length(vars)), 'VariableNames', vars');

    for is = loopin

        sbj_id = subject_names{is};
        dirs = InitializeDirs(project_name, sbj_id, comp_root, server_root, code_root);

        d = load(fullfile(dirs.data_root, 'originalData', sbj_id, ['global_CCEP_' sbj_id '_metadata.mat']));

        vars_ = d.ccep_table.Properties.VariableNames';
        [~,~,idx] = intersect(vars, vars_, 'stable');
        metaT0 = [metaT0; d.ccep_table(:, idx)];

    end

    %% add coordinates
    subjVarDir = fullfile(home_dir, 'Dropbox/Stanford_Matters/subjVars');
    HotT = readtable(fullfile(home_dir, 'Dropbox/Stanford_Matters/data/cohort_data/SelfProj_allStims_fullCohort.csv'));
    HotT = unique(HotT, 'rows');
    metaT = [];

    for sn = loopin
        sbj_id = subject_names{sn};
        disp(['Registering ' sbj_id '...'])
        %-----------------------------------------------
        subIdx = find(contains(metaT0.subject, sbj_id));
        subT = metaT0(subIdx,:);

        % sort subjVar
        %subjVarDir = fullfile(dirs.server_root, ['SHICEP_' sbj_name], 'Analyzed Data');

        if exist(fullfile(subjVarPool, ['subjVar_', sbj_id, '_jplabel.mat']), 'file') == 2
            subjVarFile = fullfile(subjVarPool, ['subjVar_', sbj_id, '_jplabel.mat']);
        else

            if exist(fullfile(subjVarDir, ['subjVar_', sbj_id, '_jplabel.mat']), 'file') == 2
                subjVarFile = fullfile(subjVarDir, ['subjVar_', sbj_id, '_jplabel.mat']);
            else
                subjVarFile = fullfile(subjVarDir, ['subjVar_', sbj_id, '.mat']);
                warning('Sort out JP_label first!')
            end
        end

        clear subjVar
        load(subjVarFile);
        % clean up subjVar, do not include EXCLUDE channels
        subjVar.elinfo = subjVar.elinfo(~strcmpi(subjVar.elinfo.JP_label, 'empty'),:);
        subjVar.elinfo = subjVar.elinfo(~isnan(subjVar.elinfo.MNI_coord(:,1)),:);
        % clean up subT, do not record anything that is out of the scope of subjVar
        subT = subT(ismember(lower(subT.sc1), lower(subjVar.elinfo.FS_label)) &...
            ismember(lower(subT.sc2), lower(subjVar.elinfo.FS_label)) & ...
            ismember(lower(subT.rc1), lower(subjVar.elinfo.FS_label)) & ...
            ismember(lower(subT.rc2), lower(subjVar.elinfo.FS_label)) ...
            ,:);

        %-----------------------------------------------
        subT.aSubID         = cell(size(subT,1),0); % anonymous subID
        subT.record_hot     = nan(size(subT,1), 1);
        subT.stim_hot       = nan(size(subT,1), 1);

        subT.JP_label_out1   = cell(size(subT,1),0);
        subT.JP_label_out2   = cell(size(subT,1),0);
        subT.JP_label_in1    = cell(size(subT,1),0);
        subT.JP_label_in2    = cell(size(subT,1),0);

        subT.Yeo7_out1       = cell(size(subT,1),0);
        subT.Yeo7_out2       = cell(size(subT,1),0);
        subT.Yeo7_in1        = cell(size(subT,1),0);
        subT.Yeo7_in2        = cell(size(subT,1),0);

        subT.MNIout_coord      = nan(size(subT,1), 3);
        subT.MNIin_coord       = nan(size(subT,1), 3);
        subT.LEPTOout_coord    = nan(size(subT,1), 3);
        subT.LEPTOin_coord     = nan(size(subT,1), 3);
        subT.NATIVEout_coord   = nan(size(subT,1), 3);
        subT.NATIVEin_coord    = nan(size(subT,1), 3);

        subT.eudDist         = nan(size(subT,1), 1);
        %--------------------------------------------------------------------
        chanCell = reshape([subT.sc1 subT.sc2 subT.rc1 subT.rc2],[],1);
        jplabelCell = cell(size(chanCell));
        Yeo7Cell = cell(size(chanCell));
        idMat = zeros(size(chanCell));
        listhot = nan(size(chanCell));
        recordhot = nan(size(chanCell));

        clear i
        parfor i = 1:size(idMat,1)
            idMat(i) = find(strcmpi(subjVar.elinfo.FS_label, chanCell{i}));
            jplabelCell{i,1} = subjVar.elinfo.JP_label{idMat(i)};
            Yeo7Cell{i,1} = subjVar.elinfo.Yeo7{idMat(i)};
            % SORT EBS HOT: involving using another table
            index = find(strcmpi(HotT.sbj_name, sbj_id) &...
                strcmpi(HotT.FS_label, chanCell{i}) );
            [~, value_returned] = revalueArray( HotT.Hot_Cold_all, index ) ;
            listhot(i) = value_returned;
        end

        jplabelCell = reshape(jplabelCell, [], 4);
        subT.JP_label_out1 = jplabelCell(:,1);
        subT.JP_label_out2 = jplabelCell(:,2);
        subT.JP_label_in1  = jplabelCell(:,3);
        subT.JP_label_in2  = jplabelCell(:,4);

        Yeo7Cell = reshape(Yeo7Cell, [], 4);
        subT.Yeo7_out1    = Yeo7Cell(:,1);
        subT.Yeo7_out2    = Yeo7Cell(:,2);
        subT.Yeo7_in1     = Yeo7Cell(:,3);
        subT.Yeo7_in2     = Yeo7Cell(:,4);

        listhot = reshape(listhot, [], 4);
        nobothnan = ~isnan(listhot(:,1)) & ~isnan(listhot(:,2));
        subT.stim_hot(nobothnan) = sum(listhot(nobothnan,[1 2]), 2,'omitnan');
        nobothnan = ~isnan(listhot(:,3)) & ~isnan(listhot(:,4));
        subT.record_hot(nobothnan) = sum(listhot(nobothnan,[3 4]), 2,'omitnan');
        %--------------
        idMat = reshape(idMat,[],4);

        MNIout_coord    = nan(size(idMat,1), 3);
        NATIVEout_coord = nan(size(idMat,1), 3);
        LEPTOout_coord  = nan(size(idMat,1), 3);

        MNIin_coord    = nan(size(idMat,1), 3);
        NATIVEin_coord = nan(size(idMat,1), 3);
        LEPTOin_coord  = nan(size(idMat,1), 3);

        eudDist = nan(size(idMat,1), 1); %between stim and record centroids

        for i = 1:size(idMat,1)

            id1 = idMat(i,1); id2 = idMat(i,2); id3 = idMat(i,3); id4 = idMat(i,4);

            coord1 = subjVar.elinfo.MNI_coord(id1,:);
            coord2 = subjVar.elinfo.MNI_coord(id2,:);
            MNIout_coord(i,:) = (coord1+coord2)./2;

            coord1 = subjVar.elinfo.LEPTO_coord(id1,:);
            coord2 = subjVar.elinfo.LEPTO_coord(id2,:);
            LEPTOout_coord(i,:) = (coord1+coord2)./2;

            coord1 = subjVar.elinfo.ScannerNativeRAS_coord(id1,:);
            coord2 = subjVar.elinfo.ScannerNativeRAS_coord(id2,:);
            NATIVEout_coord(i,:) = (coord1+coord2)./2;

            coord3 = subjVar.elinfo.MNI_coord(id3,:);
            coord4 = subjVar.elinfo.MNI_coord(id4,:);
            MNIin_coord(i,:) = (coord3+coord4)./2;

            coord3 = subjVar.elinfo.LEPTO_coord(id3,:);
            coord4 = subjVar.elinfo.LEPTO_coord(id4,:);
            LEPTOin_coord(i,:) = (coord3+coord4)./2;

            coord3 = subjVar.elinfo.ScannerNativeRAS_coord(id3,:);
            coord4 = subjVar.elinfo.ScannerNativeRAS_coord(id4,:);
            NATIVEin_coord(i,:) = (coord3+coord4)./2;

            % euclidean distance between stim and record centroids
            eudDist(i) = sqrt(sum((NATIVEout_coord(i,:) - NATIVEin_coord(i,:)).^2));
        end

        subT.MNIout_coord    = MNIout_coord; clear MNIout_coord
        subT.MNIin_coord     = MNIin_coord;  clear MNIin_coord
        subT.LEPTOout_coord  = LEPTOout_coord; clear LEPTOout_coord
        subT.LEPTOin_coord   = LEPTOin_coord; clear LEPTOin_coord
        subT.NATIVEout_coord = NATIVEout_coord; clear NATIVEout_coord
        subT.NATIVEin_coord  = NATIVEin_coord; clear NATIVEin_coord
        subT.eudDist         = eudDist; clear eudDist
        % create an anonymous ID
        subKeys = strsplit(sbj_id, '_');
        subT.aSubID = cellstr(repmat(sprintf('S%02d_%s', sn, subKeys{2}), size(subT,1),1));
        
        metaT = [metaT; subT]; clear subT
    end

    %% filter metaT:
    % 1. get rid of stimulation/recording sites with both EXCLUDE label
    % 2. get rid of stim_chan == record_chan
    excludeIdx = ((strcmpi(metaT.JP_label_in1, 'EXCLUDE') & strcmpi(metaT.JP_label_in2, 'EXCLUDE')) | ...
        (strcmpi(metaT.JP_label_out1, 'EXCLUDE') & strcmpi(metaT.JP_label_out2, 'EXCLUDE'))) | ...
        cellfun(@strcmpi, metaT.stim_chan, metaT.record_chan);
    metaT = metaT(~excludeIdx,:);
   
    %% sort metaT anatomical label
    % 1. load the csv with JP labels cleaned up,when there is an exclusion in
    % the pair, use the non-exluded ROI as the ROI label
    JP_label_in = metaT.JP_label_in1;
    JP_label_in(strcmpi(metaT.JP_label_in1, 'EXCLUDE')) = ...
        metaT.JP_label_in2(strcmpi(metaT.JP_label_in1, 'EXCLUDE'));
    JP_label_out = metaT.JP_label_out1;
    JP_label_out(strcmpi(metaT.JP_label_out1, 'EXCLUDE')) = ...
        metaT.JP_label_out2(strcmpi(metaT.JP_label_out1, 'EXCLUDE'));
    % 2. decide what is cross border pair, and write as a separate column
    % stim
    noExclusion = ~(strcmpi(metaT.JP_label_out1, 'EXCLUDE') | strcmpi(metaT.JP_label_out2, 'EXCLUDE'));
    sCrossBorder = zeros(size(metaT,1),1);
    sCrossBorder(noExclusion & ~cellfun(@isequal, metaT.JP_label_out1, metaT.JP_label_out2)) = 1;
    % record
    noExclusion = ~(strcmpi(metaT.JP_label_in1, 'EXCLUDE') | strcmpi(metaT.JP_label_in2, 'EXCLUDE'));
    rCrossBorder = zeros(size(metaT,1),1);
    rCrossBorder(noExclusion & ~cellfun(@isequal, metaT.JP_label_in1, metaT.JP_label_in2)) = 1;
    % 3. sort naming
    metaT.JP_label_out = ListSortAnatLabel_THAL(JP_label_out, 1);
    metaT.JP_label_in = ListSortAnatLabel_THAL(JP_label_in, 1);
    % 4. add crossBorder column to the metaT
    metaT.sCrossBorder = sCrossBorder;
    metaT.rCrossBorder = rCrossBorder;
    % reorder the subtable so that aSubID is in the front (1st = subject original ID)
     vars = metaT.Properties.VariableNames;
     idIdx = find(strcmpi(vars, 'aSubID'));
     allIdx = 1:length(vars);
     allIdx(idIdx) = [];
     reorderIdx = [allIdx(1),idIdx,allIdx(2:end)];
     metaT = metaT(:, reorderIdx);

    %% save metaT
    writetable(metaT, fullfile(cohort_result_folder, 'table_CCEPnewpipOutput_wholebrain_anatomical_info.csv'))
    save(fullfile(cohort_result_folder, 'metaT.mat'), 'metaT')
    clear metaT0
    toc
end


%% metaT group-level explore

if explore_metaT == 1
  
    metaT = readtable(fullfile(cohort_result_folder, 'table_CCEPnewpipOutput_wholebrain_anatomical_info.csv'));
    subject_IDs  = unique(metaT.aSubID);
    electrodes0 = [];  electrodes1 = [];
    Labels0 = {};
    % for electrodes frequency table
   TBL_stim = table(); TBL_rec = table();

    for sn = 1:length(subject_IDs)

        sbj_id = subject_IDs{sn};
        subIdx = find(contains(metaT.aSubID, sbj_id));
        subT = metaT(subIdx,:);
        [StimChans, StimIdx, ~] = unique(subT.stim_chan);
        [RecChans, RecIdx, ~]   = unique(subT.record_chan);
        %
        % pass cross border filter
        noCrossBorder_in  = find(subT.rCrossBorder == 0);
        noCrossBorder_out = find(subT.sCrossBorder == 0);

        StimIdx = intersect(StimIdx, noCrossBorder_out);
        RecIdx  = intersect(RecIdx, noCrossBorder_in);

        %% sort electrodes and frequency table

        electrodes0 = [electrodes0; subT.MNIin_coord_1(RecIdx) subT.MNIin_coord_2(RecIdx) subT.MNIin_coord_3(RecIdx)];
        electrodes1 = [electrodes1; subT.MNIin_coord_1(StimIdx) subT.MNIin_coord_2(StimIdx) subT.MNIin_coord_3(StimIdx)];
%
        JPlabel_stim = subT.JP_label_out(intersect(StimIdx, find(~strcmpi(subT.JP_label_out, 'EXCLUDE'))));
        JPlabel_record = subT.JP_label_in(intersect(RecIdx, find(~strcmpi(subT.JP_label_in, 'EXCLUDE'))));

        tbl = tabulate(JPlabel_stim);
        stbl             = table();
        stbl.subID       = cellstr(repmat(sbj_id, size(tbl,1) ,1));
        stbl.brainArea   = tbl(:,1);
        stbl.elecCount   = cell2num(tbl(:,2));
        stbl.elecPercSbj = cell2num(tbl(:,3));

        rtbl_rec  = tabulate(JPlabel_stim);
        rtbl             = table();
        rtbl.subID       = cellstr(repmat(sbj_id, size(tbl,1) ,1));
        rtbl.brainArea   = tbl(:,1);
        rtbl.elecCount   = cell2num(tbl(:,2));
        rtbl.elecPercSbj = cell2num(tbl(:,3));
        

        if sn == 1
            TBL_stim  = stbl;
            TBL_rec   = rtbl;
        else
            TBL_stim = [TBL_stim; stbl];
            TBL_rec  = [TBL_rec; rtbl];
        end

        % save
        writetable(TBL_stim, fullfile(cohort_result_folder, 'CCEPstimPair_noCrossRegion_JPlabel_frequency_table.csv'))
        writetable(TBL_rec, fullfile(cohort_result_folder, 'CCEPrecPair_noCrossRegion_JPlabel_frequency_table.csv'))
%}
%% sort labels
        chanOut = subT.stim_chan(StimIdx);
        chanIn  = subT.record_chan(RecIdx);

        strs = strsplit(sbj_id,'_');
        labels = cellstr([repmat(strs{2}, length(chanIn),1), char(chanIn)]);
        Labels0 = [Labels0; labels];

    end

    %% plot coverage
    close all
    pers = {'lateral', 'medial'};
    hemi = {'left', 'right'};

    for h = 1:2

        hem = hemi{h};
        [cmcortex.vert, cmcortex.tri]=read_surf(fullfile(dirs.fsDir_local , 'surf',[hem(1) 'h.pial']));

        for ps = 1:2
            per = pers{ps};
            figure('Position', [562   667   711   540])

            % pass filter for selecting the hemisphere
            if h == 1
                % record electrodes
                electrodes_0 = electrodes0(electrodes0(:,1) <= 0,:);
                % stim electrodes
                electrodes_1 = electrodes1(electrodes1(:,1) <= 0,:);
            elseif h == 2
                % record electrodes
                electrodes_0 = electrodes0(electrodes0(:,1) > 0,:);
                % stim electrodes
                electrodes_1 = electrodes1(electrodes1(:,1) > 0,:);
            end

            %Labels = Labels0(hemi_Idx);
            hold on;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [~, elec_surf_0] = ctmr_gauss_plot2(cmcortex, electrodes_0, ones(size(electrodes_0,1),1), hem, per);
            elec_surf_1 = elec2surf(cmcortex, electrodes_1); % stim Chan

            alpha(0.75)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % plot recording channels
            elec = elec_surf_0; %electrodes;
            p1 = scatter3(elec(:,1),elec(:,2), elec(:,3),  35, ...
                'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', [0.6 0.6 0.6]);
            p1.MarkerFaceAlpha = 0.8;
%{
            % plot stimulation channels
            elec = elec_surf_1; %electrodes;
            p2 = scatter3(elec(:,1)-0.001,elec(:,2), elec(:,3),  40, 'x', 'LineWidth',2);
            p2.MarkerEdgeColor = [0.95 0.4 0.1];
            p2.MarkerFaceAlpha = 1;

            % add labels
            % text(elec(:,1),elec(:,2), elec(:,3), Labels, 'FontSize',6);
%}
            figname = sprintf('Total_RecElecCoverage_Recording_%s_%s_surf3.jpg', per, hem);
            exportgraphics(gcf, fullfile(figpath, figname),'Resolution', 300)

        end
    end
    hold off
    %     %% check network affiliation
    %
    %     % 1. stimulated from PMC
    %     pmcT  = T(strcmpi(T.JP_label_out_1, 'PMC') | strcmpi(T.JP_label_out_2, 'PMC'),:);% | ...
    %             %contains(lower(T.Destr_long_out), 'PCC'),:);
    %     thalT = T(contains(lower(T.JP_label_out_1), 'thalam') | contains(lower(T.JP_label_out_2), 'thalam'),:);%| ...
    %              % contains(lower(T.Destr_long_out), 'thalam'), :);
    %
    %
    toc
end


%% organize data according to metaT

if flatten_CCEP == 1
    tic
    if ~exist('metaT', 'var')
        metaT = readtable(fullfile(cohort_result_folder, 'table_CCEPnewpipOutput_wholebrain_anatomical_info.csv'));
    end
    vnames = metaT.Properties.VariableNames;
    block_names = unique(metaT.block_name);
    CCEP_flat =  nan(size(metaT,1), 2201); % each row corresponding to the metaT; using trial average
    CCEP_flat2 =  nan(size(metaT,1), 2201); % each row corresponding to the metaT; using trial median
    block_before = 'dummy';

    parfor ir = 1:size(metaT,1)

        block = metaT.block_name{ir};

        %if ~strcmp(block_before, block)

        sbj_id = char(unique(metaT.subject(strcmp(metaT.block_name, block))));
        ccep_dir = fullfile(comp_root, 'computed_data/CCEP', sbj_id);
        m = load(fullfile(comp_root, 'neuralData/originalData', sbj_id, ['global_CCEP_' sbj_id '_metadata.mat']));

        fname = fullfile(ccep_dir, [block '_CCEP.mat']);

        d = load(fname); % CCEP

        % end

        block_before = block;

        rec_idx = find(strcmpi(d.CCEP.recordingchannel, metaT.record_chan{ir}));

        wave = squeeze(d.CCEP.wave(:,rec_idx,:))';
        %wave_clean = nan(size(wave));

        %% exclude bad trials
        reject_trials = table2array(metaT(ir,contains(vnames, 'reject_trials')))';
        reject_trials(isnan(reject_trials)) = [];

        wv_meanTr  = mean(abs(wave),1, 'omitnan');
        thr_data   = mean(wv_meanTr, 'omitnan') + 4.*std(wv_meanTr, 'omitnan');
        reject_trials  = [reject_trials; find(wv_meanTr >= thr_data)'];
        % [rb, cb]       = find(wave>100); % containing any extreme values bigger than 100 z-scores
        % reject_trials  = [reject_trials; unique(cb)];
        reject_trials  = unique(reject_trials);

        % detected from the pipeline
        trial_No = 1:size(d.CCEP.wave, 1);
        trial_No(ismember(trial_No, reject_trials)) = [];

        wave_clean = wave(:, trial_No);
        wave_avg   = mean(wave_clean, 2, 'omitnan');
        wave_med   = median(wave_clean, 2, 'omitnan');

        % put in the big CCEP mat
        CCEP_flat(ir,:) = wave_avg';
        CCEP_flat2(ir,:) = wave_med';

    end

    %% save flat CCEP
    save(fullfile(data_result_folder, 'CCEP_all_flat_meanTr.mat'), 'CCEP_flat', '-v7.3')
    save(fullfile(data_result_folder, 'CCEP_all_flat_medianTr.mat'), 'CCEP_flat2', '-v7.3')

    %% clean CCEP
    [badChan, zccep_clean]   = detectCCEP_badChan(CCEP_flat, 1);
    save(fullfile(data_result_folder, 'CCEP_all_flat_meanTr_cleaned.mat'), 'badChan', 'zccep_clean', '-v7.3');
    
    [badChan, zccep_clean2]   = detectCCEP_badChan(CCEP_flat2, 1);
    save(fullfile(data_result_folder, 'CCEP_all_flat_medianTr_cleaned.mat'), 'badChan', 'zccep_clean2', '-v7.3');

    toc
end

%%
if redo_activation == 1
    % use multiple detection rules
[activated_redo, act_loc] = redetect_activation(metaT, CCEP_flat); %may take a few minutes

activated = metaT.activated;
CECs_activation = metaT.CECS_activation;

%% activation stratege: democratic voting
cd(cohort_result_folder)
diary activationRedo.log

act_vote = 0+(sum([activated_redo, activated, CECs_activation], 2) >= 2) ;
% default-act
evaluation = makeConfusionMat(activated, act_vote);
fprintf(...
    '\ndefault-act has error rate = %.2f, accuracy = %.2f, \nsensitivity = %.2f, specificity = %.2f, \nprecision = %.2f, false positive rate = %.2f, \ncorrelation with the observation is %.2f\n', ...
    evaluation.ERR, evaluation.ACC, evaluation.SN, evaluation.SP, ...
    evaluation.PREC, evaluation.FPR, evaluation.MCC)
% simpRule-act
evaluation = makeConfusionMat(activated_redo, act_vote);
fprintf(...
    '\nsimpRule-act has error rate = %.2f, accuracy = %.2f, \nsensitivity = %.2f, specificity = %.2f, \nprecision = %.2f, false positive rate = %.2f, \ncorrelation with the observation is %.2f\n', ...
    evaluation.ERR, evaluation.ACC, evaluation.SN, evaluation.SP, ...
    evaluation.PREC, evaluation.FPR, evaluation.MCC)
% CEC-act
evaluation = makeConfusionMat( CECs_activation, act_vote);
fprintf(...
    '\nCEC-act has error rate = %.2f, accuracy = %.2f, \nsensitivity = %.2f, specificity = %.2f, \nprecision = %.2f, false positive rate = %.2f, \ncorrelation with the observation is %.2f\n', ...
    evaluation.ERR, evaluation.ACC, evaluation.SN, evaluation.SP, ...
    evaluation.PREC, evaluation.FPR, evaluation.MCC)

diary off

%% add the act_vote to meta table, and reorganize

varnames = metaT.Properties.VariableNames;
pkIdx = contains(varnames, 'pks_time_');
pktimes = table2array(metaT(:,pkIdx));
pks_time = cell(length(act_loc),1);

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
    else
        pks_time{i} = nan;
    end
end

metaT.activated_default = activated;
metaT.activated_SimRule = activated_redo;
metaT.activated = act_vote;
metaT(:, pkIdx) = [];
metaT.pks_time = pks_time;
writetable(metaT, fullfile(cohort_result_folder, 'table_CCEPnewpipOutput_wholebrain_anatomical_info_activationRedone.csv'));
disp('Activation redone is written into metaT, saved in the hard disk.')
end