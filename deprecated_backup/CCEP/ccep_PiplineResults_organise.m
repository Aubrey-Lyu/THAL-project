%% sort the metat table which records the CCEP newpipeline output/results and the anatomical information
% Before running it, to use the sortSubjVar_toLocal.m script to sort out the
% ScannerNativeRAS_coord and JP_labels.
%---------------------------------------------------------------------
% the output table will be used to concatenate the CCEP epochs, the order of the
% table should be corresponded to the CCEP matrix.
clear
addpath('~/Dropbox/scripts/external')
home_dir = '/home/dian';%getuserdir;
cmp_node = 'exo-lbcn'; % 'Dian_MacBook'; %'LBCN_iMac_Pro'; %
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

tic
%% sort path

project_name = 'CCEP';
center = 'Stanford';
freesurfer_home = '/Applications/freesurfer/7.2.0';
result_folder = fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/COHORT');
if exist(result_folder, 'dir') == 0; mkdir(result_folder); end
[server_root, comp_root, code_root] = AddPaths(cmp_node, project_name);

%parpool(12) % initialize number of cores

%% what analyses?
sort_metaT    = 1;
explore_metaT = 0;


%% Retrieve subject information
% Read the google sheets
subject_names = {'S22_178_AF', 'S21_172_KS','S21_166_TM', 'S20_152_HT', 'S21_167_MQ', 'S21_169_BH', 'S21_171_MM',...
    'S19_137_AF', 'S21_165_WN', 'S22_176_LB', 'S22_177_JM', 'S22_182_DH', 'S22_183_CR'...
    'S22_185_TW', 'S22_188_CB', 'S22_189_LMA', 'S22_190_AS', 'S22_191_KM', 'S22_192_LG', 'S22_193_AM','S23_194_PS',  'S23_195_MZ'};%, , 'S22_181_CB'}; %,

Cvec = brewermap(length(subject_names), 'Set3');
sbj_name_example = 'S23_195_MZ'; % short, no ITPC
dirs = InitializeDirs(project_name, sbj_name_example, comp_root, server_root, code_root);
loopin = 1:length(subject_names);

%% sort metaT
if sort_metaT == 1

    d = load(fullfile(dirs.data_root, 'originalData', sbj_name_example, ['global_CCEP_' sbj_name_example '_metadata.mat']));
    vars = d.ccep_table.Properties.VariableNames'; clear d
    subjVarPool = fullfile(home_dir, 'Dropbox/Stanford_Matters/subjVars');

    %% loop through subjects
    metaT0 = cell2table(cell(0,length(vars)), 'VariableNames', vars');

    for is = loopin

        sbj_name = subject_names{is};
        dirs = InitializeDirs(project_name, sbj_name, comp_root, server_root, code_root);

        d = load(fullfile(dirs.data_root, 'originalData', sbj_name, ['global_CCEP_' sbj_name '_metadata.mat']));

        vars_ = d.ccep_table.Properties.VariableNames';
        idx = find(ismember(vars, vars_));
        metaT0 = [metaT0; d.ccep_table(:, idx)];

    end

    %% add coordinates
    %subjVarDir = fullfile(home_dir, 'Dropbox/Stanford_Matters/subjVars');
    HotT = readtable(fullfile(home_dir, 'Dropbox/Stanford_Matters/data/cohort_data/SelfProj_allStims_fullCohort.csv'));
    HotT = unique(HotT, 'rows');
    metaT = [];

    for sn = loopin
        sbj_name = subject_names{sn};
        disp(['Registering for ' sbj_name '...'])
        %-----------------------------------------------
        subIdx = find(contains(metaT0.subject, sbj_name));
        subT = metaT0(subIdx,:);

        % sort subjVar
        subjVarDir = fullfile(dirs.server_root, ['SHICEP_' sbj_name], 'Analyzed Data');

        if exist(fullfile(subjVarPool, ['subjVar_', sbj_name, '_jplabel.mat']), 'file') == 2
            subjVarFile = fullfile(subjVarPool, ['subjVar_', sbj_name, '_jplabel.mat']);
        else

            if exist(fullfile(subjVarDir, ['subjVar_', sbj_name, '_jplabel.mat']), 'file') == 2
                subjVarFile = fullfile(subjVarDir, ['subjVar_', sbj_name, '_jplabel.mat']);
            else
                subjVarFile = fullfile(subjVarDir, ['subjVar_', sbj_name, '.mat']);
                warning('Sort out JP_label first!')
            end
        end

        clear subjVar
        load(subjVarFile);
        % clean up subjVar, do not include EXCLUDE channels
        subjVar.elinfo = subjVar.elinfo(~strcmp(subjVar.elinfo.JP_label, 'EXCLUDE'),:);

        % clean up subT, do not record anything that is out of the scope of subjVar
        subT = subT(ismember(subT.sc1, subjVar.elinfo.FS_label) &...
            ismember(subT.sc2, subjVar.elinfo.FS_label) & ...
            ismember(subT.rc1, subjVar.elinfo.FS_label) & ...
            ismember(subT.rc2, subjVar.elinfo.FS_label) ...
            ,:);

        %-----------------------------------------------
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
            index = find(strcmp(HotT.sbj_name, sbj_name) &...
                strcmp(HotT.FS_label, chanCell{i}) );
            [~, value_returned] = revalueCell( HotT.Hot_Cold_all, index ) ;
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

            id1 = idMat(i,1); id2 = idMat(i,2); id3 = idMat(i,3); id4 = idMat(i,4); %#ok<PFBNS>

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


        metaT = [metaT; subT]; clear subT
    end
    metaT.subject = strrep(metaT.subject, 'SHICEP_', '');
    writetable(metaT, fullfile(result_folder, 'table_CCEPnewpipOutput_wholebrain_anatomical_info.csv'))
    save( fullfile(result_folder, 'metaT.mat'), 'metaT')
end
toc
%% metaT group-level explore

if explore_metaT == 1

    metaT = readtable(fullfile(result_folder, 'table_CCEPnewpipOutput_wholebrain_anatomical_info.csv'));
    electrodes0 = [];
    Labels0 = {};
    for sn = loopin%length(subject_names)

        sbj_name = subject_names{sn};
        subIdx = find(contains(metaT.subject, sbj_name));
        subT = metaT(subIdx,:);
        [StimChans, StimIdx, ~] = unique(subT.stim_chan);
        [RecChans, RecIdx, ~]   = unique(subT.record_chan);

        % pass EXLUDE filter
        JP_exclude_in  = find(~strcmpi(subT.JP_label_in1, 'EXCLUDE'));
        JP_exclude_out = find(~strcmpi(subT.JP_label_out1, 'EXCLUDE'));

        StimIdx = intersect(StimIdx, JP_exclude_out);
        RecIdx  = intersect(RecIdx, JP_exclude_in);

        %% plot


        %                 p2 = scatter3(subT.MNIin_coord_1(RecIdx),subT.MNIin_coord_2(RecIdx), subT.MNIin_coord_3(RecIdx),  30, ...
        %                     'MarkerEdgeColor', Cvec(sn,:), 'MarkerFaceColor', Cvec(sn,:), 'DisplayName', ['recChan in ' strrep(sbj_name, '_', '\_')]);
        %                 p2.MarkerFaceAlpha = .5;
        %  legend(p2, strrep(sbj_name, '_', '\_'))
        electrodes0 = [electrodes0; subT.MNIin_coord_1(RecIdx) subT.MNIin_coord_2(RecIdx) subT.MNIin_coord_3(RecIdx)];

        chanOut = subT.stim_chan(StimIdx);
        chanIn  = subT.record_chan(RecIdx);

        strs = strsplit(sbj_name,'_');
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
                 hemi_Idx = find(electrodes0(:,1) <= 0);
            elseif h == 2
                 hemi_Idx = find(electrodes0(:,1) > 0);
            end

            electrodes = electrodes0(hemi_Idx,:);
            Labels = Labels0(hemi_Idx);

            hold on;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [~, elec_surf] = ctmr_gauss_plot2(cmcortex,electrodes, ones(size(electrodes,1),1), hem, per);
            alpha(0.75)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elec = elec_surf; %electrodes; 
            p1 = scatter3(elec(:,1),elec(:,2), elec(:,3),  35, ...
                'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerFaceColor', [0.6 0.6 0.6]);
            p1.MarkerFaceAlpha = 0.8;

            % add labels
           % text(elec(:,1),elec(:,2), elec(:,3), Labels, 'FontSize',6);

            figname = sprintf('Total_RecElecCoverage_Recording_%s_%s_surf.jpg', per, hem);
            exportgraphics(gcf, fullfile(figpath, figname),'Resolution', 300)

        end
        %if ps == 2 && h == 2
        %             legend('show');
        %end
    end
    hold off
    %     %% check network affiliation
    %
    %     % 1. stimulated from PMC
    %     pmcT  = T(strcmp(T.JP_label_out_1, 'PMC') | strcmp(T.JP_label_out_2, 'PMC'),:);% | ...
    %             %contains(lower(T.Destr_long_out), 'PCC'),:);
    %     thalT = T(contains(lower(T.JP_label_out_1), 'thalam') | contains(lower(T.JP_label_out_2), 'thalam'),:);%| ...
    %              % contains(lower(T.Destr_long_out), 'thalam'), :);
    %
    %

end

toc