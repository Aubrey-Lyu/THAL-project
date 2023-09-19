% This script will generate plots for evoked potentials of the
% adjacent Rchans to the Schans, for visual inspection of stimulation artifacts.
% The adjacency will be determined by 10 mm
% radius in the Euclidean space. For now the MNI scanner RAS coordinates
% are used for this task. When available, the native scanner RAS
% coordinates should be used.

% 28 June, Dian Lyu, dl577@stanford.edu

%=======================================================
clear
dir_base = '/Volumes/G-DRIVE-Dian/DIAN';
home_dir = getuserdir;
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/CECS_Pipeline/CCEP_pipeline'))
addpath(fullfile(home_dir, 'Dropbox/scripts/external'))
mf = dir(fullfile(dir_base, 'neuralData/originalData', '**', 'global_CCEP_*_metadata.mat'));
metaFolders = {mf.folder}';
metaFiles = {mf.name}';

for fi = 8:length(mf)
    
    mfolder = metaFolders{fi};
    mfile = metaFiles{fi};
    if strcmp(mfile, 'global_CCEP_S22_182_DH_metadata.mat')
        continue;
    end
    load(fullfile(mfolder, mfile));
    sbj = ccep_metadata.subject{1};
    fprintf('Inspecting file: %s for %s ... \n', mfile, sbj )
    % project_name = 'CCEP';
    % [server_root, comp_root, code_root] = AddPaths('LBCN_iMac_Pro', project_name);
    % center = 'Stanford';
    % dirs = InitializeDirs(project_name, sbj, comp_root, server_root, code_root);
    if strcmp(sbj, 'S21_170_JL')
        continue;
    end

    if strcmp(sbj, 'S19_137_AF')
        subjVarF = fullfile(home_dir, 'Dropbox/subjVars', ['subjVar_' sbj '.mat']); % assuming the subjVar is saved in the same folder
    load(subjVarF)
        Edist = squareform(pdist(subjVar.elinfo.MNI_coord));
    else
        subjVarF = fullfile(mfolder, ['subjVar_' sbj '.mat']); % assuming the subjVar is saved in the same folder
    load(subjVarF)
    Edist = squareform(pdist(subjVar.elinfo.ScannerNativeRAS_coord));
    end
    %% plot
%     figure; set(gcf, 'Position', [239   412   972   863])
%     imagesc(Edist); colorbar; %caxis([1,20])
%     set(gca, 'XTick', subjVar.elinfo.chan_num(1:2:end), 'XTickLabel', subjVar.elinfo.FS_label(1:2:end),...
%         'YTick', subjVar.elinfo.chan_num(1:2:end), 'YTickLabel', subjVar.elinfo.FS_label(1:2:end))
%     title('Euclidean distance between contacts (mm)')
% 
%     output_dir = fullfile(strrep(mfolder, 'neuralData/originalData', 'Results/CCEP'), 'validation');
%     saveas(gcf,fullfile(output_dir, 'Euclidean distance between contacts.png'))
%     close all
    %% pass threshold = 10 mm to check artifacts
    % assign r = stim, c = record
    [r,c] = find(Edist < 10);
    CCEP_adjacent = struct();
    for i = 1:length(r)
        ri = r(i);
        ci = c(i);
        dist = round(Edist(ri,ci));
        sc = subjVar.elinfo.FS_label{ri};
        rc = subjVar.elinfo.FS_label{ci};
        block_name = unique(ccep_table.block_name(strcmp(ccep_table.sc1, sc)| strcmp(ccep_table.sc2, sc)));
        for ib = 1:length(block_name)
            block = block_name{ib};
            %----------------------------------
            CCEP_data_folder = strrep(mfolder, 'neuralData/originalData', 'computed_data/CCEP');
            load(fullfile(CCEP_data_folder, [block, '_CCEP.mat']));
            [~, rchans]=regexp(CCEP.recordingchannel','-','match','split');
            rchans = reshape([rchans{:}], 2, [])';
            idx = find(strcmp(rchans(:,1),  rc) | strcmp(rchans(:,2),  rc));
%             outdir = fullfile(output_dir, sprintf('ccep_inspect_%s_%s_%s-%s_%dmm.png', sbj, block, sc, rc, dist));
%             plot_CCEP_SpecChan(CCEP,sbj, idx, 1, outdir)
            CCEP_adjacent(ib).block_name = block;
            CCEP_adjacent(ib).stimChan   = 

        end
    end
end