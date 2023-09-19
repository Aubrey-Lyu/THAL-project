% this script will reduce the dimension of the spectrum for UMAP learn
% part1: average trials to see if there is a general differences between
% activated vs non-activated, cor vs thal categorties
% version 2 will use UMAP semi-supervised learned activation scores.
% reduction scheme is 3.
% The key variable used in this script (UMAP identified activation) is the
% output of the prepare_spectrumData_for_umap_part2.m file (using the
% reduce 1 scheme) and the python scripts (UMAP_actLeanr_sbj*) in the folder
% "/data/dian/Dropbox/Stanford_Matters/data/THAL/UMAP"

home_dir = '/data/dian';
addpath(genpath(fullfile(home_dir, 'Dropbox/scripts/external/cbrewer2/cbrewer2')));
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/personal_validation'))
addpath('~/scripts/my_functions')
addpath(genpath('~/scripts/Stanford/ThalamocoricalLoop-project'))
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/CECS_Pipeline/CCEP_pipeline/HFB plot/subfunctions'));

data_folder = '~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore5_locked';
umap_outdir = '~/Dropbox/Stanford_Matters/data/THAL/UMAP/WITHINSBJ_semisupervise';
result_folder = '~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore5_locked';

data_folder_intermediate = fullfile(data_folder,'render2');
if ~exist(data_folder_intermediate,'dir'); mkdir(data_folder_intermediate);end

plot_folder = '~/Dropbox/Stanford_Matters/data/THAL/CCEP/Plots/explore5/render2-original';
if ~exist(plot_folder, 'dir'); mkdir(plot_folder);end
% load metaTable
metaT = readtable('~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore5_locked/table_CCEPnewpipOutput_wholebrain_anatomical_info_activationRedone.csv');
%
% load data
%load('~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore5_locked/CCEP_all_flat_meanTr_cleaned.mat');%, 'CCEP_flat' each row corresponding to the metaT
%

afns = dir(fullfile(umap_outdir, 'act*.txt'));
sblist = cell(length(afns),1);
for i = 1:length(afns)
    keys = strsplit(afns(i).name, '_');
    sb = strjoin([keys(4:5), {keys{6}(1:end-4)}], '_');
    sblist{i} = sb;
end

% sblist = unique(metaT.subject);
% vars = metaT.Properties.VariableNames;
% ttpIdx = find(contains(vars, 'pks_time'));
% peakTimeMat = table2array(metaT(:,ttpIdx));

load(fullfile(data_folder, 'CCEP_all_flat_meanTr_cleaned.mat'), 'badChan');

goodChan = ones(size(metaT,1),1);
goodChan(badChan) = 0;
cd(fullfile(data_folder, 'spectCCEP'));

%% separating CCEP within and between hemispheres
hemi_CAT = {'ipsi', 'contr'};
connTypes = {...
    'THAL', 'COR';...
    'COR', 'THAL';...
    'THAL', 'THAL';...
    'COR', 'COR';...
    };
%-----------------------------------------------
for ih = 1:length(hemi_CAT)

    hcat = hemi_CAT{ih};
    if ih == 1
        hemiIdx = (metaT.MNIout_coord_1 .* metaT.MNIin_coord_1 > 0);  % ipsilateral
    else
        hemiIdx = (metaT.MNIout_coord_1 .* metaT.MNIin_coord_1 <= 0);  % contralateral
    end
    % RECORD FOR GROUP LEVEL
    PW = nan(59, 501, 2, length(connTypes), length(sblist));
    PH = nan(59, 501, 2, length(connTypes), length(sblist));
    PC = nan(59, 501, 2, length(connTypes), length(sblist));
    metaT_filterIdx = [];
    UMAP_Activation = [];

    for is = 1:length(sblist)

        sb = sblist{is};

        power = load(sprintf('spectCCEP_power_%s.mat', sb));
        phase = load(sprintf('spectCCEP_phase_%s.mat', sb));
        itpc = load(sprintf('spectCCEP_itpc_%s.mat', sb));

        prefilterIdx =  find(goodChan & ...
            hemiIdx &... %metaT.MNIout_coord_1 .* metaT.MNIin_coord_1 > 0 &... % ipsilateral
            metaT.eudDist>5 & ...
            ~ismember(metaT.JP_label_in, {'', 'empty', 'NAN', 'NA'}) & ...
            ~ismember(metaT.JP_label_out, {'', 'empty', 'NAN', 'NA'}) & ...
            metaT.sCrossBorder == 0 & metaT.rCrossBorder == 0)  ;

        prefilterIdx_act =  find(goodChan & ...
            ~ismember(metaT.JP_label_in, {'', 'empty', 'NAN', 'NA'}) & ...
            ~ismember(metaT.JP_label_out, {'', 'empty', 'NAN', 'NA'}))  ; % refer to the script 'prepare_spectrumData_for_umap_part2.m'

        [filterIdx, preIdxsub, ~] = intersect(phase.idx_in_metaT, prefilterIdx);
        [filterIdx_act, ~, ~]     = intersect(phase.idx_in_metaT, prefilterIdx_act);

        subT = metaT(filterIdx,:);
        powerDat = power.CCEP_power(preIdxsub,:,:);
        phaseDat = phase.CCEP_phase(preIdxsub,:,:);
        itpcDat  = itpc.CCEP_itpc(preIdxsub,:,:);

        [~, ~, reIdx_actLab] = intersect(filterIdx, filterIdx_act,  'stable');

        % load umap activation labels
        afn = fullfile(umap_outdir,sprintf('actLearn_2clusters_anot_%s.txt', sb));

        actlab0 = importdata(afn);
        actlab0 = actlab0-1; % in this labelling, 2=activated, 1=inactivated, -1=exclude
        actlab = actlab0(reIdx_actLab); % actIdx0 follows the filterIdx_act order

        %% flip phase based on CCEP waveforms
        toflip=[];

        if ~exist(['flipCCEP_' sb '.mat'], 'file')
            [phaseflipped, toflip] = ...
                flipPhase(phaseDat, zccep_clean(filterIdx,:), peakTimeMat(filterIdx,:), toflip);
            %======================================
            save(['flipCCEP_' sb '.mat'], 'toflip');
            %======================================
        else
            load(['flipCCEP_' sb '.mat'], 'toflip')
            phaseflipped = ...
                flipPhase(phaseDat, [], [], toflip);
        end


        %% plot to see different categories' spectra
        freq = cell(length(power.freqs),1);
        for i = 1:length(power.spectrum)
            freq{i} = sprintf('%.1f', power.freqs(i));
        end

        time = power.time;

        %% sort stim pairs
        THAL = {'antTH', 'midTH', 'pstTH'};
        SUBCOR = {'antTH', 'midTH', 'pstTH', 'BG', 'AMY', 'CLT'};

        fromTHAL = ismember(subT.JP_label_out, THAL) ;
        fromCOR = ~ismember(subT.JP_label_out, SUBCOR)  ;
        toTHAL = ismember(subT.JP_label_in, THAL)  ;
        toCOR = ~ismember(subT.JP_label_in, SUBCOR)  ;

        %% For different types of connectivity

        %-----------------------------------------------
        close all
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for it = 1:4

            type = connTypes(it,:);

            eval(sprintf('from = from%s;', type{1}));
            eval(sprintf('to = to%s;', type{2}));
           
            figure('Position', [1          41        1874         984]);
            sgtitle(sprintf('spectral CCEP from %s to %s of %s (%s)',type{1}, type{2}, strrep(sb, '_', '\_'), hcat));
            % separate between ACTIVATION == 1 AND == 0
            %-----------------------------------------------
            for ia = 1:3

                act = ia-2;

                %                 actIdx = (subT.CECS_activation == act & ...
                %                     subT.activated_default == act &...
                %                     subT.activated_SimRule == act);
                actIdx = actlab == act;
                UMAP_Activation = [UMAP_Activation; actlab(actlab==act)];
                %-------------------------------------
                % slice data
                sub_filterIdx = from & to & actIdx;
                metaT_filterIdx = [metaT_filterIdx; filterIdx(sub_filterIdx)];

                % log transform to power
                pw = squeeze(mean(log10(abs(powerDat(sub_filterIdx, :,:))), 1, 'omitnan'))';
                % no transform to phase
                ph = squeeze(mean(phaseflipped(sub_filterIdx, :,:),  1, 'omitnan'))';
                % sqrt transform to itpc
                pc = squeeze(mean(sqrt(itpcDat(sub_filterIdx, :,:)),  1, 'omitnan'))';

                %------------------------
                if act == 0
                    act_cat = 'umap-nonactivated';
                    pw0 = pw;
                    ph0 = ph;
                    pc0 = pc;
                elseif act == 1
                    act_cat = 'umap-activated';
                elseif act == -1
                    act_cat = 'umap-excluded';
                end

                %------------------------
                % RECORD FOR GROUP
                PW(:,:, ia, it, is) = pw; % freq*time(concatenated over ia, it, is)
                PH(:,:, ia, it, is) = ph;
                PC(:,:, ia, it, is) = pc;
                %------------------------
                %%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%
                % power
                n1 = ia;
                subplot(3,4,n1)

                d = pw;
                ttl = sprintf('power (%s)', act_cat);
                cmap = myColors('plasma',[],[],[],0);
                cmap = cmap.cm;
                labelstring = 'zscore-log(power)';
                subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring)

                % phase
                subplot(3,4,n1+4)
                d = ph;
                ttl = sprintf('phase (%s)', act_cat);
                cmap = flipud(cbrewer2('RdBu'));
                labelstring = 'zscore-phase flipped';
                subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring)

                % itpc
                subplot(3,4,n1+4*2)
                d = pc;
                ttl = sprintf('inter-trial phase coherence (%s)', act_cat);
                cmap = myColors('viridis',[],[],[],0);
                cmap = cmap.cm;
                labelstring = 'zscore-sqrt(itpc)';
                subplotCCEPspectrum(d, ttl,time,freq(1:end-1),cmap,labelstring)

            end

            % add one subplot panel for contrast
            % power
            n1 = 4;
            subplot(3,4,n1)

            d = pw-pw0;
            ttl = sprintf('power (%s)', 'act-nonact');
            cmap = myColors('plasma',[],[],[],0);
            cmap = cmap.cm;
            labelstring = 'zscore-log(power)';
            subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring)

            % phase
            n = n+1;
            subplot(3,4,n1+4)
            d = ph-ph0;
            ttl = sprintf('phase (%s)', 'act-nonact');
            cmap = flipud(cbrewer2('RdBu'));
            labelstring = 'zscore-phase flipped';
            subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring)

            % itpc
            n = n+1;
            subplot(3,4,n1+4*2)
            d = pc-pc0;
            ttl = sprintf('inter-trial phase coherence (%s)', 'contrast');
            cmap = myColors('viridis',[],[],[],0);
            cmap = cmap.cm;
            labelstring = 'zscore-sqrt(itpc)';
            subplotCCEPspectrum(d, ttl,time,freq(1:end-1),cmap,labelstring)

            figname = sprintf('%s_sbjCCEP_anatomyCAT_spectrogram_%s_%s-%s.png',hcat, sb, type{1}, type{2});
            exportgraphics(gcf, fullfile(plot_folder, figname))
        end % for type
    end % for is
    % save group-level spectrum
    save(fullfile(data_folder_intermediate, ...
        ['original_', hcat, '_CCEP_anatomyCAT_spectrogram_group.mat']), ...
        'PW','PH','PC','metaT_filterIdx','UMAP_Activation', '-v7.3');

    %% SORT GROUP LEVEL
    disp('Plot for group level spectrogram....')

    % take average
    PW_m = squeeze(mean(PW, 5, 'omitnan'));
    PH_m = squeeze(mean(PH, 5, 'omitnan'));
    PC_m = squeeze(mean(PC, 5, 'omitnan'));

    % plot for group level spectrogram
    for it = 1:length(connTypes)
        type = connTypes(it,:); 

        figure('Position', [1          41        1874         984]);
        sgtitle(sprintf('spectral CCEP from %s to %s of %s (%s)',type{1}, type{2}, 'Group average', hcat));

        for ia = 1:3
            pw = squeeze(PW_m(:,:,ia,it));
            ph = squeeze(PH_m(:,:,ia,it));
            pc = squeeze(PC_m(:,:,ia,it));
            act = ia-2;
            if act == 0
                act_cat = 'umap-nonactivated';
                pw0 = pw;
                ph0 = ph;
                pc0 = pc;
            elseif act == 1
                act_cat = 'umap-activated';
            elseif act == -1
                act_cat = 'umap-excluded';

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % power
            
            n1 = ia;
            subplot(3,4,n1)

            d = pw;
            ttl = sprintf('power (%s)', act_cat);
            cmap = myColors('plasma',[],[],[],0);
            cmap = cmap.cm;
            labelstring = 'zscore-log(power)';
            subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring)

            % phase
            subplot(3,4,n1+4)
            d = ph;
            ttl = sprintf('phase (%s)', act_cat);
            cmap = flipud(cbrewer2('RdBu'));
            labelstring = 'zscore-phase flipped';
            subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring)

            % itpc
            subplot(3,4,n1+2*4)
            d = pc;
            ttl = sprintf('inter-trial phase coherence (%s)', act_cat);
            cmap = myColors('viridis',[],[],[],0);
            cmap = cmap.cm;
            labelstring = 'zscore-sqrt(itpc)';

            subplotCCEPspectrum(d, ttl,time,freq(1:end-1),cmap,labelstring)

        end

        % add one subplot panel for contrast
        % power
        n1 = 4;
        subplot(3,4,n1)

        d = pw-pw0;
        ttl = sprintf('power (%s)', 'act-nonact');
        cmap = myColors('plasma',[],[],[],0);
        cmap = cmap.cm;
        labelstring = 'zscore-log(power)';
        subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring)

        % phase
        subplot(3,4,n1+4)
        d = ph-ph0;
        ttl = sprintf('phase (%s)', 'act-nonact');
        cmap = flipud(cbrewer2('RdBu'));
        labelstring = 'zscore-phase flipped';
        subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring)

        % itpc
        n = n+1;
        subplot(3,4,n1+2*4)
        d = pc-pc0;
        ttl = sprintf('inter-trial phase coherence (%s)', 'act-nonact');
        cmap = myColors('viridis',[],[],[],0);
        cmap = cmap.cm;
        labelstring = 'zscore-sqrt(itpc)';
        subplotCCEPspectrum(d, ttl,time,freq(1:end-1),cmap,labelstring)

        figname = sprintf('%s_groupCCEP_anatomyCAT_spectrogram_%s-%s.png',hcat, type{1}, type{2});
        exportgraphics(gcf, fullfile(plot_folder, figname))
    end % for type
end
