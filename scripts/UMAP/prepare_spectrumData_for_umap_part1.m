% this script will reduce the dimension of the spectrum for UMAP learn
% part1: average trials to see if there is a general differences between
% activated vs non-activated, cor vs thal categorties
home_dir = '/data/dian';
addpath(genpath(fullfile(home_dir, 'Dropbox/scripts/external/cbrewer2/cbrewer2')));
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/personal_validation'))
addpath('~/scripts/my_functions')
addpath(genpath('~/scripts/Stanford/ThalamocoricalLoop-project'))
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/CECS_Pipeline/CCEP_pipeline/HFB plot/subfunctions'));

data_folder = '~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore5_locked';
result_folder = '~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore5_locked';
if ~exist(fullfile(data_folder,'render1'),'dir'); mkdir (fullfile(data_folder,'render1'));end
plot_folder = '~/Dropbox/Stanford_Matters/data/THAL/CCEP/Plots/explore5/render1-reduced';
if ~exist(plot_folder, 'dir'); mkdir(plot_folder);end
% load metaTable
metaT = readtable('~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore5_locked/table_CCEPnewpipOutput_wholebrain_anatomical_info_activationRedone.csv');
%
% load data
load('~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore5_locked/CCEP_all_flat_meanTr_cleaned.mat');%, 'CCEP_flat' each row corresponding to the metaT
%
sblist = unique(metaT.subject);
vars = metaT.Properties.VariableNames;
ttpIdx = find(contains(vars, 'pks_time'));
peakTimeMat = table2array(metaT(:,ttpIdx));

%filterIdx = load('~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore3_new2Sbjs/CCEP_all_flat_meanTr_cleaned.mat', 'badChan');
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
    PW = nan(11, 65, 2, length(connTypes), length(sblist));
    PH = nan(11, 65, 2, length(connTypes), length(sblist));
    PC = nan(10, 65, 2, length(connTypes), length(sblist));

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

        [filterIdx, preIdxsub, ~] = intersect(phase.idx_in_metaT, prefilterIdx);
        subT = metaT(filterIdx,:);
        powerDat = power.CCEP_power(preIdxsub,:,:);
        phaseDat = phase.CCEP_phase(preIdxsub,:,:);
        itpcDat  = itpc.CCEP_itpc(preIdxsub,:,:);

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
        % freqSeg = [0.5 5; 5 8; 8 15; 15 30; 30 70; 70 256];
        % timeSeg = [-30 10; 10 117; 117 290; 282 800];
        % kfreq = [2, 2, 2, 2, 2, 1];ktime = [5 25 20 15];
        [~, rfreq, rtime] = reduceSpectrum([], power.freqs,power.time*1000,...
            [2, 2, 2, 2, 2, 1], [5 25 20 15],0);

        freq = cell(length(rfreq),1);
        for i = 1:length(rfreq)
            freq{i} = sprintf('%.1f', rfreq(i));
        end

        time = rtime/1000;

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

            figure('Position', [206           5        1664         946]);
            sgtitle(sprintf('reduced spectral CCEP from %s to %s of %s (%s)',type{1}, type{2}, strrep(sb, '_', '\_'), hcat));
            % separate between ACTIVATION == 1 AND == 0
            %-----------------------------------------------
            for ia = 1:2

                act = ia-1;

                actIdx = (subT.CECS_activation == act & ...
                    subT.activated_default == act &...
                    subT.activated_SimRule == act);

                sub_filterIdx = from & to & actIdx;

                % log transform to power
                pw = squeeze(mean(log10(abs(powerDat(sub_filterIdx, :,:))), 1, 'omitnan'))';
                % no transform to phase
                ph = squeeze(mean(phaseflipped(sub_filterIdx, :,:),  1, 'omitnan'))';
                % sqrt transform to itpc
                pc = squeeze(mean(sqrt(itpcDat(sub_filterIdx, :,:)),  1, 'omitnan'))';
                %-------------------------------------
                % reduce dimension of the spectrogram
                freqs = power.freqs;
                times = power.time;

                freq_seg = [2, 2, 2, 2, 2, 1];
                time_seg = [5 25 20 15];
                freq_seg2 = [2, 2, 2, 2, 2, 0];

                rpw = reduceSpectrum(pw, ...
                    freqs, times*1000, freq_seg, time_seg);
                rph = reduceSpectrum(ph, ...
                    freqs, times*1000, freq_seg, time_seg);
                rpc = reduceSpectrum(pc, ...
                    freqs, times*1000, freq_seg2, time_seg);
               
                %------------------------
                if act == 0
                    act_cat = 'non-activated';
                    pw0 = rpw;
                    ph0 = rph;
                    pc0 = rpc;
                else
                    act_cat = 'activated';
                end

                %------------------------
                % RECORD FOR GROUP
                PW(:,:, ia, it, is) = rpw; % freq*time(concatenated over ia, it, is)
                PH(:,:, ia, it, is) = rph;
                PC(:,:, ia, it, is) = rpc;
                %------------------------
                %%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%
                % power
                subplot(3,3,3*act+1)

                d = rpw;
                ttl = sprintf('power (%s)', act_cat);
                cmap = myColors('plasma',[],[],[],0);
                cmap = cmap.cm;
                labelstring = 'zscore-log(power)';
                subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring,5,1)

                % phase
                subplot(3,3,3*act+2)
                d = rph;
                ttl = sprintf('phase (%s)', act_cat);
                cmap = flipud(cbrewer2('RdBu'));
                labelstring = 'zscore-phase flipped';
                subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring,5,1)

                % itpc
                subplot(3,3,3*act+3)
                d = rpc;
                ttl = sprintf('inter-trial phase coherence (%s)', act_cat);
                cmap = myColors('viridis',[],[],[],0);
                cmap = cmap.cm;
                labelstring = 'zscore-sqrt(itpc)';
                subplotCCEPspectrum(d, ttl,time,freq(1:end-1),cmap,labelstring,5,1)

            end

            % add one subplot panel for contrast
            % power
            subplot(3,3,7)

            d = rpw-pw0;
            ttl = sprintf('power (%s)', 'contrast');
            cmap = myColors('plasma',[],[],[],0);
            cmap = cmap.cm;
            labelstring = 'zscore-log(power)';
            subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring,5,1)

            % phase
            subplot(3,3,8)
            d = rph-ph0;
            ttl = sprintf('phase (%s)', 'contrast');
            cmap = flipud(cbrewer2('RdBu'));
            labelstring = 'zscore-phase flipped';
            subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring,5,1)

            % itpc
            subplot(3,3,9)
            d = rpc-pc0;
            ttl = sprintf('inter-trial phase coherence (%s)', 'contrast');
            cmap = myColors('viridis',[],[],[],0);
            cmap = cmap.cm;
            labelstring = 'zscore-sqrt(itpc)';
            subplotCCEPspectrum(d, ttl,time,freq(1:end-1),cmap,labelstring,5,1)

            figname = sprintf('%s_sbjCCEP_anatomyCAT_spectrogram_%s_%s-%s.png',hcat, sb, type{1}, type{2});
            exportgraphics(gcf, fullfile(plot_folder, figname))
        end % for type
    end % for is
    % save group-level spectrum
    save(fullfile(data_folder, 'render1',['reduced_', hcat, '_CCEP_anatomyCAT_spectrogram_group.mat']), 'PW','PH','PC', '-v7.3');

    %% SORT GROUP LEVEL
    disp('Plot for group level spectrogram....')

    % take average
    PW_m = squeeze(mean(PW, 5, 'omitnan'));
    PH_m = squeeze(mean(PH, 5, 'omitnan'));
    PC_m = squeeze(mean(PC, 5, 'omitnan'));

    % plot for group level spectrogram
    for it = 1:length(connTypes)
        type = connTypes(it,:);

        figure('Position', [206           5        1664         946]);
        sgtitle(sprintf('reduced spectral CCEP from %s to %s of %s (%s)',type{1}, type{2}, 'Group average', hcat));

        for ia = 1:2
            rpw = squeeze(PW_m(:,:,ia,it));
            rph = squeeze(PH_m(:,:,ia,it));
            rpc = squeeze(PC_m(:,:,ia,it));
            act = ia-1;
            if act == 0
                act_cat = 'non-activated';
                pw0 = rpw;
                ph0 = rph;
                pc0 = rpc;
            else
                act_cat = 'activated';

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % power
            subplot(3,3,3*act+1)

            d = rpw;
            ttl = sprintf('power (%s)', act_cat);
            cmap = myColors('plasma',[],[],[],0);
            cmap = cmap.cm;
            labelstring = 'zscore-log(power)';
            subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring,5,1)

            % phase
            subplot(3,3,3*act+2)
            d = rph;
            ttl = sprintf('phase (%s)', act_cat);
            cmap = flipud(cbrewer2('RdBu'));
            labelstring = 'zscore-phase flipped';
            subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring,5,1)

            % itpc
            subplot(3,3,3*act+3)
            d = rpc;
            ttl = sprintf('inter-trial phase coherence (%s)', act_cat);
            cmap = myColors('viridis',[],[],[],0);
            cmap = cmap.cm;
            labelstring = 'zscore-sqrt(itpc)';

            subplotCCEPspectrum(d, ttl,time,freq(1:end-1),cmap,labelstring,5,1)

        end

        % add one subplot panel for contrast
        % power
        subplot(3,3,7)

        d = rpw-pw0;
        ttl = sprintf('power (%s)', 'contrast');
        cmap = myColors('plasma',[],[],[],0);
        cmap = cmap.cm;
        labelstring = 'zscore-log(power)';
        subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring,5,1)

        % phase
        subplot(3,3,8)
        d = rph-ph0;
        ttl = sprintf('phase (%s)', 'contrast');
        cmap = flipud(cbrewer2('RdBu'));
        labelstring = 'zscore-phase flipped';
        subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring,5,1)

        % itpc
        subplot(3,3,9)
        d = rpc-pc0;
        ttl = sprintf('inter-trial phase coherence (%s)', 'contrast');
        cmap = myColors('viridis',[],[],[],0);
        cmap = cmap.cm;
        labelstring = 'zscore-sqrt(itpc)';
        subplotCCEPspectrum(d, ttl,time,freq(1:end-1),cmap,labelstring,5,1)

        figname = sprintf('%s_groupCCEP_anatomyCAT_spectrogram_%s-%s.png',hcat, type{1}, type{2});
        exportgraphics(gcf, fullfile(plot_folder, figname))
    end % for type
end
