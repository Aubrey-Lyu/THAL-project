% this script will take the CCEP multi-freq power, phase and ITPC as input,
% to plot and compare between different brain areas.
home_dir = '/data/dian';
addpath(genpath(fullfile(home_dir, 'Dropbox/scripts/external/cbrewer2/cbrewer2')));
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/personal_validation'))
addpath('~/scripts/my_functions')
addpath(genpath('~/scripts/Stanford/ThalamocoricalLoop-project'))

result_folder = '~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore4_spectrum';
plot_folder = '~/Dropbox/Stanford_Matters/data/THAL/CCEP/Plots/explore4/render3-ActRelabel';
if ~exist(plot_folder, 'dir'); mkdir(plot_folder); end

% load metaTable
metaT = readtable('~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore3_new2Sbjs/table_CCEPnewpipOutput_wholebrain_anatomical_info_activationRedone.csv');
% load data
load('~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore3_new2Sbjs/CCEP_all_flat_meanTr_cleaned.mat');%, 'CCEP_flat' each row corresponding to the metaT
 
sblist = unique(metaT.subject);
 vars = metaT.Properties.VariableNames;
ttpIdx = find(contains(vars, 'pks_time'));
peakTimeMat = table2array(metaT(:,ttpIdx));

%filterIdx = load('~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore3_new2Sbjs/CCEP_all_flat_meanTr_cleaned.mat', 'badChan');
goodChan = ones(size(metaT,1),1);
goodChan(badChan) = 0;
cd(fullfile(result_folder, 'spectCCEP'));

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

if exist(['flipCCEP_' sb '.mat'], 'file')
    load(['flipCCEP_' sb '.mat'], 'toflip')
    phaseflipped = ...
        flipPhase(powerDat, zccep_clean(filterIdx,:), peakTimeMat(filterIdx,:), toflip);
else
    [phaseflipped, toflip] = ...
        flipPhase(powerDat, zccep_clean(filterIdx,:), peakTimeMat(filterIdx,:), toflip);
    save(['flipCCEP_' sb '.mat'], 'toflip', 'filterIdx');
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

            figure('Position', [206           5        1664         946]);
            sgtitle(sprintf('CCEP from %s to %s of %s (%s)',type{1}, type{2}, strrep(sb, '_', '\_'), hcat));
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
                ph = squeeze(mean(phaseflipped(sub_filterIdx, :,:), 1, 'omitnan'))';
                % sqrt transform to itpc
                pc = squeeze(mean(sqrt(itpcDat(sub_filterIdx, :,:)), 1, 'omitnan'))';
                %------------------------
                if act == 0
                    act_cat = 'non-activated';
                    pw0 = pw;
                    ph0 = ph;
                    pc0 = pc;
                else
                    act_cat = 'activated';
                end

                %------------------------
                % RECORD FOR GROUP
                PW(:,:, ia, it, is) = pw; % freq*time(concatenated over ia, it, is)
                PH(:,:, ia, it, is) = ph;
                PC(:,:, ia, it, is) = pc;
                %------------------------
                %%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%
                % power
                subplot(3,3,3*act+1)

                d = pw;
                ttl = sprintf('power (%s)', act_cat);
                cmap = myColors('plasma',[],[],[],0);
                cmap = cmap.cm;
                labelstring = 'log(power)';
                subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring)

                % phase
                subplot(3,3,3*act+2)
                d = ph;
                ttl = sprintf('phase (%s)', act_cat);
                cmap = flipud(cbrewer2('RdBu'));
                labelstring = 'phase flipped';
                subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring)

                % itpc
                subplot(3,3,3*act+3)
                d = pc;
                ttl = sprintf('inter-trial phase coherence (%s)', act_cat);
                cmap = myColors('viridis',[],[],[],0);
                cmap = cmap.cm;
                labelstring = 'sqrt(itpc)';
                subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring)

            end

            % add one subplot panel for contrast
            % power
            subplot(3,3,7)

            d = pw-pw0;
            ttl = sprintf('power (%s)', 'contrast');
            cmap = myColors('plasma',[],[],[],0);
            cmap = cmap.cm;
            labelstring = 'log(power)';
            subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring)

            % phase
            subplot(3,3,8)
            d = ph-ph0;
            ttl = sprintf('phase (%s)', 'contrast');
            cmap = flipud(cbrewer2('RdBu'));
            labelstring = 'phase flipped';
            subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring)

            % itpc
            subplot(3,3,9)
            d = pc-pc0;
            ttl = sprintf('inter-trial phase coherence (%s)', 'contrast');
            cmap = myColors('viridis',[],[],[],0);
            cmap = cmap.cm;
            labelstring = 'sqrt(itpc)';
            subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring)

            figname = sprintf('%s_sbjCCEP_anatomyCAT_spectrogram_%s_%s-%s.png',hcat, sb, type{1}, type{2});
            exportgraphics(gcf, fullfile(plot_folder, figname))
        end % for type
    end % for is
  % save group-level spectrum
    save(fullfile(result_folder, [hcat, '_CCEP_anatomyCAT_spectrogram_group.mat']), 'PW','PH','PC', '-v7.3');

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
        sgtitle(sprintf('CCEP from %s to %s of %s (%s)',type{1}, type{2}, 'Group average', hcat));

        for ia = 1:2
            pw = squeeze(PW_m(:,:,ia,it));
            ph = squeeze(PH_m(:,:,ia,it));
            pc = squeeze(PC_m(:,:,ia,it));
            act = ia-1;
            if act == 0
                act_cat = 'non-activated';
                pw0 = pw;
                ph0 = ph;
                pc0 = pc;
            else
                act_cat = 'activated';

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % power
            subplot(3,3,3*act+1)

            d = pw;
            ttl = sprintf('power (%s)', act_cat);
            cmap = myColors('plasma',[],[],[],0);
            cmap = cmap.cm;
            labelstring = 'log(power)';
            subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring)
            
            % phase
            subplot(3,3,3*act+2)
            d = ph;
            ttl = sprintf('phase (%s)', act_cat);
            cmap = flipud(cbrewer2('RdBu'));
            labelstring = 'phase flipped';
            subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring)

            % itpc
            subplot(3,3,3*act+3)
            d = pc;
            ttl = sprintf('inter-trial phase coherence (%s)', act_cat);
            cmap = myColors('viridis',[],[],[],0);
            cmap = cmap.cm;
            labelstring = 'sqrt(itpc)';

            subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring)
           
        end

        % add one subplot panel for contrast
        % power
        subplot(3,3,7)

        d = pw-pw0;
        ttl = sprintf('power (%s)', 'contrast');
        cmap = myColors('plasma',[],[],[],0);
        cmap = cmap.cm;
        labelstring = 'log(power)';
        subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring)

        % phase
        subplot(3,3,8)
        d = ph-ph0;
        ttl = sprintf('phase (%s)', 'contrast');
        cmap = flipud(cbrewer2('RdBu'));
        labelstring = 'phase flipped';
       subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring)

        % itpc
        subplot(3,3,9)
        d = pc-pc0;
        ttl = sprintf('inter-trial phase coherence (%s)', 'contrast');
        cmap = myColors('viridis',[],[],[],0);
        cmap = cmap.cm;
        labelstring = 'sqrt(itpc)';
        subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring)

        figname = sprintf('%s_groupCCEP_anatomyCAT_spectrogram_%s-%s.png',hcat, type{1}, type{2});
        exportgraphics(gcf, fullfile(plot_folder, figname))

    end % for type

end % for hemisphere