% This code is inherited from prelim_CCEP_explore5_2.m of the SELF-project
%dl577@stanford.edu; Mar 7, 2023
% The input data are partly imported from SELF-project, and partly from
% _CCEP_Behzad_Neda

home_dir = 'C:\Users\lvdia'; %getuserdir;

server_root = '/mnt/neurology_jparvizi$';

addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/Cohort_Organization'))
addpath(fullfile(home_dir, 'Dropbox/scripts/my_functions'));
addpath(fullfile(home_dir, 'Dropbox/scripts/external'));
addpath(fullfile(home_dir, 'Dropbox/scripts/external/cbrewer2/cbrewer2'));
addpath(fullfile(home_dir, 'Dropbox/scripts/external/ColorBrewer'));
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/lbcn_preproc/vizualization'));
addpath(fullfile(home_dir, 'Dropbox/scripts/external/pca_ica'));
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/CECS_Pipeline/CCEP_pipeline/tools'));
addpath(fullfile(home_dir, 'Dropbox/scripts/JAKE_2/lbcn_preproc/freesurfer'))
addpath(fullfile(home_dir, 'Dropbox/scripts/external/tight_subplot'))
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/personal_validation'))

comp_root = '/data/dian/Working_projects/data';
result_folder = fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore1_dataSort_Denoise');
prelim_plot_folder = fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/CCEP/Plots/explore1');

sublist = {'S22_178_AF', 'S21_172_KS','S21_166_TM', 'S20_152_HT', 'S21_167_MQ', 'S21_169_BH', 'S21_171_MM',...
    'S19_137_AF', 'S21_165_WN', 'S22_176_LB', 'S22_177_JM', 'S22_182_DH', 'S22_183_CR', ...
    'S22_185_TW', 'S22_188_CB', 'S22_189_LMA', 'S22_190_AS', 'S22_191_KM', 'S22_192_LG', 'S22_193_AM', 'S23_194_PS', 'S23_195_MZ'};%, , 'S22_181_CB'};


%%
time = -199:500;
%%
if ~exist(result_folder, 'dir'); mkdir(result_folder); end
if ~exist(prelim_plot_folder, 'dir'); mkdir(prelim_plot_folder); end

% Make sure your are connected to CISCO and logged in the server

%% what analyses
flatten_CCEP        = 0;
plot_regional_CCEP  = 1;
plot_brain          = 0;

%% Load metatable
% load the csv with JP labels cleaned up
metaT = readtable(fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/COHORT/table_CCEPnewpipOutput_wholebrain_anatomical_info_importedData.csv'));
JP_label_out = metaT.JP_label_out1;
JP_label_in = metaT.JP_label_in1;

JP_label_out = ListSortAnatLabel_THAL(JP_label_out, 1);
JP_label_in = ListSortAnatLabel_THAL(JP_label_in, 1);

metaT.JP_label_in = JP_label_in;
metaT.JP_label_out = JP_label_out;

vnames = metaT.Properties.VariableNames;

%%
tic
if flatten_CCEP == 1

    block_names = unique(metaT.block_name);
    CCEP_flat =  nan(size(metaT,1), 2201); % each row corresponding to the metaT; using trial average
    CCEP_flat2 =  nan(size(metaT,1), 2201); % each row corresponding to the metaT; using trial median
    block_before = 'dummy';

    parfor ir = 1:size(metaT,1)

        block = metaT.block_name{ir};

        %if ~strcmp(block_before, block)

        sbj_name = char(unique(metaT.subject(strcmp(metaT.block_name, block))));
        ccep_dir = fullfile(comp_root, 'computed_data/CCEP', sbj_name);
        m = load(fullfile(comp_root, 'neuralData/originalData', sbj_name, ['global_CCEP_' sbj_name '_metadata.mat']));

        fname = fullfile(ccep_dir, [block '_CCEP.mat']);

        d = load(fname); % CCEP

        % end

        block_before = block;

        rec_idx = find(strcmp(d.CCEP.recordingchannel, metaT.record_chan{ir}));

        wave = squeeze(d.CCEP.wave(:,rec_idx,:))';
        wave_clean = nan(size(wave));

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
    save(fullfile(result_folder, 'CCEP_all_flat_meanTr.mat'), 'CCEP_flat', '-v7.3')
    save(fullfile(result_folder, 'CCEP_all_flat_medianTr.mat'), 'CCEP_flat2', '-v7.3')

end
toc

%%
%%   ==================================== PLOT REGION AVERATED CCEP ====================================
%% loop over self hot and self cold elecgtrodes
if plot_regional_CCEP == 1

    %load(fullfile(result_folder, 'CCEP_all_flat_meanTr.mat'));%, 'CCEP_flat' each row corresponding to the metaT
  %  [badChan, zccep_clean]   = detectCCEP_badChan(CCEP_flat, 1);
    zccep_Tr = smoothdata(zccep_clean', 'gaussian', 10)';
    goodChan = ~(isnan(mean(zccep_Tr,2, 'omitnan')));

    %% filter data

    filterIdx = find( goodChan &...
        metaT.MNIout_coord_1 .* metaT.MNIin_coord_1 > 0 &... % ipsilateral
        metaT.eudDist>5 & ...
        (metaT.min_pk_time > 10 ) &...
        ~ismember(metaT.JP_label_in, {'', 'empty', 'NAN', 'EXCLUDE', 'NA'}) & ...
        ~ismember(metaT.JP_label_out, {'', 'empty', 'NAN', 'EXCLUDE', 'NA'}) & ...
        (ismember(metaT.JP_label_in, {'antTH', 'midTH', 'pstTH'}) | ...
        ismember(metaT.JP_label_out, {'antTH', 'midTH', 'pstTH'}) ...
        ));

    T = metaT(filterIdx,:);
    jp_labels = unique(T.JP_label_in);

    vars = T.Properties.VariableNames;
    ttpIdx = find(contains(vars, 'pks_time_'));

    min_ttp_all = [];
    % sort out the minimum time to peak
    for i = 1:size(T,1)
        spkIdx = table2array(T(i, ttpIdx));
        spkIdx = spkIdx(~isnan(spkIdx))+200;

        if isempty(spkIdx)
            min_ttp_all(i,1) = nan;
        else
            min_ttp_all(i,1) = min(spkIdx)-200;
        end
    end

    %% sort CCEP waves
    peakTimeMat = table2array(T(:,ttpIdx));
    zccepall_flip   = flipSignROI(zccep_Tr(filterIdx,:), peakTimeMat);
    zccep = zccepall_flip';

    %% divide thamalus divisions
    pairTypes = {'antTH', 'brain'; 'midTH', 'brain'; 'pstTH','brain'};

    for p = 1:size(pairTypes,1)
        pair1 = pairTypes{p,1};
        pair2 = pairTypes{p,2};
        %% manually select partial data to present
        %  jp_labels = jp_labels(ismember(jp_labels, jp_include));
        ttls = {['outflow pathways of ' pair1], ['inflow pathways of ' pair1]};

        cmap0 = cbrewer2('Set2', length(jp_labels));
        %% plot
        figure('Position',[543.0000   50.1429  950.0000  969.7143]);
        for pathway_id = 1:2
            %% group based on JP labels

            jp_IDX      = cell(size(jp_labels));
            meanCCEP_JP = [];
            min_ttp_JP  = [];
            jp_ccep_collect = {};
            %
            for ij = 1:length(jp_IDX)

                if pathway_id == 1 % outflow
                    pairType = strjoin({pair1, pair2, '-'});
                    jp_IDX{ij} = find(strcmp(T.JP_label_in, jp_labels{ij}) & ...
                        strcmp(T.JP_label_out, pair1) ...
                        );
                elseif pathway_id == 2 % inflow
                    pairType = strjoin({pair2, pair1, '-'});
                    jp_IDX{ij} = find(strcmp(T.JP_label_out, jp_labels{ij}) & ...
                        strcmp(T.JP_label_in, pair1) ...
                        );
                end
                jp_ccep = zccep(jp_IDX{ij}, 1:length(time));

                peakTimeMat = table2array(T(jp_IDX{ij},ttpIdx));
                %group ROI, to flip sign
                [zccep_flip,C] = flipSignROI(jp_ccep, peakTimeMat);
                zccep_flip = zccep_flip';
                jp_ccep_collect{ij,1} = zccep_flip;
                % ==============================

                mn = mean(zccep_flip, 1,  'omitnan');
                meanCCEP_JP = [meanCCEP_JP; mn];
                se = std(zccep_flip,1,  'omitnan')./sqrt(sum(~isnan(zccep_flip))); % standard error

                min_ttp_JP = [min_ttp_JP; mean(min_ttp_all(jp_IDX{ij}), 'omitnan')];

            end
            %}
            %% order by the time of the first peak
            pkloc = [];
            for i = 1:size(meanCCEP_JP,1)
                [pks, locs] = findpeaks(meanCCEP_JP(i,:),'MinPeakHeight',1., 'MinPeakProminence',0.5);
                quickest = find(locs == min(locs));
                if isempty(quickest)
                    pkloc = [pkloc; nan, nan];
                else
                    pkloc = [pkloc; pks(quickest), locs(quickest)];
                end
            end
            % sort labels by time to first peak
            [~,id_ordered]=sort(pkloc(:,2));
            jp_fast = jp_labels(id_ordered);
            pkloc   = pkloc(id_ordered,:);
            cmap    = cmap0(id_ordered,:);
            ccep_jp = jp_ccep_collect(id_ordered);
            %{
        %% order based on the min_ttp from the meta table
        [~,id_ordered] = sort(min_ttp_JP);
        jp_fast = jp_labels(id_ordered);
        min_ttp = min_ttp_JP(id_ordered);
        cmap    = cmap0;%(id_ordered,:);
            %}

            %% replot
            jp_IDX = cell(size(jp_fast));

            subplot(1,2, (-1.*pathway_id)+3)
            lift_par = 8;

            for ij = 1:length(jp_IDX)

                if size(zccep_flip,1) <= 3
                    continue; 
                end

                zccep_flip = ccep_jp{ij};
                mn = mean(zccep_flip, 1,  'omitnan') + lift_par.*(ij-1);
                %     meanCCEP_JP = [meanCCEP_JP; mn];

                se = std(zccep_flip, 1,  'omitnan')./sqrt(sum(~isnan(zccep_flip))); % standard error

                hold on
                % for i = 1:size(jp_ccep,1)
                % plot(time, jp_ccep(i,:), 'Color', [cmap(ij,:), 0.2]);
                shadedErrorBar(time, mn, se,...
                    'lineProps',{'Color',cmap(ij,:),'LineWidth', 1.35});
                % ====== decorations ==============
                set(gca,'ytick',[])
                ylims = get(gca, 'ylim');
                xlims = get(gca, 'xlim');
                hold on

                % add annotation
                if ~isnan(pkloc(ij,1))
                    % if rely on the zccep data itself
                    x = pkloc(ij,2) + min(time);
                    y = pkloc(ij,1) + lift_par.*(ij-1) - 1.5;
                    % if rely on the meta table
                    %  x = min_ttp(ij);
                    %  y = mn(round(x - min(time))) - 0.5;
                    % add marker
                    plot(x,y+2,'v','MarkerFaceColor',cmap(ij,:)/1.1,'MarkerEdgeColor','none');
                else
                    x = 500;
                    y = mean(mn)-1;

                end
                % add anatomical labeling
                text(x,y, jp_fast{ij}, 'Color', cmap(ij,:)/1.2, 'Fontsize', 12)
                %add yline
                plot([xlims(1) xlims(2)], [lift_par.*(ij-1),lift_par.*(ij-1)], ':','Color','k')
                title(ttls{pathway_id})
                %end
                %plot(time, mean(jp_ccep,1, 'omitnan'), 'Color', cmap(ij,:));
                % add onset

            end % for jp label
            % ====== decorations ==============
            set(gca,'ytick',[])
            ylims = get(gca, 'ylim');
            % add xline
            plot([0 0],[ylims(1) ylims(2)],'Color','k')
            % add ruler
            plot([400 400],[-30 -35],'Color',[0.2 0.2 0.2], 'LineWidth', 1.1)
            plot([380 420],[-35 -35],':','Color',[0.2 0.2 0.2])
            text(415,-32 ,'= 5 z-scores', 'Color',[ 0.2 0.2 0.2],'fontsize',10)
            hold off
            xlabel('Time (ms)')
            ylabel('Z score')
            fname = sprintf('%s_time2peaks', pairType);
            saveas(gcf, fullfile(prelim_plot_folder, [fname '.fig']))
            saveas(gcf, fullfile(prelim_plot_folder, [fname '.png']))
        end
    end
end

%% plot brain
if plot_brain == 1
    %% ==================PREPARE DATA==================
      filterIdx = find( goodChan &...
        metaT.MNIout_coord_1 .* metaT.MNIin_coord_1 > 0 &... % ipsilateral
        metaT.eudDist>5 & ...
        (metaT.min_pk_time > 10 ) &...
        ~ismember(metaT.JP_label_in, {'', 'empty', 'NAN', 'EXCLUDE', 'NA'}) & ...
        ~ismember(metaT.JP_label_out, {'', 'empty', 'NAN', 'EXCLUDE', 'NA'}) & ...
        (ismember(metaT.JP_label_in, {'antTH', 'midTH', 'pstTH'}) | ...
        ismember(metaT.JP_label_out, {'antTH', 'midTH', 'pstTH'}) ...
        ));

    T = metaT(filterIdx,:);
    CCEP_flat_ = CCEP_flat(filterIdx,:);
   %  [badChan, zccep_clean]   = detectCCEP_badChan(CCEP_flat, 1);
    zccep_Tr = smoothdata(zccep_clean', 'gaussian', 10)';
    %-- preselect the stim-rec pairs on the same hemisphere
    idx_LL = find(T.MNIout_coord_1 <= 0 & T.MNIin_coord_1 <= 0);
    idx_RR = find(T.MNIout_coord_1 > 0 & T.MNIin_coord_1 > 0);

    T          = T([idx_LL; idx_RR],:);
    zccep      = zccep_Tr([idx_LL; idx_RR],:);
    % project left to right hemisphere
    T.MNIin_coord_1(T.MNIin_coord_1 <= 0)   = -1.*T.MNIin_coord_1(T.MNIin_coord_1 <= 0);
    T.MNIout_coord_1(T.MNIout_coord_1 <= 0) = -1.*T.MNIout_coord_1(T.MNIout_coord_1 <= 0);

    %----

    output_folder = fullfile(prelim_plot_folder, 'render1');
    if ~exist(output_folder, 'dir'); mkdir(output_folder); end

    SbjElec_rec =  cellstr(string(T.subject)+'-'+string(T.rc1));
    SbjElec_stim = cellstr(string(T.subject)+'-'+string(T.sc1));

    [~, iu, ~] = unique(SbjElec_rec,'legacy');

    mni_in  = [T.MNIin_coord_1(iu),...
        T.MNIin_coord_2(iu),...
        T.MNIin_coord_3(iu)];
    jplabel_in = T.JP_label_in(iu);
    %     mni_out = [T.MNIout_coord_1(iu), T.MNIout_coord_2(iu), T.MNIout_coord_3(iu)];

    %---- filter with "stim chan == hot/cold" indexing (referring to 'T') -----------
    hotIdx  = find(T.stim_hot == 1);
    coldIdx = find(T.stim_hot == -1);
    % reindexing
    istim = [hotIdx; coldIdx];
    hotreIdx  = 1:length(hotIdx);
    coldreIdx = length(hotIdx)+1 : length(istim);

    SbjElec_rec_fromPMC = SbjElec_rec(istim);
    % ccep_data = ccep_icItp; % use alternative metric
    % {
    ccep_data = abs(zscore(zccep(istim,1:800),0,'all')); % ccep_icItp

    %===== prepare the zccep values =====
    % delete apparent artifact
    ccep_data(ccep_data(:,(500-min(time)) )>10,:) = 0;
    ccep_data(max(ccep_data(:, (0-min(time)):(10-min(time)) )) >5,:) = 0;
    % zero non-active : threshold is 3
    ccep_data(ccep_data<3) = 0;
    %ccep_data(ccep_data>5) = log(ccep_data(ccep_data>5));
    %}
    %-------  create mean ccep over stim channels of each condition -----
    ccep_data_grouped = zeros(size(mni_in, 1), max(time)-min(time)+1, 2); % the thrid dimension is for condition (hot=1 & cold=2)

    for ie = 1:2
        if ie == 1
            data = ccep_data(hotreIdx,:); reclist = SbjElec_rec_fromPMC(hotreIdx);
        else
            data = ccep_data(coldreIdx,:); reclist = SbjElec_rec_fromPMC(coldreIdx);
        end
        for ri = 1:length(iu)
            u = iu(ri);
            rec_unique = SbjElec_rec{u};
            ccep_data_grouped(ri,:, ie) = mean(data(strcmp(reclist,rec_unique),:), 1, 'omitnan'); % using mean
            %{
            % using exponential 
            d = data(strcmp(reclist,rec_unique),:);
            % using exp(mean)
           d = mean(d, 1, 'omitnan');
           d(d>0) = exp(d(d>0))-1;
           ccep_data_grouped(ri,:, ie) = d;
            %}
            %             if isempty(find(d~=0, 1)) || length(find(d~=0)) == 1 || size(d,1) == 1
            %                 ccep_data_grouped(ri,:, ie) = zeros(1, size(d,2));
            %             else
            %
            %             ccep_data_grouped(ri,:, ie) = exp(mean(d, 1, 'omitnan')./(var(d, [], 1, 'omitnan')./sqrt(size(d,1)))); % using tscore of the logarithic value;
            %             end
        end

    end
    %
    clear data reclist d
    % brain surface data
    if strcmp(home_dir, '/Users/dianlyu')
        [cmcortexL.vert, cmcortexL.tri]=read_surf(fullfile('/Applications/freesurfer/subjects/fsaverage' , 'surf', 'lh.pial'));
        [cmcortexR.vert, cmcortexR.tri]=read_surf(fullfile('/Applications/freesurfer/subjects/fsaverage' , 'surf', 'rh.pial'));
    elseif strcmp(home_dir, '/Users/dl577')
        [cmcortexL.vert, cmcortexL.tri]=read_surf(fullfile('/Applications/freesurfer/7.2.0/subjects/fsaverage' , 'surf', 'lh.pial'));
        [cmcortexR.vert, cmcortexR.tri]=read_surf(fullfile('/Applications/freesurfer/7.2.0/subjects/fsaverage' , 'surf', 'rh.pial'));
    end
    %--------------------------------------------------------------
    %% ======================== PLOTTING ========================
    tic
    disp('Plotting brain activation in time....')
    %============================================================
    % local electrode transparency for each perspective
    pers_rule = readtable(fullfile(home_dir, 'Dropbox/Stanford_Matters/data/cohort_data/ELECTRODE_VISUALIZATION_RULE.csv'));
    %clim = [min(reshape(ccep_icItp,[],1)), max(reshape(ccep_icItp,[],1))];
    ccepdist_hot  = reshape(ccep_data_grouped(:,:,1), [],1);
    ccepdist_cold = reshape(ccep_data_grouped(:,:,2), [],1);
    % pure difference
    ccepdist_diff = ccepdist_hot - ccepdist_cold;
    %     % relative change to cold CCEP
    %     ccepdist_diff = (ccepdist_hot - ccepdist_cold) ./ abs(ccepdist_cold) ;

    clim = [0, min([quantile(ccepdist_hot(ccepdist_hot>0), 0.9),...
        quantile(ccepdist_cold(ccepdist_cold>0), 0.9)])];

    % clim_diff = [-0.6 .* max(abs(ccepdist_diff), [], 'omitnan'), ...
    %    0.6 .* max(abs(ccepdist_diff), [], 'omitnan')];
    clim_diff_max = max(abs(quantile(ccepdist_diff(ccepdist_diff<0), 0.8)),...
        quantile(ccepdist_diff(ccepdist_diff>0), 0.8));

    clim_diff = [-1.*clim_diff_max, clim_diff_max];

    %  clear ccepdist_hot ccepdist_cold ccepdist_diff

    tic

    for timepoint = -18:500 % in ms
        close
        if exist(fullfile(output_folder, sprintf('BrainHeatMap_Outflow_%dms.jpg', timepoint)),'file') == 2
            continue;
        end

        val1  = reshape(ccep_data_grouped(:, timepoint-min(time), 1), [],1);
        val2  = reshape(ccep_data_grouped(:, timepoint-min(time), 2), [],1);

        % add a third column of contrast
        vals = [val1 val2 val1-val2];
        % zero the small difference values
        %         val3 = val1-val2;
        %         val3( mean(val3, 'omitnan')-std(val3, 'omitnan') < vals | vals < mean(val3, 'omitnan') + std(val3, 'omitnan'))

        clear var1 var2

        hemi = {'right'}; %hemi = {'left', 'right'};
        pers = {'lateral', 'medial', 'ventral'};
        nRow  = 3;
        nCol  = length(hemi).*length(pers); % nCol for visual contents
        nColP = 5;%*4-2; % nCol for subplots to iterate through
        % sort display order
        desired_order = {...
            'right-medial1', 'right-lateral1', 'right-ventral1';...
            'right-medial2', 'right-lateral2', 'right-ventral2';...
            'right-medial3', 'right-lateral3', 'right-ventral3'};
        % iter_order = cell(nRow, nCol);
        plot_order = [ ...
            {[1 2]},          { [3 4]},         {[5]};...
            {[1 2]+nColP},    { [3 4]+nColP},   {[5]+nColP};...
            {[1 2]+2*nColP},  { [3 4]+2*nColP}, {[5]+2*nColP}...
            ];
        in = 0;
        N = cell(numel(plot_order),1);
        for iv = 1:3
            for h = 1:length(hemi)
                for ps = 1:length(pers)
                    in = in+1;
                    [ir,ic] = find(strcmp([hemi{h} '-' pers{ps} num2str(iv)], desired_order));
                    N{in} = plot_order{ir,ic};
                end
            end
        end

        % ---------------------------------------------------
        figure('Position',[84   533   720   794]);
        %figure;
        %[ha, pos] = tight_subplot(nRow,nColP);

        n=0;
        for iv = 1:size(vals,2)

            if iv < 3
                clear cmap_
                cmap_ = cbrewer2('seq','RdPu', 64, 'spline');
                clim_ = clim;
            else
                clear cmap_
                cmap_ = flipud(cbrewer2('div','RdBu', 64, 'spline'));
                clim_ = clim_diff;
            end

            % n=0;figure;
            for h = 1:length(hemi)

                hem = hemi{h};

                if strcmp(hemi, 'left')
                    hemiIdx = find(mni_in(:,1) < 0);
                    cmcortex = cmcortexL;
                elseif strcmp(hemi, 'right')
                    hemiIdx = find(mni_in(:,1) >= 0);
                    cmcortex = cmcortexR;
                end
                %hemiIdx = find((mni_in(:,1) .* ((h.*2)-3))>0);

                electrodes = mni_in(hemiIdx,:);
                JP_list    = jplabel_in(hemiIdx,:);
                weight     = vals(hemiIdx,iv);
                weight(isnan(weight) ) = 0;
                %weight = ones(size(weight));

                for ps = 1:length(pers)

                    per = pers{ps};
                    %figure('Position', [562   667   711   540])
                    n = n + 1;
                    subplot(nRow,nColP,N{n})
                    % axes(ha(N{n}));
                    % figure;
                    %
                    if ismember([hem '-' per num2str(iv)], desired_order(:,end))
                        show_colorbar = 1;
                    else
                        show_colorbar = 0;
                    end
                    %                    show_colorbar = 0;
                    %figure;
                    covgAdj = 1/2; % when two hemispheres are projected to one
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    [elec_surf] = ctmr_gauss_plot3(cmcortex, electrodes, weight, cmap_, hem, per, 1, 80, clim_, 1, show_colorbar, covgAdj);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if strcmp([hem '-' per num2str(iv)], desired_order{1,1})
                        text(0,100, 70, 'Hot', 'FontSize',18)
                    elseif strcmp([hem '-' per num2str(iv)], desired_order{2,1})
                        text(0,100, 70, 'Cold', 'FontSize',18)
                    elseif strcmp([hem '-' per num2str(iv)], desired_order{3,1})
                        text(0,100, 70, 'Contrast', 'FontSize',18)
                    end
                    alpha(0.75)

                    % sort the explicit electrodes for this perspective
                    jpSHOW    = pers_rule.JP_label_std(strcmpi(pers_rule.perspective, per));
                    elec_surfSHOW  = elec_surf(ismember(JP_list, jpSHOW),:);

                    %% plot individual electrodes
                    hold on
                    recElecDotSize = 15;
                    %p0 = scatter3(electrodes(:,1), electrodes(:,2), electrodes(:,3),  5, ...
                    p0 = scatter3(elec_surf(:,1),elec_surf(:,2), elec_surf(:,3),  recElecDotSize, ...
                        'MarkerEdgeColor', [0.25 0.25 0.25], 'MarkerFaceColor', [0.65 0.65 0.65]);
                    p0.MarkerFaceAlpha = 0.0; %.5;
                    p0.MarkerEdgeAlpha = 0.0; %.85;
                    %---- make the local electrodes on the correct perspective more
                    %transparent ----
                    p0_1 = scatter3(elec_surfSHOW(:,1),elec_surfSHOW(:,2), elec_surfSHOW(:,3),  recElecDotSize, ...
                        'MarkerEdgeColor', [0.25 0.25 0.25], 'MarkerFaceColor', [0.65 0.65 0.65]);
                    p0_1.MarkerFaceAlpha = .7; %.5;
                    p0_1.MarkerEdgeAlpha = .95; %.85;
                    %=========================== activated elec ===================================
                    % litup_colors = %[[255, 145, 172]/256; [255, 145, 172]/256]; % pink
                    litup_colors = [[166,27,41]/256; [255, 203, 1]/256;  [35, 118, 183]/256];% dark red[160, 28, 52],yellow, blue,
                    %litup_colors = [ [188, 9, 45]/256; [188, 9, 45]/256];
                    %---> if the ccep_data is ccep_icItp
                    % reweights = weights; reweights(reweights<0.0001) = 0;
                    % % ---> if the ccep_data is zccep
                    implicit_FaceAlpha = 0.0;
                    implicit_EdgeAlpha = 0.0;
                    explicit_FaceAlpha = 0.25;
                    explicit_EdgeAlpha = 0.65;

                    maxdotsize  = 250;
                    weightSHOW = weight(ismember(JP_list, jpSHOW),:);

                    if iv ~= 3
                        scp = (max(clim_)- min(clim_)) ./maxdotsize ;
                        reweight = (weight - min(clim_))/scp + 0.01;
                        reweightSHOW = reweight(ismember(JP_list, jpSHOW),:);

                        p1 = scatter3(elec_surf(:,1), elec_surf(:,2), elec_surf(:,3), reweight,  ...
                            'MarkerEdgeColor', litup_colors(1,:), ...%repmat(litup_color, size(weights,1),1) weights(:,1)] ,...
                            'MarkerFaceColor', litup_colors(1,:), ...
                            'AlphaData', reweight/0.5);%[repmat(litup_color, size(weights,1),1) weights(:,1)] );
                        p1.MarkerFaceAlpha = implicit_FaceAlpha;
                        p1.MarkerEdgeAlpha = implicit_EdgeAlpha;

                        %-------
                        p1_1 = scatter3(elec_surfSHOW(:,1), elec_surfSHOW(:,2), elec_surfSHOW(:,3), reweightSHOW,  ...
                            'MarkerEdgeColor', litup_colors(1,:), ...%repmat(litup_color, size(weights,1),1) weights(:,1)] ,...
                            'MarkerFaceColor', litup_colors(1,:), ...
                            'AlphaData', reweightSHOW/0.5);%[repmat(litup_color, size(weights,1),1) weights(:,1)] );
                        p1_1.MarkerFaceAlpha = explicit_FaceAlpha;
                        p1_1.MarkerEdgeAlpha = explicit_EdgeAlpha;
                    else

                        scp = clim_(2) ./maxdotsize ;
                        reweight = abs(weight)/scp + 0.01;
                        reweightSHOW = reweight(ismember(JP_list, jpSHOW),:);

                        p1 = scatter3(elec_surf(weight>0, 1), elec_surf(weight>0, 2), elec_surf(weight>0, 3), reweight(weight>0),  ...
                            'MarkerEdgeColor', litup_colors(2,:), ...%repmat(litup_color, size(weights,1),1) weights(:,1)] ,...
                            'MarkerFaceColor', litup_colors(2,:), ...
                            'AlphaData', reweight/0.5);%[repmat(litup_color, size(weights,1),1) weights(:,1)] );
                        p1.MarkerFaceAlpha = implicit_FaceAlpha;
                        p1.MarkerEdgeAlpha = implicit_EdgeAlpha;
                        %-----------------------
                        p1_1 = scatter3(elec_surfSHOW(weightSHOW>0, 1), elec_surfSHOW(weightSHOW>0, 2), elec_surfSHOW(weightSHOW>0, 3), reweightSHOW(weightSHOW>0),  ...
                            'MarkerEdgeColor', litup_colors(2,:), ...%repmat(litup_color, size(weights,1),1) weights(:,1)] ,...
                            'MarkerFaceColor', litup_colors(2,:), ...
                            'AlphaData', reweightSHOW/0.5);%[repmat(litup_color, size(weights,1),1) weights(:,1)] );
                        p1_1.MarkerFaceAlpha = explicit_FaceAlpha;
                        p1_1.MarkerEdgeAlpha = explicit_EdgeAlpha;
                        %=========================================================
                        p2 = scatter3(elec_surf(weight<0, 1), elec_surf(weight<0, 2), elec_surf(weight<0, 3), reweight(weight<0),  ...
                            'MarkerEdgeColor', litup_colors(3,:), ...%repmat(litup_color, size(weights,1),1) weights(:,1)] ,...
                            'MarkerFaceColor', litup_colors(3,:), ...
                            'AlphaData', reweight/0.5);%[repmat(litup_color, size(weights,1),1) weights(:,1)] );
                        p2.MarkerFaceAlpha = implicit_FaceAlpha;
                        p2.MarkerEdgeAlpha = implicit_EdgeAlpha;
                        %-----------------------
                        p2_1 = scatter3(elec_surfSHOW(weightSHOW<0, 1), elec_surfSHOW(weightSHOW<0, 2), elec_surfSHOW(weightSHOW<0, 3), reweightSHOW(weightSHOW<0),  ...
                            'MarkerEdgeColor', litup_colors(3,:), ...%repmat(litup_color, size(weights,1),1) weights(:,1)] ,...
                            'MarkerFaceColor', litup_colors(3,:), ...
                            'AlphaData', reweightSHOW/0.5);%[repmat(litup_color, size(weights,1),1) weights(:,1)] );
                        p2_1.MarkerFaceAlpha = explicit_FaceAlpha;
                        p2_1.MarkerEdgeAlpha = explicit_EdgeAlpha;

                    end
                    %                    reweights(reweights<5) = 0;
                    %reweights(reweights>5) = log(reweights(reweights>5));

                    sgtitle(sprintf('       Outflow pathway (%d ms)\n   ', timepoint), 'FontSize',12)

                    hold off
                end
            end
        end
        exportgraphics(gcf,fullfile(output_folder, sprintf('BrainHeatMap_Outflow_%dms.jpg', timepoint)),'Resolution',320)
    end % for time point

end
toc
%ctmr_gauss_plot(cmcortex,[0 0 0], 0, hem, per)

%
%         alpha(0.8)
%
%         hold off
%
%         figname = sprintf('Total_elec_coverage_SelfCohort_%s_%s.jpg', per, hem);
%exportgraphics(gcf, fullfile(figpath, figname),'Resolution', 300)

%if ps == 2 && h == 2
%             legend('show');
%end

%     %% check network affiliation
%
%     % 1. stimulated from PMC
%     pmcT  = T(strcmp(T.JP_label_out_1, 'PMC') | strcmp(T.JP_label_out_2, 'PMC'),:);% | ...
%             %contains(lower(T.Destr_long_out), 'PCC'),:);
%     thalT = T(contains(lower(T.JP_label_out_1), 'thalam') | contains(lower(T.JP_label_out_2), 'thalam'),:);%| ...
%              % contains(lower(T.Destr_long_out), 'thalam'), :);
%
%