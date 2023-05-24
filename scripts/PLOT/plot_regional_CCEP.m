function  plot_regional_CCEP(zccep_clean, metaT, ROIs, filterIdx, plot_folder, cfg)
% dl577@stanford.edu, March 31 2023
% 1. the plotting function requires cleaned zccep matrix with a row,
% the figure has two subplot, the subplot(1,2,1) has the inflow pathway,
% and the subplot(1,2,2) has the outflow pathway
% number (pair number) same as the meta table (metaT)
% 2. metaT needs to have metaT.JP_label_in = JP_label_in;
%metaT.JP_label_out colmns sorted.
% 3. dimension of zccep_clean : pairsID * timepoints
% 4.The ROIs can be a cell or characters, it will be the region of interest to
% examine their inflow and outflow pathway among all the JP_labelled
% regions given in metaT.
% 5. cfg is a structure with some parameters to manipulate the data
% cfg.thr is the threshold used to detect ROI-CCEP peaks
% 6. dependencies:
%{
addpath(fullfile(home_dir, 'Dropbox/scripts/external/cbrewer2/cbrewer2'));
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/lbcn_preproc/vizualization'));
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/CECS_Pipeline/CCEP_pipeline/tools'));
addpath(fullfile(home_dir, 'Dropbox/scripts/external/tight_subplot'))
addpath(fullfile(home_dir, 'Dropbox/scripts/Stanford/CECS_Pipeline_COPY/personal_validation'))
%}
%% define default
zccep_Tr = smoothdata(zccep_clean', 'gaussian', 10)';
goodChan = ~(isnan(mean(zccep_Tr,2, 'omitnan')));

if ischar(ROIs)
    ROIs = strcell(ROIs);
end
pairTypes = cell(length(ROIs),2);
for ir = 1:length(ROIs)
    pairTypes(ir,:) = {ROIs{ir}, 'brain'};
end

vars = metaT.Properties.VariableNames;
ttpIdx = find(contains(vars, 'pks_time'));

if nargin<4 || isempty(filterIdx)
    artpeakearlytime = 10; % assumption
    % filter data: use activated ccep, Euclidean distance>5, ipsilateral pairs
    disp('Use activated CCEP (with time-to-peaks>10 ms), Euclidean distance>5mm, ipsilateral pairs...')
    filterIdx = find( goodChan &...
        metaT.MNIout_coord_1 .* metaT.MNIin_coord_1 > 0 &... % ipsilateral
        metaT.eudDist>5 & ...
        (sum(table2array(metaT(:, ttpIdx)), 2, 'omitnan') ) > artpeakearlytime   &... % exclude case when the only peak happens ealier than 10 ms
        ~ismember(metaT.JP_label_in, {'', 'empty', 'NAN', 'EXCLUDE', 'NA'}) & ...
        ~ismember(metaT.JP_label_out, {'', 'empty', 'NAN', 'EXCLUDE', 'NA'}) & ...
        (ismember(metaT.JP_label_in, {'antTH', 'midTH', 'pstTH'}) | ...
        ismember(metaT.JP_label_out, {'antTH', 'midTH', 'pstTH'}) ...
        ));
end

if nargin<5
    warning('No output directory is given, using the current directory')
    plot_folder = '.';
end

default_peakThr = 1.2;
default_promThr = 0.5;

if nargin<6
    cfg.lift_par = 12;
    cfg.markerAddY = 2.5;
    cfg.peakThr = default_peakThr;
    cfg.promThr = default_promThr;
    cfg.time = -199:500;
    cfg.art_win = [-200:artpeakearlytime, 800:2000]; % default: not counting peaks before 10 ms and after 800 ms
    fprintf('\n>>> Using roi-ccep peak detect thresholds: peak > %.2f, prominance > %.2f\n',cfg.peakThr,cfg.promThr )
else
    if isfield('cfg', 'time')
        cfg.time = -199:500;
    end
    if isfield('cfg', 'lift_par')
        cfg.lift_par = 12;
    end
    if isfield('cfg', 'markerAddY')
        cfg.markerAddY = 4;
    end
    if isfield('cfg', 'peakThr')
        cfg.peakThr = default_peakThr;
    end
    if isfield('cfg', 'promThr')
        cfg.promThr = default_promThr;
    end
end

%% main function

time = cfg.time;

T = metaT(filterIdx,:);
jp_labels0 = unique(T.JP_label_in);

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
zccep   = flipSignROI(zccep_Tr(filterIdx,:), peakTimeMat);


%% divide thamalus divisions
%   pairTypes = {'antTH', 'brain'; 'midTH', 'brain'; 'pstTH','brain'};

for p = 1:size(pairTypes,1)
    pair1 = pairTypes{p,1};
    pair2 = pairTypes{p,2};
    %% manually select partial data to present
    %  jp_labels = jp_labels(ismember(jp_labels, jp_include));
    ttls = {['outflow pathways of ' pair1], ['inflow pathways of ' pair1]};

    cmap0 = cbrewer2('Set2', length(jp_labels0));
    %% plot
    figure('Position',[543.0000   50.1429  950.0000  969.7143]);
    for pathway_id = 1:2
        %% group based on JP labels

        jp_IDX      = cell(size(jp_labels0));
        jp_labels   = {};
        meanCCEP_JP = []; seCCEP_JP = [];
        min_ttp_JP  = [];
        %
        for i = 1:length(jp_IDX)

            if pathway_id == 1 % outflow
                stimchan = pair1; 
                recordchan = pair2;
                pair_type = strjoin({stimchan, recordchan, '-'});
                jp_IDX{i} = find(strcmp(T.JP_label_in, jp_labels0{i}) & ...
                    strcmp(T.JP_label_out, stimchan) ...
                    );
            elseif pathway_id == 2 % inflow
                stimchan = pair2; 
                recordchan = pair1;
                pair_type = strjoin({stimchan, recordchan, '-'});
                jp_IDX{i} = find(strcmp(T.JP_label_out, jp_labels0{i}) & ...
                    strcmp(T.JP_label_in, recordchan) ...
                    );
            end
            if isempty(jp_IDX{i} )||length(jp_IDX{i})<3 || ...
              strcmp(jp_labels0{i}, pair1)   % does not allow self connection 
                continue; 
            end

            jp_labels{end+1,1} = jp_labels0{i};
            jp_ccep = zccep(jp_IDX{i}, time+200);

            peakTimeMat = table2array(T(jp_IDX{i},ttpIdx));
            %group ROI, to flip sign
            [zccep_flip,~] = flipSignROI(jp_ccep, peakTimeMat);

            % ==============================

            mn = mean(zccep_flip, 1,  'omitnan');
            meanCCEP_JP = [meanCCEP_JP; mn];
            se = std(zccep_flip,1,  'omitnan')./sqrt(sum(~isnan(zccep_flip))); % standard error
            seCCEP_JP = [seCCEP_JP; se];
            min_ttp_JP = [min_ttp_JP; mean(min_ttp_all(jp_IDX{i}), 'omitnan')];

        end
        %}
        %% order by the time of the first peak
        pkloc = [];  peakThr = cfg.peakThr; promThr = cfg.promThr;
        artwin = cfg.art_win - min(cfg.time);
        for i = 1:size(meanCCEP_JP,1)
            % when not flip
            [pks, locs] = findpeaks(meanCCEP_JP(i,:),'MinPeakHeight',peakThr, 'MinPeakProminence',promThr);
            locs(ismember(locs, artwin)) = []; pks(ismember(locs, artwin)) = []; %disregard the peaks happening in the first 15ms
            quickest = find(locs == min(locs, [], 'omitnan'));
            if isempty(quickest)
                pkloc_tmp = [nan, nan];
            else
                pkloc_tmp = [pks(quickest), locs(quickest)];
            end
            % when flipped
            [pks_, locs_] = findpeaks(-1 .* meanCCEP_JP(i,:),'MinPeakHeight',peakThr, 'MinPeakProminence',promThr);
            locs_(ismember(locs, artwin)) = []; pks_(ismember(locs, artwin)) = []; %disregard the peaks happening in the first 15ms
            quickest_ = find(locs_ == min(locs_, [], 'omitnan'));

            if isempty(quickest_)
                pkloc_tmp_ = [nan, nan];
            else
                pkloc_tmp_ = [pks_(quickest_), locs_(quickest_)];
            end

            pkloc(i,:) = pkloc_tmp; % default;
            if (pkloc_tmp_(2) < pkloc_tmp(2)) && ...
                    (~isnan(pkloc_tmp(2)) && ~isnan(pkloc_tmp_(2)))
                % use the flipped CCEP
                pkloc(i,:) = pkloc_tmp_;
                meanCCEP_JP(i,:) = -1 .* meanCCEP_JP(i,:);
            elseif isnan(pkloc_tmp(2)) && ~isnan(pkloc_tmp_(2))
                pkloc(i,:) = pkloc_tmp_;
                meanCCEP_JP(i,:) = -1 .* meanCCEP_JP(i,:);
            elseif isnan(pkloc_tmp(2)) && isnan(pkloc_tmp_(2))
                pkloc(i,:) = [nan, nan];
            end
        end
        % sort labels by time to first peak
        [~,id_ordered]=sort(pkloc(:,2));
        jp_fast = jp_labels(id_ordered);
        pkloc   = pkloc(id_ordered,:);
        meanCCEP_JP = meanCCEP_JP(id_ordered,:);
        seCCEP_JP = seCCEP_JP(id_ordered,:);
        % sort colors
        [~,~,cid_ordered] = intersect(jp_fast, jp_labels0, 'stable');
        cmap    = cmap0(cid_ordered,:);

        %{
        %% order based on the min_ttp from the meta table
        [~,id_ordered] = sort(min_ttp_JP);
        jp_fast = jp_labels(id_ordered);
        min_ttp = min_ttp_JP(id_ordered);
        cmap    = cmap0;%(id_ordered,:);
        %}

        %% replot

        subplot(1,2, (-1.*pathway_id)+3)
        lift_par= cfg.lift_par;
        for i = 1:size(jp_fast,1)

            mn = meanCCEP_JP(i,:)+lift_par.*(i-1);
            se = seCCEP_JP(i,:);

            hold on
            % for i = 1:size(jp_ccep,1)
            % plot(time, jp_ccep(i,:), 'Color', [cmap(ij,:), 0.2]);
            shadedErrorBar(time, mn, se,...
                'lineProps',{'Color',cmap(i,:),'LineWidth', 1.35});
            % ====== decorations ==============
            set(gca,'ytick',[])
            ylims = get(gca, 'ylim');
            xlims = get(gca, 'xlim');
            hold on

            % add annotation
            if ~isnan(pkloc(i,1))
                % if rely on the zccep data itself
                x = pkloc(i,2) + min(time);
                y = pkloc(i,1) + lift_par.*(i-1) + cfg.markerAddY ;
                % if rely on the meta table
                %  x = min_ttp(ij);
                %  y = mn(round(x - min(time))) - 0.5;
                % add marker
                plot(x,y,'v','MarkerFaceColor',cmap(i,:)/1.1,'MarkerEdgeColor','none');
            else
                x = 500;
                y = mean(mn)-1;

            end
            % add anatomical labeling
            text(x+5,y, jp_fast{i}, 'Color', cmap(i,:)/1.2, 'Fontsize', 12)
            %add yline
            plot([xlims(1) xlims(2)], [lift_par.*(i-1),lift_par.*(i-1)], ':','Color','k')
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
        fname = sprintf('%s_time2peaks', pair_type);

    end
    saveas(gcf, fullfile(plot_folder, [fname '.fig']))
    saveas(gcf, fullfile(plot_folder, [fname '.png']))
end