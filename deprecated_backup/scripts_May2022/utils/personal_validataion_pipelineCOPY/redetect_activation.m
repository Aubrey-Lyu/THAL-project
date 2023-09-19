function [activated_redo, act_time] = redetect_activation(metaT, CCEP_meanTr, threshold)
% this function will detect activation based on hard threshold based on the
% cleaned trial-averaged dataset
%activated_redo is a vector with 1 = activation detected, 0 = none
% metaT is the organized dataset from the meta table generated from the
% pipeline
% CCEP_meanTr is the preprocessed CCEP data structure concatenated based on
% the group level metatable, with the dimension = sites * timepoints
% threshold is a structure with the following default

if nargin < 3
    onset = 200; % assume onset = 200;
    threshold.earlypeak = 5;
    threshold.earlypeakwidth = 10;
    threshold.earlytime = 1:500 + onset;
    threshold.earlyprominence = 3;
    threshold.earlyzscore = 7;

    threshold.latepeak = 3;
    threshold.latepeakwidth = 50;
    threshold.lateprominence = 1.5;
    threshold.latezscore = 3;
    threshold.latetime = 200:800 + onset;
    threshold.art_timerange = [1:15, 600:2000] + onset;
end
threshold.art_timerange(threshold.art_timerange > size(CCEP_meanTr,2)) = [];

activated_redo = zeros(size(metaT,1),1);
act_time = cell(size(metaT,1),1);

parfor i = 1:size(metaT,1) % all data input takes about one minute
    warning('off','all')
    peaklocs = [];
    ccep = smoothdata(CCEP_meanTr(i, :),2, 'movmean', 10, 'omitnan');
    
    %% FIRST DETECT PEAK EARLY
    [~, locs] = findpeaks(ccep, 'MinPeakHeight', threshold.earlypeak,...
          'MinPeakProminence', threshold.earlyprominence, 'MinPeakWidth', threshold.earlypeakwidth);
     peaklocs = [peaklocs, locs];
     % TO FIP AND DETECT EARLY
     [~, locs] = findpeaks(-1.*ccep, 'MinPeakHeight', threshold.earlypeak,...
          'MinPeakProminence', threshold.earlyprominence, 'MinPeakWidth', threshold.earlypeakwidth);
     peaklocs = [peaklocs, locs];
     %{
     % USE HARD THRESHOLD
     earlyccep = zeros(size(ccep));
     earlyccep(1,threshold.earlytime) = ccep(1,threshold.earlytime);
     locs = find(abs(earlyccep) == max(abs(earlyccep),[],'omitnan') &...
         abs(earlyccep) >= threshold.earlyzscore);
    peaklocs = [peaklocs, locs];
     %}
     % PASS EARLY TIME
     peaklocs  = peaklocs(ismember(peaklocs,threshold.earlytime));

    %% DETECT LATE PEAKS
    peaklocs_late = []
    [~, locs] = findpeaks(ccep, 'MinPeakHeight', threshold.latepeak,...
          'MinPeakProminence', threshold.lateprominence, 'MinPeakWidth', threshold.latepeakwidth);
     peaklocs_late  = [peaklocs_late , locs];
    
     % TO FIP AND DETECT EARLY
     [~, locs] = findpeaks(-1.*ccep, 'MinPeakHeight', threshold.latepeak,...
          'MinPeakProminence', threshold.lateprominence, 'MinPeakWidth', threshold.latepeakwidth);
     peaklocs_late  = [peaklocs_late, locs];
     %{
     % USE HARD THRESHOLD
     lateccep = zeros(size(ccep));
     lateccep(1,threshold.latetime) = ccep(1,threshold.latetime);
     locs = find(abs(lateccep) == max(abs(lateccep),[],'omitnan') &...
         abs(lateccep) >= threshold.latezscore);
     peaklocs_late = [peaklocs_late, locs];
     %}
    % PASS LATE TIME
     peaklocs_late = peaklocs_late(ismember(peaklocs_late,threshold.latetime));
     % add up
     peaklocs = unique([peaklocs, peaklocs_late]);

     %% disregard artifacts
     peaklocs(ismember(peaklocs, threshold.art_timerange)) = [];
     peaktime = peaklocs-onset;
     peaktime(peaktime<=0) = [];
     act_time{i} = peaktime;
end

for i = 1:size(metaT,1)
    peaklocs =  act_time{i};
    if ~isempty(peaklocs)
        activated_redo(i) = 1;
    end
end