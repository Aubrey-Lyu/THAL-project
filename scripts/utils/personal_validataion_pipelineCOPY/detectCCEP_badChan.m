function [badChan, zccepTr_cleaned] = detectCCEP_badChan(zccepTr, do_linearScrub, time_onset,fs)
% the function detect CCEP badChan based on the trial average signl
% zccepTr: [chan * timepoints] is trial-averaged zccep at an electrode
% level
% badChan return badChan index corresponding to CCEP data
% zccep_cleaned has the same dimension as CCEP
if nargin < 2 || isempty(do_linearScrub) % optional
    do_linearScrub = 1;
end

if nargin < 3
    time_onset = 200; %assuming the stim happens at the 200 ms
end

if nargin < 4
    fs = 1000; %assuming the samping rate is 1000
end

badChan = [];
ccepTr_baseline = zccepTr(:, 1:(round(time_onset.*fs/1000)-1));
ccepTr_after = zccepTr(:, round((time_onset).*fs/1000):end);

%% to do simple linear trend scrub
linear_trend = (1: size(ccepTr_after,2))';
if do_linearScrub == 1
    parfor i = 1:size(ccepTr_after,1)
        ccep1 = ccepTr_after(i,:);
        gl = fitlm(linear_trend, ccep1);
        yhat = predict(gl, linear_trend);
        res = ccep1 - yhat';
        ccepTr_after(i,:) = res;
    end
end
%% detect big spikes in the second derivative > 20ms after onset
ccepTr_late = ccepTr_after(:, round((20).*fs/1000):end);
ccepTr_diff = abs(normalize(diff(diff(ccepTr_late')),2, "zscore")); % note that diff will take derivative of column vectors
%  decide the big spike is when the second derivative is bigger than 10
[~,badChan_1] = find(ccepTr_diff>10);
badChan = [badChan; badChan_1];

%% filter out channels with zscore >= 7 after time > 1000ms
ccepTr_late = ccepTr_after(:, round((1000).*fs/1000):end);
[badChan_2, ~] = find(abs(ccepTr_late)>=7);
badChan = [badChan; badChan_2];
%{
% plot to check
badChan_1=unique(badChan_1);
badChan_2=unique(badChan_2);
badChan_ = unique(badChan_2(~ismember(badChan_2, badChan_1)));
close all
figure; plot(zccep(unique(badChan_),:)'); legend(CCEP.recordingchannel(unique(badChan_))); title('raw ccep')
figure; plot(ccep_after(unique(badChan_),:)'); legend(CCEP.recordingchannel(unique(badChan_))); title ('detrended ccep')
%}
badChan = unique(badChan);

%% clean ccep to check the effect
zccepTr_cleaned = [ccepTr_baseline ccepTr_after];
zccepTr_cleaned(badChan,:) = nan;
%{
figure;
subplot(1,2,1);plot(zccep'); title('raw zccep');
ylims = [min(zccepTr,[],'all'), max(zccepTr,[],'all')];
ylims = [-100 100];
ylim(ylims)
subplot(1,2,2); plot(zccepTr_cleaned'); title('cleaned zccep')
ylim(ylims)
%}