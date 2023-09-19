function [rdata, rfreq, rtime] = reduceSpectrum(data,freq,time, kfreq, ktime, tozscore)
% input: data has the dimention of freq * time
% nfreq is a vector each element is the downsampled data kept for each
% frequency band, in the following order: delta, theta, alpha, beta, gamma,
% and high
% ntime is a vector with each element is the number of time points to
% downsample
% for the following time window : -5-10ms (artifact), 10-117ms (1st peak),
% 118-282ms (2nd peak), 283-700ms (3rd peak)
%(assuming the input time starts from -0.5 ms, sampling rate is 200)

% output: reduced spectrum matrix of a single CCEP with the dimension of (sum(kfreq) * sum(ktime));

%freq = power.freqs;
% (assuming the following frequency range)
%time = power.time*1000 (in ms);

if isempty(data)
    data = zeros(59,501);
end

if nargin < 4
    kfreq = [2, 2, 2, 2, 2, 1];
end

if nargin < 5
    ktime = [5 25 20 15];
end

if nargin < 6
    tozscore=1;
end
%%
freqSeg = [0.5 5; 5 8; 8 15; 15 30; 30 70; 70 256];
timeSeg = [-30 10; 10 117; 117 280; 280 800];

% reduce on the frequency domain
rdata0 = nan(sum(kfreq), length(time));
rfreq = nan(sum(kfreq),1);

for ik = 1:length(kfreq)
    k = kfreq(ik);

    if k==0
        continue;
    end
    
    fS_L = abs(freq-freqSeg(ik,1));
    fS_U = abs(freq-freqSeg(ik,2));
    iL = find(fS_L==min(fS_L));
    iU = find(fS_U==min(fS_U));
    d = data(iL:iU, :);
    iseg = segment(1:size(d,1),k);
    counter_kcum = sum(kfreq(1:ik-1));
    for ig = 1:size(iseg,1)
        rdata0(ig+counter_kcum,:) = mean(d(iseg(ig,1):iseg(ig,2), :), 1, 'omitnan');
        rfreq(ig+counter_kcum) = freq(round(0.5*(iseg(ig,1)+iseg(ig,2)))+iL-1);
    end
end

% reduce on the time domain
rdata = nan(sum(kfreq), sum(ktime));
rtime = nan(sum(ktime),1);

for ik = 1:length(ktime)
    k = ktime(ik);
    if k==0
        continue;
    end
    
    tS_L = abs(time-timeSeg(ik,1));
    tS_U = abs(time-timeSeg(ik,2));
    iL = find(tS_L==min(tS_L));
    iU = find(tS_U==min(tS_U));
    d = rdata0(:, iL:iU);
    iseg = segment(1:size(d,2),k);
    counter_kcum = sum(ktime(1:ik-1));
    for ig = 1:size(iseg,1)
        rdata(:,ig+counter_kcum) = mean(d(:, iseg(ig,1):iseg(ig,2)), 2, 'omitnan');
        rtime(ig+counter_kcum) = time(round(0.5*(iseg(ig,1)+iseg(ig,2)))+iL-1);
    end
end

%% zscore data
if tozscore == 1
    zrdata = zscore(rdata, 0, [1 2]);
    rdata = zrdata;
end