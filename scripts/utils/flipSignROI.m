function [ccep_flipped,toflip,C] = flipSignROI(ccep, peakTimeMat, onset, corr_thr)
% ccep dimension = stimPair * T
% assuming samping rate is 1000
% peakTimeMat = data array with each row indicating the time point when an
% T = metaT(filterIdx1,:);
% vars = T.Properties.VariableNames;
% ttpIdx = find(contains(vars, 'pks_time_'));
% peakTimeMat = table2array(T(:,ttpIdx));
% activation is detected for one pair, baseline = 200 points
% if ccep data is baseline-trimmed, then the onset = 0;
% corr_thr: the threshold of the mean correlation value to determine to flip
% or not
% ccep_flipped: flipped CCEP to insure the maximum correlation with a
% template;
% C records the mean correlation value along the iteration
if size(ccep,1) == size(peakTimeMat,1)
    ccep = ccep';
else
    if size(ccep,2) ~= size(peakTimeMat,1)
        error('Data dimension and peakTimeMat dimension do not match.')
    end
    % otherwise, the ccep has the desired dimention, do nothing.
end


if nargin < 3
    % assuming the onset starts from the 200 time point
    onset = 200;
end

if nargin < 4
    corr_thr = 0.2;
end
corr_thr = abs(corr_thr);
% time used to do calculation
calc_time = (1+onset):size(ccep,1) ;


%% initial flip
actIdx = ~isnan(mean(peakTimeMat, 2, 'omitnan'));
ccep_act = ccep(calc_time,actIdx);
peaks = peakTimeMat(actIdx,:);
hottimezone = [16:80 200:500];

for i = 1:size(ccep_act,2)
    peak = peaks(i,:);
    peak(~isnan(peak) | peak>max(calc_time)) = [];
    if isempty(peak)
        val = max(ccep_act(hottimezone,i), [], 1, 'omitnan');
    else
        pk = [peak(1) peak(ismember(peak,hottimezone))]; % assumption: the first (16-80 ms) and third (200-500) peaks are positive; consider time before 15ms is artifact
        if isempty(pk) || isnan(pk)
            val = max(ccep_act(hottimezone,i), [], 1, 'omitnan');
        else
            val = max(ccep_act(pk,i), [], 1, 'omitnan'); % find the biggest value in the range
        end
    end
    if val < 0
        ccep(:,i) = -1.*ccep(:,i);
    end
end

%% iterations
% create template
tokeep = 1:size(ccep,2); % initial: to use all

Nitr = size(ccep,2);
C = nan(Nitr,1);

for c = 1:Nitr

    template = mean(ccep(calc_time,tokeep), 2, 'omitnan');

    corr = corrcoef([template, ccep(calc_time,:)]);
    cor = corr(2:end,1);
    avcor0 = mean(cor, 'omitnan');

    toflip = cor < (-1*corr_thr); % negatively correlated
    tokeep = cor > corr_thr;

    % to do flip tentatively
    ccep1 = ccep; toflip1 = toflip;
    ccep1(:,toflip1) = -1.* ccep1(:,toflip1);
    corr = corrcoef([template, ccep1(calc_time,:)]);
    cor = corr(2:end,1);
    avcor1 = mean(cor, 'omitnan');

    %-- evaluate ---
    improvement = avcor1-avcor0;
    C(c) = improvement;

    if improvement <= 0
        break
    else
        ccep = ccep1;
        toflip = toflip1;
    end

end

ccep_flipped = ccep';
C(isnan(C)) = [];

%
%
%
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [E,S]=kt_pca(X)
% % function [E,S]=kt_pca(X)
% % This is an implementation of the linear kernel PCA method ("kernel trick")
% % described in "Kernel PCA Pattern Reconstruction via Approximate Pre-Images"
% % by Schölkopf et al, ICANN, 1998, pp 147-15
% % See also course lectures by K.R. Muller (TU Berlin)
% %
% % Inputs:
% % X(T,N) - Matrix of data in. Only need this trick if T>>N
% %
% % Outputs:
% % E(T,N) - Columns of E are estimated Eigenvectors of X, in descending order
% % S(1,N) - Corresponding estimated Eigenvalues of X, in descending order
% %
% % kjm, june 2020
%
% %% use the "kernel trick" to estimate eigenvectors of this cluster of pair groups
%     [F,S2]=eig(X.'*X); % eigenvector decomposition of (covariance of transpose) - may want to replace this with cov function so mean is subtracted off
%     [S2,v_inds]=sort(sort(sum(S2)),'descend'); F=F(:,v_inds); %reshape properly
%
%     S=S2.^.5; % estimated eigenvalues of both X.'*X and X*X.'
%
%     %     ES=X*F.';
%     ES=X*F; % kernel trick
%     E=ES./(ones(size(X,1),1)*S); % divide through to obtain unit-normalized eigenvectors