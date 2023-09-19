function [phase_flipped, toflip] = flipPhase(phaseDat, zccep_Tr, peakTimeMat, toflip)
% phase if flipped based on the maximally uniform waveforms of the
% averaged CCEP
if ~exist('toflip', 'var')||isempty(toflip)
[~, toflip,~]   = flipSignROI(zccep_Tr, peakTimeMat);
end
n = size(phaseDat,1);
ntime = size(phaseDat,2);
nfreq = size(phaseDat,3);

flipsign = ones(1,n);
flipsign(toflip) = -1;
flipsign_mat = permute(reshape(repmat(flipsign, ntime*nfreq,1),ntime, nfreq,[]), [3 1 2]);
phase_flipped = flipsign_mat.*phaseDat;


%{
if nargin<2
seed=66;
end
rng(seed);

n = size(phase,1);
ntime = size(phase,2);
nfreq = size(phase,3);

phase_vec = reshape(phase,size(phase,1),[],1)';

toflip = 2*(randi([0 1],n,1))'-1;
toflip_mat = permute(reshape(repmat(toflip, ntime*nfreq,1),ntime, nfreq,[]), [3 1 2]);
phase_flip = toflip_mat.*phase;
fase_mean = squeeze(mean(phase_flip,1,'omitnan'));
fase_mean_vec = reshape(fase_mean,[],1);
allph = [fase_mean_vec  phase_vec];

R = corrcoef(allph);
r = mean(R(2:end,1), 'omitnan');
%}