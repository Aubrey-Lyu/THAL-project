function subplotCCEPspectrum(d, ttl,time,freq,cmap,labelstring, nstep_x, nstep_y)

if nargin < 7
nstep_x = 30;
end

if nargin < 8
nstep_y = 10;
end

imagesc(d); xlabel('time (s)'); ylabel('Frequency')
title(ttl)
xticklabels = time;
xticks = 1:size(d, 2);

yticklabels = freq;
yticks = 1:size(d, 1);

set(gca, ...
    'XTick', xticks(1:nstep_x:end), 'XTickLabel', xticklabels(1:nstep_x:end),...
    'YTick', yticks(1:nstep_y:end), 'YTickLabel', yticklabels(1:nstep_y:end),...
    'YDir','normal')

colormap(gca, cmap);
a = colorbar;
a.Label.String = labelstring;

freqBands = {'delta', 'theta', 'alpha', 'beta', 'gamma', 'high'};
fbs = [0.5 4 8 15 30 70];
freq_num = cellfun(@str2num, freq);

for f = 1:length(freqBands)
    band = freqBands{f};
    freq_diff = abs(freq_num - fbs(f));
    idx  = find(freq_diff == min(freq_diff));
    if idx <= size(d,1)
    yline(idx,'--', band,'FontSize',8);
    end
end
 xline( find(abs(time)==(min(abs(time)))), ':', 'LineWidth',2)