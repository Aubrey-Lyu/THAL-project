dir_base = '/data/dian/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore1_dataSort_Denoise/BPC';
plotDir = '/data/dian/Dropbox/Stanford_Matters/data/THAL/CCEP/Plots/explore1';
%ty = 'COR-COR';


types = {'THAL-COR', 'COR-THAL', 'THAL-THAL'};
close all
for t = 1:length(types)
    ty = types{t};
    fn = dir(fullfile(dir_base, ty, 'BPC_output*.mat'));
fns = {fn.name}';
cd(fullfile(dir_base, ty))
%plot
figure('Position', [675        -585        2033        1536]);
nCol = 4; nSubplot = length(fns);
for ss = 1:length(fns)
    keys = strsplit(fns{ss}, '_');
    tl = strjoin(keys(3:end), '_');
    tl = tl(1:end-4);
    d = load(fns{ss});
    nRow = ceil(nSubplot/nCol);
    subplot(nRow, nCol, ss)
    dat = [d.B.curve];
    plot(dat(1:1200,:)); title(strrep(tl,'_','\_'))
    hold on
    xline(200, 'k-')
    
end

saveas(gcf, fullfile(plotDir, ['BPC_allsbjs_' ty '.png']))
end