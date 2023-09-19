clear
home_dir = '/data/dian';
%% plot CCEP time domain
ccep_folder = fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore5_locked');
sepc_folder = fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/UMAP/ALLDATA_semisupervise');
outfigDir = fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/PLOTS/UMAP/supervisedLearn_Group');

%% load cleaned ccep data
metaT = readtable(fullfile(ccep_folder, 'table_CCEPnewpipOutput_wholebrain_anatomical_info_activationRedone2.csv'));
load(fullfile(sepc_folder, 'vectCCEP_featured.mat'));
load(fullfile(ccep_folder, 'CCEP_all_flat_meanTr_cleaned.mat'));%, 'CCEP_flat' each row corresponding to the metaT
zccep_clean(badChan,:) = nan;

%% plot
tic
close all
figure('Position', [990   503   755   796]);

for i = 1:4%length(vectCCEP_featured)
    subplot(2,2,i)
    feature = vectCCEP_featured(i).feature;
    r = randi([1 256],1, size(vectCCEP_featured(i).idx_in_metaT, 1));
    g = randi([1 256],1, size(vectCCEP_featured(i).idx_in_metaT, 1));
    b = randi([1 256],1, size(vectCCEP_featured(i).idx_in_metaT, 1));

    ccep_featured = zccep_clean(vectCCEP_featured(i).idx_in_metaT,1:1200);
    % flip ccep, takes ~ half a min to run 
    ccep_featured = flipSignROI(ccep_featured', [], [], [], 200:500); %hottimezone = 200:500;
   
    meanccep = mean(ccep_featured, 1, 'omitnan');
    
    baselineSTD = std(ccep_featured(:, 1:180), [] ,[1 2], 'omitnan');
    zccep_ = meanccep/baselineSTD;
    
    hold on
    for ir = 1:size(ccep_featured,1)
        ccep1 = ccep_featured(ir,100:end)';
        col = [r(ir),g(ir),b(ir)]/256;
        plot(ccep1, 'Color', [col, 0.05])
    end
  
    plot(zccep_(:,100:end)', 'Color', 'w', 'LineWidth',0.5)
    hold off
    title(feature)
    ylim([-40, 60])
    set(gca,'Color','k')
    ax = gca; % Get handle to current axes.
    ax.XColor = 'w'; 
    ax.YColor = 'w';
end
set(gcf, 'color', 'k'); 
box off
sgtitle('Brain stimulation evoked potential (z-scored)')
exportgraphics(gcf, fullfile(outfigDir, 'prelim_ccep_flipped_trials_colored.pdf') ,'Resolution',300)
toc


