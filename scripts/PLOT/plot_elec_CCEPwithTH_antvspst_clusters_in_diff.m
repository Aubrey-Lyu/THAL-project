% this script will take the electrodes of one brain area that have distinct CECS disfferences
% with antTH and pstTH, and plot their locations in the brain surface
addpath('~/scripts/my_functions')
addpath('~/scripts/external')
addpath(genpath('~/scripts/Stanford/lbcn_preproc/vizualization'))
addpath(genpath('~/scripts/external/surfAnalysis-master'))
addpath('~/scripts/external/tight_subplot')

normBrainDir = '~/Dropbox/scripts/external/brainSurfer-main/brains';
wbDir = fullfile('/data/dian', 'workbench');
m = load(fullfile(wbDir,'affineMsurf2FSLR_group.mat'));

subjVar_dir = '~/Dropbox/Stanford_Matters/subjVars';

% load metaT
%metaT = readtable('~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore3_new2Sbjs/table_CCEPnewpipOutput_wholebrain_anatomical_info_activationRedone.csv');
result_folder = '~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore3_new2Sbjs/R_stats/intermediate';
plot_folder = '~/Dropbox/Stanford_Matters/data/THAL/CCEP/Plots/explore3/R_stats';
fn = 'TH_SM_cscwithclusters.csv';

t = readtable(fullfile(result_folder, fn));
elecID = t.elec_ID;

sbj_elec = regexp(elecID, ' ', 'split');
sbj_elec = vertcat(sbj_elec{:});

subjects = sbj_elec(:,1);
elec     = sbj_elec(:,2);

sblist = unique(subjects, 'stable');
cords_surf = [];

for is = 1:length(sblist)

    subject = sblist{is};
    skeys = strsplit(subject, '_');
    sbj_short =  strjoin(skeys(1:2), '_');

    fslabels = elec(strcmp(subjects, subject));
    % load subjVar
    clear subjVar
    load(fullfile(subjVar_dir, ['subjVar_', subject, '.mat']));
    lepto_cord = [];
    for ip = 1:length(fslabels)
        fskeys = strsplit(fslabels{ip},'-');
        el1 = fskeys{1};
        el2 = fskeys{2};
        el1_cord = subjVar.elinfo.LEPTO_coord(strcmpi(subjVar.elinfo.FS_label, el1),:);
        el2_cord = subjVar.elinfo.LEPTO_coord(strcmpi(subjVar.elinfo.FS_label, el2),:);
        el_cord = (el1_cord+el2_cord)/2;
        lepto_cord = [lepto_cord; el_cord];
    end
    % find correct affine matrix to use
    im = find(contains(m.sblist, sbj_short));
    Msurf2space = m.affineMsurf2FSLR_group{im};
    
    nativeBrainDir = fullfile(wbDir, sbj_short);
    % project to the right hemisphere
    [~,~,sfs_cords] = convertLepto2FS_LR(lepto_cord, Msurf2space,  nativeBrainDir, normBrainDir, 999);
    cords_surf = [cords_surf; sfs_cords];
end

%% separate cords based on clusters
c1_cords = cords_surf(t.cluster==1,:);
c2_cords = cords_surf(t.cluster==2,:);

fn =  fullfile(normBrainDir, 'S900.R.midthickness_MSMAll.32k_fs_LR.surf.gii');
this = gifti(fn);
cortex.vert = this.vertices;
cortex.tri = this.faces;

%% plot electrodes on the brain
color1 = [209 26 45]/266; % red
color2 = [32 161 98]/266; % purple

close; figure('Position', [ 670   579   843   360]);
[ha, pos] = tight_subplot(1,2);
sgtitle('within-SM 2 clusters')
% subplot 1
axes(ha(1)); 
ctmr_gauss_plot2(cortex,[0 0 0],[], 'right','medial');alpha(0.6)
hold on
scatter3(c1_cords(:,1),c1_cords(:,2),c1_cords(:,3), 'filled',...
    'MarkerFaceColor', color1, 'MarkerEdgeColor', color1, 'MarkerFaceAlpha',0.5)
scatter3(c2_cords(:,1),c2_cords(:,2),c2_cords(:,3), 'filled',...
    'MarkerFaceColor', color2, 'MarkerEdgeColor', color2, 'MarkerFaceAlpha',0.5)
hold off
% subplot 2
axes(ha(2)); 
ctmr_gauss_plot2(cortex,[0 0 0],[], 'right','lateral');alpha(0.6)
hold on
scatter3(c1_cords(:,1),c1_cords(:,2),c1_cords(:,3), 'filled',...
    'MarkerFaceColor', color1, 'MarkerEdgeColor', color1, 'MarkerFaceAlpha',0.5)
scatter3(c2_cords(:,1),c2_cords(:,2),c2_cords(:,3), 'filled',...
    'MarkerFaceColor', color2, 'MarkerEdgeColor', color2, 'MarkerFaceAlpha',0.5)
hold off
exportgraphics(gcf, fullfile(plot_folder, 'SM_CCEPwiththal_2clusters.png'))



