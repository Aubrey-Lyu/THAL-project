
% This scripts will project freesurfer surface space to the FS_LR space
% The FS_LR has symetrical left and right hemisphere.

toolboxDir = '~/Dropbox/scripts/external/surfAnalysis-master';
addpath(genpath(toolboxDir))
addpath(genpath('~/Dropbox/scripts/external/matlab_GIfTI-master'))
addpath('~/Dropbox/scripts/external')
addpath(genpath('/Users/dl577/Dropbox/scripts/external/brainSurfer-main'))

fsDir = getenv('SUBJECTS_DIR');
brainDir = '~/Dropbox/scripts/external/brainSurfer-main/brains';

outputdir = fullfile('/Volumes/Lacie', 'workbench');
if ~exist(outputdir, 'dir'); mkdir(outputdir);end

sblist = dir(fullfile(fsDir, 'S2*'));
sblist = {sblist.name};

affineMsurf2FSLR_group = cell(length(sblist),1);
for ss = 1:length(sblist)
    subj_name = sblist{ss};
    disp('-------------------------------')
    fprintf('Processing %s...\n', subj_name)
    % Resampling individual surfaces into fs_LR space
    Msurf2space = surf_resliceFS2WB(subj_name,fsDir,toolboxDir, outputdir ,'resolution','32k') ;
    affineMsurf2FSLR_group{ss} = Msurf2space;
end
save(fullfile(outputdir, 'affineMsurf2FSLR_group.mat'), 'affineMsurf2FSLR_group','sblist')

%%
%{
t1 = fullfile(fsDir, subj_name, 'elec_recon/T1.nii.gz');
outdir = fullfile(fsDir, subj_name, 'elec_recon/foci');
if ~exist(outdir, 'dir'); mkdir(outdir);end

% create point
bashcmd = sprintf('fslmaths %s -mul 0 -add 1 -roi %d 1 %d 1 %d 1 0 1 %s -odt float',...
    t1, coord(1), coord(2), coord(3), fullfile(outdir,[fslabel, '_point']) );
disp(bashcmd); system(bashcmd,'-echo');
% create foci (will take some time)
radius = 5;
bashcmd = sprintf('fslmaths %s -kernel sphere %d -fmean %s -odt float',...
    fullfile(outdir,[fslabel, '_point']), radius, fullfile(outdir,[fslabel, '_sphere']) );
disp(bashcmd); system(bashcmd,'-echo');
% create binary mask
bashcmd = sprintf('fslmaths %s -thr 0.001 -bin %s ',...
    fullfile(outdir,[fslabel, '_sphere']), fullfile(outdir,[fslabel, '_sphereMask']));
disp(bashcmd); system(bashcmd,'-echo');
%}

%% project LEPTO freesurfer space to FS_LR space
ss = 24;
subj_name = sblist{ss};

fn = dir(sprintf('~/Dropbox/subjVars/subjVar_%s*_jplabel.mat', subj_name));
clear subjVar
load(sprintf('~/Dropbox/subjVars/%s', fn(1).name));
lepto_cords = subjVar.elinfo.LEPTO_coord;
fscords = lepto_cords;

% use affine transform
Msurf2space = affineMsurf2FSLR_group{ss};
[fscords(:,1),fscords(:,2),fscords(:,3)]=surf_affine_transform(lepto_cords(:,1),lepto_cords(:,2),lepto_cords(:,3),Msurf2space);

% [cortex1.vert, cortex1.tri] = read_surf(fullfile(fsDir,  subj_name, 'surf', 'lh.pial.T1'));
% [cortex2.vert, cortex2.tri] = read_surf(fullfile(fsDir,  subj_name, 'surf', 'rh.pial.T1'));
% if find(cortex1.tri==0); cortex1.tri = cortex1.tri+1;end
% if find(cortex2.tri==0); cortex2.tri = cortex2.tri+1;end
% % close;
% % ctmr_gauss_plot2(cortex1,[0 0 0],[],'left','medial'); hold on;
% % ctmr_gauss_plot2(cortex2,[0 0 0],[],'right','medial')
% % alpha(0.75);
% % scatter3(lepto_cords(:,1),lepto_cords(:,2),lepto_cords(:,3),...
% %     'MarkerFaceColor', 'red')
% % hold off

% FS_LR space - individual space
fn1 = fullfile(outputdir, subj_name, [subj_name, '.L.pial.T1.32k.surf.gii']);
this = gifti(fn1);
cortexFS1.vert = this.vertices;
cortexFS1.tri = this.faces;
fn2 = fullfile(outputdir, subj_name, [subj_name, '.R.pial.T1.32k.surf.gii']);
this = gifti(fn2);
cortexFS2.vert = this.vertices;
cortexFS2.tri = this.faces;
% figure;
% ctmr_gauss_plot2(cortexFS1,[0 0 0],[],'left','medial'); hold on;
% ctmr_gauss_plot2(cortexFS2,[0 0 0],[],'right','medial');
% alpha(0.9);
% scatter3(fscords(:,1),fscords(:,2),fscords(:,3),...
%     'MarkerFaceColor', 'blue')
% FS_LR standard space
fn = fullfile(brainDir, 'S900.L.midthickness_MSMAll.32k_fs_LR.surf.gii');
this = gifti(fn);
cortexFS0_1.vert = this.vertices;
cortexFS0_1.tri = this.faces;


fn =  fullfile(brainDir, 'S900.R.midthickness_MSMAll.32k_fs_LR.surf.gii');
this = gifti(fn);
cortexFS0_2.vert = this.vertices;
cortexFS0_2.tri = this.faces;

% %% use MNI space
% mni_cords = subjVar.elinfo.MNI_coord;
%  [cmcortex1.vert, cmcortex1.tri] = read_surf(fullfile(fsDir, 'fsaverage' , 'surf', 'lh.pial'));
%  [cmcortex2.vert, cmcortex2.tri] = read_surf(fullfile(fsDir, 'fsaverage' , 'surf', 'rh.pial'));
% if find(cmcortex1.tri==0); cmcortex1.tri = cmcortex1.tri+1;end
% if find(cmcortex2.tri==0); cmcortex2.tri = cmcortex2.tri+1;end
% figure;
% ctmr_gauss_plot2(cmcortex1,[0 0 0],[],'left','medial'); hold on;
% ctmr_gauss_plot2(cmcortex2,[0 0 0],[],'right','medial')
% alpha(0.75);
% scatter3(mni_cords(:,1),mni_cords(:,2),mni_cords(:,3),...
%     'MarkerFaceColor', 'green')
% hold off

%% find the index of the nearest brain surface (radius = 5)

vert = [cortexFS1.vert; cortexFS2.vert];
vIofElec = nan(size(fscords,1),1);

for i = 1:size(fscords,1)
    e = fscords(i,:); % elec
    ee = repmat(e, size(vert,1),1);
    ed = sqrt(sum((e-vert).^2, 2));
    I = find(ed == min(ed));
    if ed(I) < 5
        vIofElec(i) = I;
    end
end
vIofElec_L = vIofElec(vIofElec<=size(vert,1)/2);
vIofElec_R = vIofElec(vIofElec>size(vert,1)/2);
vIofElec_R = vIofElec_R - size(vert,1)/2;
% find the corresponding coordinates in the standard space
vert0 = [cortexFS0_1.vert; cortexFS0_2.vert];
sfs_cords = nan(size(fscords));
sfs_cords(~isnan(vIofElec),:) = vert0(vIofElec(~isnan(vIofElec)),:);


% Tc = table( );
% Tc.FScord = sfs_cords;
% Tc.JPlabel = subjVar.elinfo.JP_label;
% 
%% plot
figure;
ctmr_gauss_plot2(cortexFS0_1,[0 0 0],[],'right','medial');
alpha(0.7);hold on
ctmr_gauss_plot2(cortexFS0_2,[0 0 0],[],'right','medial');


scatter3(sfs_cords(:,1),sfs_cords(:,2),sfs_cords(:,3),...
   55, 'MarkerFaceColor', 'blue')
text(sfs_cords(:,1),sfs_cords(:,2),sfs_cords(:,3),...
    Tc.JPlabel)
 hold off
%  
%  %% plot MNI coord
%  figure;
% ctmr_gauss_plot2(cortexFS0_1,[0 0 0],[],'right','medial');
% alpha(0.7);hold on
% ctmr_gauss_plot2(cortexFS0_2,[0 0 0],[],'right','medial');
% 
% 
% scatter3(mni_cords(:,1),mni_cords(:,2),mni_cords(:,3),...
%    55, 'MarkerFaceColor', 'red')
% text(mni_cords(:,1),mni_cords(:,2),mni_cords(:,3),...
%     Tc.JPlabel)
%  hold off

