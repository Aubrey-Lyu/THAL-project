function [VRPW, VRPC] = get_vectCCEP(filteridx_in_metaT, inputFolder)

if exist('metaT','var')==0
    error('metaT is not loaded.')
end

if nargin < 2
    inputFolder = '~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore5_locked/UMAP_learn/resample3';
end

fIdx = filteridx_in_metaT;

sblist = metaT.subect(fIdx);

fn = dir(fullfile(inputFolder, ['*' sblist{1} '.mat']));
eg = load(fullfile(data_folder, fn.name));

VRPW = nan(length(fIdx), size(eg.Vrpw,2));
VRPC  = nan(length(fIdx), size(eg.Vrpc,2));
%filterIdx = [];% the order may be different than the input filteridx_in_metaT

    tic
    for sn = 1:length(sblist)
        sbj = sblist{sn};
        fn = dir(fullfile(inputFolder, ['*' sbj '.mat']));
        
        sIdx = find(strcmp(metaT.subject, sblist{sn}));
        d = load(fullfile(data_folder, fn.name));

        [filter_sidx,~, row2select] = intersect(fIdx(sIdx), d.filteridx_metaT,  'stable'); % +1 account for python emuerating starting from 0

        spect1 = Vrpw(row2select,:,:);
        spect2 = Vrpc(row2select,:,:);

        VRPW(fIdx(sIdx),:) = spect1;
        VRPC(fIdx(sIdx),:) = spect2;

        %filterIdx = [filterIdx; filter_sidx];

    end
    
end