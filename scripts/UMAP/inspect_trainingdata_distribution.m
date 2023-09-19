data_folder = '~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore5_locked';
result_folder =  fullfile(data_folder, 'UMAP_learn');
metaT = readtable('~/Dropbox/Stanford_Matters/data/THAL/CCEP/results/explore5_locked/table_CCEPnewpipOutput_wholebrain_anatomical_info_activationRedone.csv');
%
sblist = unique(metaT.subject);
sblist(strcmp(sblist, 'S23_196_HL')) = [];
%%

mat = [];
for ss = 1:length(sblist)
    sb = sblist{ss};
    fname = sprintf('SpecReduceCollapse_%s', sb);
    load(fullfile(result_folder,'resample3', fname), 'filteridx_metaT', 'Vrpw', 'Vrph', 'Vrpc');
    mat = [Vrpw, Vrpc];
end