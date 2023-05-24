home_dir = '/data/dian';
metaT = readtable(fullfile(home_dir, 'Dropbox/Stanford_Matters/data/THAL/COHORT/table_CCEPnewpipOutput_wholebrain_anatomical_info.csv'));
JP_label_out = metaT.JP_label_out1;
JP_label_in = metaT.JP_label_in1;
jplabels = [JP_label_out; JP_label_in];
chanO = strcat(metaT.subject,', ', metaT.sc1);
chanI = strcat(metaT.subject,', ', metaT.rc1);
chan = [chanO; chanI];
jplabel_u = unique([JP_label_out; JP_label_in])
%% organize JP_label
% check unexpected JP labels
unique(chan(strcmp(jplabels, 'INSULA/LFC'))) % only one electrode, could exclude

