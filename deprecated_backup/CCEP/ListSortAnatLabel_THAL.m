function newlist = ListSortAnatLabel(oldlist, level)
% sort inconsistent JP label, input and output are cell arrays
% level is the level of anatomical details.
% level 1: anat label
% level 2: anat label with subdivisions

if nargin < 2
    level = 2;
end

newlist = strtrim(oldlist);
newlist = revalueCell(newlist, strcmp(newlist, 'MD'), 'THALAMUS MID');

newlist = revalueCell(newlist, strcmp(newlist, 'VA-TC'), 'OFC');
newlist = revalueCell(newlist, strcmp(newlist, 'THALAMUS MD'), 'THALAMUS MID');
newlist = revalueCell(newlist, strcmp(newlist, 'AC'), 'ACC');
newlist = revalueCell(newlist, strcmp(newlist, 'VM-PFC'), 'OFC');
newlist = revalueCell(newlist, strcmp(newlist, 'MPFC'), 'MFC');
newlist = revalueCell(newlist, strcmp(newlist, 'HIPP-ANT'), 'HPC ANT');
newlist = revalueCell(newlist, strcmp(newlist, 'AMYGDLA'), 'AMY');
newlist = revalueCell(newlist, strcmp(newlist, 'FUSIFORM GYRUS FG'), 'FG');
newlist = revalueCell(newlist, strcmp(newlist, 'PRECENTRAL GYRUS'), 'preCG');
newlist = revalueCell(newlist, strcmp(newlist, 'POSTCENTRAL GYRUS'), 'postCG');
newlist = revalueCell(newlist, strcmp(newlist, 'INSULA'), 'INS');
newlist = revalueCell(newlist, strcmp(newlist, 'SMG'), 'IPL');

newlist = revalueCell(newlist, strcmp(newlist, 'CLAUSTRUM'), 'CLT');
newlist = revalueCell(newlist, strcmp(newlist, 'THALAMUS ANT'), 'antTH');
newlist = revalueCell(newlist, strcmp(newlist, 'THALAMUS POST'), 'pstTH');
newlist = revalueCell(newlist, strcmp(newlist, 'THALAMUS MID'), 'midTH');
newlist = revalueCell(newlist, strcmp(newlist, 'THALMUS POST'), 'pstTH');
newlist = revalueCell(newlist, strcmp(newlist, 'THALMUS MID'), 'midTH');
newlist = revalueCell(newlist, strcmp(newlist, 'THALMUS ANT'), 'antTH');
newlist = revalueCell(newlist, strcmp(newlist, 'ITL'), 'IPL');
newlist = revalueCell(newlist, strcmp(newlist, 'IFG'), 'LFC');
newlist = revalueCell(newlist, strcmp(newlist, 'MPFC'), 'MFC');
newlist = revalueCell(newlist, strcmp(newlist, 'SMG'), 'IPL');
 
if level == 1
    %==============================================
    %group anatomical labelling
    newlist = revalueCell(newlist, ismember(newlist, {'HPC ANT', 'HPC MID', 'HPC POST'}), 'HPC');
    newlist = revalueCell(newlist, ismember(newlist, {'postCG', 'preCG', 'SOMATOMOTOR', 'SOMATOSENSORY'}), 'SM');
%     newlist = revalueCell(newlist, ismember(newlist, {'antTH', 'midTH', 'pstTH'}), 'TH');
end

end