function wave_interp = interp_stimArtifact(wave, artifact_search, fs)
if nargin < 3
    fs = 1000;
end

 wave(isinf(wave)) = NaN;
    if sum(isnan(wave)) == length(wave)
        wave(:) = 0;
    end
    %define channel threshold difference
    upper_lim = prctile(diff(wave),99);
    lower_lim = prctile(diff(wave),1);
    pk_th = min(abs([upper_lim,lower_lim]));
    [pks,locs,width] = findpeaks(diff(wave),'MinPeakHeight',pk_th,'MinPeakDistance',fs*1);
    [pks_neg,locs_neg,width] = findpeaks(-1.*diff(wave),'MinPeakHeight',pk_th,'MinPeakDistance',fs*1.5);
    locs = cat(1,locs,locs_neg);
    if isempty(locs)
        return
    end
    pk_win = [-4,4];
    artifact_found = [];
    for l = 1:length(locs)
        artifact_found = cat(1,artifact_found,...
            ((locs(l) + pk_win(1)):(locs(l) + pk_win(2)))');
    end
    artifact_ind = intersect(artifact_search,artifact_found);

    %% High thresh artifact
    %define channel threshold difference
    upper_lim = prctile(diff(wave),99.9);
    lower_lim = prctile(diff(wave),.01);
    pk_th = min(abs([upper_lim,lower_lim]));
    [pks,locs,width] = findpeaks(diff(wave),'MinPeakHeight',pk_th);
    [pks_neg,locs_neg,width] = findpeaks(-1.*diff(wave),'MinPeakHeight',pk_th,'MinPeakDistance',fs*1.5);
    locs = cat(1,locs,locs_neg);
    if ~isempty(locs)
        
        artifact_found = [];
        for l = 1:length(locs)
            artifact_found = cat(1,artifact_found,...
                ((locs(l) + pk_win(1)):(locs(l) + pk_win(2)))');
        end
        artifact_ind = union(artifact_ind,artifact_found);
    end
    artifact_ind(artifact_ind<1 | artifact_ind>length(wave)) = [];
    wave(artifact_ind) = NaN;
    wave_interp = fillgaps(wave,200,25);