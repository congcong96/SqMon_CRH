function batch_stimulus_to_tmf_smf(stimfiles, nlags)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('nlags','var')
    nlags = 20;
end


for i = 1:length(stimfiles)
    
    fprintf('\nProcessing stimulus file %d of %d...\n', i, length(stimfiles))
    
    load(stimfiles{i}, 'stim_mat')
    if ~exist('stim_mat', 'var')
        load(stimfiles{i}, 'stim')
        stim_mat = stim;
    end
    if contains(stimfiles{i}, 'matrix')
        paramfile = regexprep(stimfiles{i}, 'matrix', 'param');
    else
        paramfile = regexprep(stimfiles{i}, 'stim', 'param');
    end
    load(paramfile, 'taxis', 'faxis', 'MaxFM', 'MaxRD');
    
    [sprtmf, sprsmf, tmfaxis, smfaxis] = stimulus_to_tmf_smf(stim_mat, nlags, taxis, faxis, MaxFM, MaxRD);
    
    outfile = regexprep(paramfile, 'param', 'mtf');
    
    save(outfile, 'sprtmf', 'sprsmf', 'tmfaxis', 'smfaxis');

end

