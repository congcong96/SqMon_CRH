function batch_stimulus_to_tmf_smf(stimfiles, nlags)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% batch process stimulus to get moment by moment tmf and smf
% INPUT:
%   stimfiles: cell containing address of stimuli files. 
%               The stimuli are stored as matlab matrix. 
%   nlags: number of time bins to get a snapshot of stimulus at each time 
%           point and calculate tmf and smf
% tmf and smf will be stored as matlab files at current directory with
% suffix -mtf.mat. varaibles in the output file:
%   sprtmf: moment by moment tmf
%   sprsmf: moment by moment tmf
%   tmfaxis: taxis in 2-D Fourier transform
%   smfaxis: faxis in 2-D Fourier transform

if ~exist('nlags','var')
    nlags = 20;
end

load(stimfiles, 'stim_mat')
if ~exist('stim_mat', 'var')
    load(stimfiles, 'stim')
    stim_mat = stim;
end
findfile = strfind(stimfiles, '\');
stimfile = stimfiles(findfile(end)+1:end);
if contains(stimfiles, 'matrix')
    paramfile = regexprep(stimfile, 'matrix', 'param');
else
    paramfile = regexprep(stimfile, 'stim', 'param');
end
paramfiles = [stimfiles(1:findfile(end)) paramfile];
load(paramfiles, 'taxis', 'faxis', 'MaxFM', 'MaxRD');

[sprphase, sprsmf, taxis, smfaxis] = stimulus_to_tmf_smf(stim_mat, nlags, taxis, faxis, MaxFM, MaxRD);
%[sprphase, sprtmf, saxis, tmfaxis] = stimulus_to_tmf_smf(stim_mat, nlags, taxis, faxis, MaxFM, MaxRD);

outfile = regexprep(paramfiles, 'param', 'mtf');

save(outfile,'sprphase', 'sprsmf', 'taxis', 'smfaxis');
%save(outfile,'sprphase', 'sprtmf', 'saxis', 'tmfaxis');

% for i = 1:length(stimfiles)
%     
%     fprintf('\nProcessing stimulus file %d of %d...\n', i, length(stimfiles))
%     
%     load(stimfiles{i}, 'stim_mat')
%     if ~exist('stim_mat', 'var')
%         load(stimfiles{i}, 'stim')
%         stim_mat = stim;
%     end
%     if contains(stimfiles{i}, 'matrix')
%         paramfile = regexprep(stimfiles{i}, 'matrix', 'param');
%     else
%         paramfile = regexprep(stimfiles{i}, 'stim', 'param');
%     end
%     load(paramfile, 'taxis', 'faxis', 'MaxFM', 'MaxRD');
%     
%     [sprtmf, sprsmf, tmfaxis, smfaxis] = stimulus_to_tmf_smf(stim_mat, nlags, taxis, faxis, MaxFM, MaxRD);
%     
%     outfile = regexprep(paramfile, 'param', 'mtf');
%     
%     save(outfile, 'sprtmf', 'sprsmf', 'tmfaxis', 'smfaxis');
% 
% end

