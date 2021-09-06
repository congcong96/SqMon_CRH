function display_dmr_stm(SpecFile, nbin)
% from Figure4Paper_SoundStimuli05 NH
% --edited by Congcong 11/05/2019
bin = 125; % ms
fileName = SpecFile;
savefilename = sprintf('%s-spr-rtf.mat',fileName(1:end-4));
if isempty(dir(savefilename))
    [stimulus,taxis,faxis] = read_entire_spr_file_2(fileName);
    [sprtmf, sprsmf, sumrtf, tmf, smf] = sm_stimulus_to_tmf_smf_NH(stimulus, bin, taxis, faxis, nbin);
    save(savefilename,'sprtmf', 'sprsmf', 'sumrtf', 'tmf', 'smf');
    
    figure
    imagesc(tmf,smf,sumrtf)
    axis xy;
end