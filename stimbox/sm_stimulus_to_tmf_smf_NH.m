function [sprtmf, sprsmf, sumrtf, tmf, xmf] = sm_stimulus_to_tmf_smf_NH(stimulus, bin, taxis, faxis, nbin)
% sm_stim_modulation_params Peak temporal/spectral modulation value in stimulus
%
%    [sprtmf, sprsmf] = sm_stimulus_to_tmf_smf(stimulus, nlags, taxis, faxis)
%    
%    stimulus : Ntrials x Ndim stimulus matrix. Each row is one stimulus trial.
%    
%    bin (ms): is to calculate nlags (number of time bins in each trial. This should match the number
%       used in sta or mid calculations). Natsumi changed from nlags to bin 7Jun18
%   
%    taxis : temporal axis of stimulus or filter.
%    
%    faxis : spectral axis of stimulus or filter.
%
%    sprtmf : peak temporal modulation value during the stimulus trial frame.
%    sprsmf : peak spectral modulation value.
%    sumrtf : averaged rtf (added by Natsumi 7Jum18)

fprintf('%s\n', mfilename);
fs = length(taxis)/taxis(end);
nlags = floor(bin/1000 * fs);

sprtmf = zeros(1, size(stimulus,2));
sprsmf = zeros(1, size(stimulus,2));

maxtm = 40; %150
maxsm = 4;

% dt = taxis(2) - taxis(1);
% doct = log2(faxis(2)/faxis(1));
% 
% maxtm = 1 / dt / 2;
% maxsm = 1 / doct / 2;
% % maxtm = ceil( 1 / dt ); % changed to make this the same as strf_parameters Natsumi 05Jun18
% % maxsm = ceil( 1 / doct );

sumrtf = 0; %sumrtf = zeros(31,61);  %sumrtf = zeros(51,101);
countsumrtf = 0;
for i = nlags:size(stimulus,2)

    if ( mod(i, 10000)==0 )
        fprintf('%s: #%.0f of %.0f\n', mfilename, i, size(stimulus,2));
    end

    frame = stimulus(:,i-nlags+1:i);
    [tmf, xmf, rtf] = sm_mtf_sta2rtf(frame, taxis, faxis, maxtm, maxsm, nbin);
    %[rtf_fold, tmf_mtf, tmtf, smf_mtf, smtf] = sm_mtf_rtf2mtf(rtf, tmf, xmf);

    [indfmax, indtmax] = find(rtf == max(rtf(:)));
    xmf_max = max(xmf(indfmax));
    tmf_max = max(tmf(indtmax));


    %fprintf('tmf: %.1f, smf: %.2f\n', tmf_max, xmf_max);
    %pause

    % Uncomment the following to see the stimulus, the rtf, and
    % the estimated peak tmf and smf
%    clf(hf);
%    subplot(2,1,1);
%    imagesc(trialmat);
%
%    subplot(2,1,2);
%    imagesc(tmf, xmf, rtf);
%    title(sprintf('TBMF = %.2f, SBMF = %.2f', tmf_max, xmf_max));
%    pause(0.05); 

    sprtmf(i) = tmf_max;
    sprsmf(i) = xmf_max;
    sumrtf = sumrtf + rtf./sum(sum(rtf)); % sum up by averaging to 1, added by Natsumi 7Jun18
    countsumrtf = countsumrtf + 1;
end

sumrtf = sumrtf/countsumrtf;

return;


