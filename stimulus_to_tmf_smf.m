function [sprtmf, sprsmf, tmfaxis, smfaxis] = stimulus_to_tmf_smf(stimulus, nlags, taxis, faxis, maxtm, maxsm)
% sm_stim_modulation_params Peak temporal/spectral modulation value in stimulus
%
%    [sprtmf, sprsmf] = sm_stimulus_to_tmf_smf(stimulus, nlags, taxis, faxis)
%    
%    stimulus : Ntrials x Ndim stimulus matrix. Each row is one stimulus trial.
%       Saved in: dmr-50flo-40000fhi-4SM-150TM-40db-96kHz-96DF-30min_DFt5_DFf8-matrix.mat
%    
%    nlags : number of time bins in each trial. This should match the number
%       used in sta or mid calculations.
%       
%    taxis : temporal axis of stimulus or filter.
%       Saved in: dmr-50flo-40000fhi-4SM-150TM-40db-96kHz-96DF-30min_DFt5_DFf8_param.mat
%    
%    faxis : spectral axis of stimulus or filter.
%       Saved in: dmr-50flo-40000fhi-4SM-150TM-40db-96kHz-96DF-30min_DFt5_DFf8_param.mat
%
%    sprtmf : peak temporal modulation value during the stimulus trial frame.
%    sprsmf : peak spectral modulation value.
%
%    sprtmf and sprsmf are saved in:
%       dmr-50flo-40000fhi-4SM-150TM-40db-96kHz-96DF-30min_DFt5_DFf8-mtf-hires.mat 
%    
fprintf('%s\n', mfilename);

sprtmf = zeros(1, size(stimulus,2));
sprsmf = zeros(1, size(stimulus,2));

% maxtm = 40;
% maxsm = 4;

% dt = taxis(2) - taxis(1);
% doct = log2(faxis(2)/faxis(1));
% 
% maxtm = 1 / dt / 2;
% maxsm = 1 / doct / 2;

for i = nlags:size(stimulus,2)

    if ( mod(i, 10000)==0 )
        fprintf('%s: #%.0f of %.0f\n', mfilename, i, size(stimulus,2));
    end

    frame = stimulus(:,i-nlags+1:i);
    
%     if mod(i, 100) == 0
%         [tmf, xmf, rtf] = sm_mtf_sta2rtf(frame, taxis, faxis, maxtm, maxsm, 'y');
%     else
        [tmfaxis, smfaxis, rtf] = sm_mtf_sta2rtf(frame, taxis, faxis, maxtm, maxsm);
%     end
    
%     [rtf_fold, tmf_mtf, tmtf, smf_mtf, smtf] = sm_mtf_rtf2mtf(rtf, tmf, xmf);

    [indfmax, indtmax] = find(rtf == max(rtf(:)));
    xmf_max = max(smfaxis(indfmax));
    tmf_max = max(tmfaxis(indtmax));


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
end

return;


