function [sprphase, sprsmf, taxis, smfaxis] = stimulus_to_tmf_smf(stimulus, nlags, taxis, faxis, maxtm, maxsm)
% stimulus_to_tmf_smf Peak temporal/spectral modulation value in stimulus
%
%    [sprtmf, sprsmf] = sm_stimulus_to_tmf_smf(stimulus, nlags, taxis, faxis)
%    
%    stimulus : Nfrequencies x Ntimebins stimulus matrix.
%       Saved in: dmr-50flo-40000fhi-4SM-150TM-40db-96kHz-96DF-30min_DFt5_DFf8-matrix.mat
%    
%    nlags : number of time bins in each trial (time point). This should match the number
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
%       dmr-50flo-40000fhi-4SM-150TM-40db-96kHz-96DF-30min_DFt5_DFf8-mtf.mat 
%    
fprintf('%s\n', mfilename);

%sprtmf = zeros(1, size(stimulus,2));
sprsmf = zeros(1, size(stimulus,2));
sprphase = zeros(1, size(stimulus,2));

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
    
    % get ripple tranfer function of the stimulus frame at the moment
    %[tmfaxis, smfaxis, rtf] = sm_mtf_sta2rtf(frame, taxis, faxis, maxtm, maxsm, 15, 15);
    [smfaxis, rtf, phase] = sm_mtf_sta2rtf(frame, taxis, faxis, maxtm, maxsm);
    %[tmfaxis, rtf, phase] = sm_mtf_sta2rtf(frame, taxis, faxis, maxtm, maxsm);
    % find the maximum modulation frequncy in rtf
    %indfmax = 5;
    %[~,indtmax] = find(rtf == max(rtf(indfmax,:)));
    
    %[indfmax, indtmax] = find(rtf == max(rtf(:)));
    %xmf_max = max(smfaxis(indfmax)); %incase there are two maximum values in rtf, take the larger modulation rate
    
    t0 = 200;
    [indfmax,~] = find(rtf == max(rtf(:,t0)));
    xmf_max = smfaxis(indfmax);
    taxis = 1:nlags;
    
    %tmf_max = max(tmfaxis(indtmax));
    %saxis = faxis;
    
    %phase_max_lst = [];
    %for j = 1:length(indtmax)
    %    phase_max_lst = [phase_max_lst phase(indfmax(j), indtmax(j))];
    %end
    %k = find(phase_max_lst == max(phase_max_lst));
    %phase_max = phase(indfmax(k), indtmax(k));

    %phase_max = phase(indfmax, indtmax);
    %phase_max = std(phase(indfmax, :));
    %phase_max = mean(phase(indfmax, :));
    phase_max = phase(indfmax, t0);

%     dmax = 0;
%     [indfmax,~] = find(rtf == max(rtf(:)));
%     indtmax = 0;
%     d_lst = [];
%     for n = 1:(size(phase,2)-1)
%         d = abs(phase(indfmax,n+1)-phase(indfmax,n));
%         if d <= 180
%             d_lst = [d_lst d];
%             if d >= dmax
%                 dmax = d;
%                 indtmax = n+1;
%             end
%         else
%             d = 360 -d;
%             d_lst = [d_lst d];
%             if d >= dmax
%                 dmax = d;
%                 indtmax = n+1;
%             end
%         end
%     end      
%     xmf_max = smfaxis(indfmax); 
%     phase_max = phase(indfmax,indtmax) - phase(indfmax,indtmax-1);
    %phase_max = max(d_lst);
    %d_lst(find(d_lst==max(d_lst)))=[];
    %phase_max = mean(d_lst);
    %phase_max = max(d_lst)/mean(d_lst);
    
%     dmax = 0;
%     [indfmax,~] = find(rtf == max(rtf(:)));
%     indtmax = 0;
%     d1_lst = [];
%     d2_lst = [];
%     for n = 1:size(phase,2)
%         if indfmax ~= 21
%             d1 = abs(phase(indfmax,n)-phase(indfmax+1,n));
%             d1_lst = [d1_lst d1];
%         else
%             d1_lst = [0];
%         end
%         if indfmax ~= 1
%             d2 = abs(phase(indfmax,n)-phase(indfmax-1,n));
%             d2_lst = [d2_lst d2];
%         else
%             d2_lst = [0];
%         end
%     end
%     xmf_max = smfaxis(indfmax); 
%     phase_max = max([mean(d1_lst) mean(d2_lst)]);
    
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

%     figure;
%     imagesc(phase);
%     clim = max(abs(phase(:)));
%     caxis([-clim clim]);
%     axis xy
%     colorbar
    
    %sprtmf(i) = tmf_max;
    sprsmf(i) = xmf_max;
    sprphase(i) = phase_max;
end

return;


