function [phase_hist, smfaxis, phaseaxis] = calculate_phase_CRH(spktrain, mtffile)

load(mtffile,'sprphase', 'sprsmf', 'taxis', 'smfaxis')
%load(mtffile,'sprphase', 'sprtmf', 'saxis', 'tmfaxis')

% CRH of the neurons
%phase_hist = zeros(size(spktrain, 1), length(taxis) * length(smfaxis));
%phase_hist = zeros(size(spktrain, 1), length(saxis) * length(tmfaxis));
phase_hist = zeros(size(spktrain, 1), 12 * length(smfaxis));
%phase_hist_s = zeros(size(spktrain, 1), length(tmfaxis) * length(smfaxis));
if length(sprphase) > length(spktrain)
    downsample = round(length(sprphase)/length(spktrain));
    sprphase = sprphase(downsample:downsample:end);
    %sprtmf = sprtmf(downsample:downsample:end);
    sprsmf = sprsmf(downsample:downsample:end);
end
diff = length(sprphase) - length(spktrain);
if abs(diff) > 1
    error('stim subsample wrong!')
elseif diff > 0
    sprphase = sprphase(1:end-diff);
    %sprtmf = sprtmf(1:end-diff);
    sprsmf = sprsmf(1:end-diff);
elseif diff < 0
    spktrain = spktrain(:, 1:end+diff);
end

% phaseaxis1 = linspace(-180, 180, length(tmfaxis));
% phaseaxis2 = linspace(-180, 180, length(smfaxis));
phaseaxis = linspace(-180, 180, 12);
%phaseaxis = linspace(-180, 180, length(taxis));
%phaseaxis = linspace(-180, 180, length(saxis));

for ii = 1:size(spktrain, 1)
    %tmf = rude(spktrain(ii,:), sprtmf);
    smf = rude(spktrain(ii,:), sprsmf);
    phase = rude(spktrain(ii,:), sprphase);
    %phaseaxis = linspace(min(abs(phase(:))), max(abs(phase(:))), length(taxis));
    %phaseaxis = linspace(0, max(phase(:)), 40);
    temp = histcounts2(smf, phase, [smfaxis, 4], [phaseaxis, 180]);%Xedges(1) 是 x 维度的第一个 bin 的第一个边界，Xedges(end) 是最后一个 bin 的外边界
    %temp = histcounts2(phase,tmf, [phaseaxis, 180], [tmfaxis, 64]);
    phase_hist(ii,:) = temp(:)';
    %temp_s = histcounts2(smf, phase, [smfaxis, 4], [phaseaxis1, 180]);
    %phase_hist_s(ii,:) = temp_s(:)';
end
