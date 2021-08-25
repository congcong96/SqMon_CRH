function plot_crh(crh, tmf, smf, colorflag)
% INPUT:
%   crh: Nsmf * Ntmf 2D matrix or linearized vector of the filter
%   tmf: tmf axis
%   smf: smf axis
%   colorflag: choose the plot to be 'r' red, 'b' blue, or 'p' purple
% need 'cbrewer' package
assert(size(crh, 1)*size(crh,2) == (length(tmf)) * (length(smf)))
if ~exist('colorflag', 'var')
    colorflag = 'r';
end
if size(crh, 1) == 1
    crh = reshape(crh, [length(smf), length(tmf)]);
end

imagesc(tmf, smf, crh)
if strcmp(colorflag, 'r')
    cmap = cbrewer('seq', 'Reds', 9);
elseif strcmp(colorflag, 'b')
    cmap = cbrewer('seq', 'Blues', 9);
elseif strcmp(colorflag, 'p')
    cmap = cbrewer('seq', 'Purples', 9);
end
colormap(gca, cmap)
clim = max(abs(crh(:)));
caxis([0 clim]);
axis xy

yticks(1:4)

xticks(-60:20:60)

xlabel('TMF (Hz)')
ylabel('SMF (oct/cyc)')