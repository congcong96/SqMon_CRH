function plot_strf(rf, n0, taxis, faxis, p, mdb, dur)
% INPUT
%   rf: strf Nfreq x Mtbins matrix
%   n0: number of spikes of the unit
%   taxis: time axis of the rf
%   faxis: frequency axis of the rf
%   p: threshold for significance of values in rf
%   mdb: modulation depth of the stimulus
%   dur: duration of the stimulus in seconds
% for proper color dispalt, download the 'cbrewer' package from matlab
fticks = [];
for ii = -1:1:5
    fticks = [fticks find(faxis/100 >= 2^ii*10, 1)];
end
tticks = [];
for ii = 0:25:100%-1:2:5
    tticks = [tticks find(taxis*1000 >= ii, 1)];
end
%set parametrs
if nargin == 4
    p = 0.002;
    mdb = 40;
    dur = 60 * 15;
end
%rf = rf(fticks(1)-20:fticks(end)-20, tticks(1):tticks(end));
[rfsig] = significant_strf(rf, p, n0, mdb, dur);
%rfsig = rfsig(:, tticks(1):tticks(end));
boundary = max(abs(rfsig(:)));

imagesc(rfsig);
cmap = cbrewer('div', 'RdBu', 100);
colormap(gca, flipud(cmap))

set(gca,'ydir', 'normal');
set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);

xlabel('time/ms')
%tticks = linspace(0, length(taxis), 6);
tidx = find(tticks > 1);
xticks(tticks(tidx))
xlabels = [0 25 50 75 100];
xticklabels(xlabels(tidx))
%xticks(tticks+1)
%xticklabels(taxis(tticks(1:end-1) + 1)*1000)

ylabel('frequency/kHz')
fidx = find(fticks > 1);

yticks(fticks(fidx))
yticklabels(0.25 * 2.^fidx)
%yticklabels([2 4 8 16])
%yticklabels([0.5 2 16 32])

