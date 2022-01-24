%% set path
addpath(genpath('/home/conghu/MatlabCodes/SqMo_cNE'))

datapath = '/data/congcong/SqMoPhys_Josh/cNE_analysis';
spkfolder = '/data/congcong/SqMoPhys_Josh/mountainsort/pydict/std4_dmr_thresh';
cd(datapath)
nefiles = dir('*-20dft.mat');
%% plot CCG of posi and neg neurons
ccgfolder = '/data/congcong/SqMoPhys_Josh/figure/cNE/ccg';
CCG = [];
tw = 100; %in ms
bin = 2; %in ms
for ii = 1:length(nefiles)
    
    % load data
    load(nefiles(ii).name, 'exp_site_nedata')
    exp = exp_site_nedata.exp;
    nedata = exp_site_nedata.nedata;
    spkfile = dir(fullfile(spkfolder, [exp, '*']));
    load(fullfile(spkfile.folder, spkfile.name), 'spktrain')
    spktrain = spktrain(:,1:end-1);
    spktrain = reshape(spktrain, [size(spktrain, 1), 2*bin, size(spktrain, 2)/2/bin]);%2ms resolution
    spktrain = squeeze(sum(spktrain,2));
    nNE = size(nedata.Patterns, 2);
    NEmembers = nedata.NEmembers;
    nfigure = 1;
    for jj = 1:nNE
        member = NEmembers{jj};
        pattern = nedata.Patterns(:, jj);
        neg_idx = member(pattern(member) < 0);
        if isempty(neg_idx)
            continue
        end
        posi_idx = member(pattern(member) > 0);
        posi_train = spktrain(posi_idx, :)';
        neg_train = spktrain(neg_idx, :)';
        for k1 = 1:length(neg_idx)
            for k2 = 1:length(posi_idx)
                nplot = (k1 - 1) * length(posi_idx) + k2;
                nplot = mod(nplot, 16);
                if nplot == 1
                    figure
                    figuresetup2savepdf(30, 30)
                elseif nplot == 0
                    nplot = 16;
                end
                xc = xcorr(neg_train(:,k1), posi_train(:,k2), tw/bin);
                c = corr(neg_train(:,k1), posi_train(:,k2));
                subplot(4,4, nplot)
                bar(-tw:bin:tw, xc, 1, 'FaceColor', 'k')
                hold on
                title(sprintf('#%d-#%d CC: %.3f', neg_idx(k1), posi_idx(k2), c))
                if max(xc) < 1
                    ylim([-1 1])
                end
                if nplot == 16
                    printPDFandPSC(gcf, fullfile(ccgfolder, sprintf('ccg_%s_cNE%d_%d', exp, jj, nfigure)))
                    nfigure = nfigure + 1;
                    close
                end
                CCG = [CCG, zscore(xc)];
            end
        end
        if nplot ~= 16
            printPDFandPSC(gcf, fullfile(ccgfolder, sprintf('ccg_%s_cNE%d_%d', exp, jj,nfigure)))
            nfigure = nfigure + 1;
            close
        end
    end
end
%%
figure
sem = std(CCG')/sqrt(size(CCG, 2));
% plot waveform
curve1 = mean(CCG, 2)' + sem;
curve2 = mean(CCG, 2)' - sem;
time = -tw:bin:tw;
hold on
a = fill([time, fliplr(time)], [curve1, fliplr(curve2)], 'k');
a.FaceColor = .5 * [1 1 1];
a.EdgeColor = 'none';
plot(time, mean(CCG, 2), 'k', 'linewidth', 2)
ylabel('z-score')
xlabel('lag(ms)')
title('correlation between positive and negative members')

%% plot CCG of posi neurons
CCG = [];
for ii = 1:length(nefiles)
    
    % load data
    load(nefiles(ii).name, 'exp_site_nedata')
    exp = exp_site_nedata.exp;
    nedata = exp_site_nedata.nedata;
    spkfile = dir(fullfile(spkfolder, [exp, '*']));
    load(fullfile(spkfile.folder, spkfile.name), 'spktrain')
    spktrain = spktrain(:,1:end-1);
    spktrain = reshape(spktrain, [size(spktrain, 1), 2*bin, size(spktrain, 2)/2/bin]);%2ms resolution
    spktrain = squeeze(sum(spktrain,2));
    nNE = size(nedata.Patterns, 2);
    NEmembers = nedata.NEmembers;
    for jj = 1:nNE
        member = NEmembers{jj};
        pattern = nedata.Patterns(:, jj);
        neg_idx = member(pattern(member) < 0);
        posi_idx = member(pattern(member) > 0);
        if isempty(neg_idx) || length(posi_idx) < 2
            continue
        end
        posi_train = spktrain(posi_idx, :)';
        cmb = nchoosek(1:length(posi_idx), 2);
        for kk = 1:size(cmb, 1)
            xc = xcorr(posi_train(:,cmb(kk,1)), posi_train(:,cmb(kk,2)), tw);
            CCG = [CCG, zscore(xc)];
        end
    end
end
%%
figure
sem = std(CCG')/sqrt(size(CCG, 2));
% plot waveform
curve1 = mean(CCG, 2)' + sem;
curve2 = mean(CCG, 2)' - sem;
time = -tw:bin:tw;
hold on
a = fill([time, fliplr(time)], [curve1, fliplr(curve2)], 'k');
a.FaceColor = .5 * [1 1 1];
a.EdgeColor = 'none';
plot(time, mean(CCG, 2), 'k', 'linewidth', 2)
ylabel('z-score')
xlabel('lag(ms)')
title('correlation among positive members')