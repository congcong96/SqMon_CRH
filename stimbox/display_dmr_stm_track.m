
function display_dmr_stm_track(SpecFile)
%  function display_dynamic_moving_ripple_track(SpecFile)
%
%  FILE NAME: DISPLAY DYNAMIC MOVING RIPPLE TRACK
%  DESCRIPTION: Displays the spectral-temporal modulation values of the stimulus in 'Specfile'
%
%	SpecFile	: Spectral Profile File (*.mat file)
%
% Natsumi 3Oct16

%load SpecFile
load(SpecFile)

% calculate the length of the stimulus (in sec) 
LengthTotal = M/Fs;
LengthPoint = LengthTotal/length(RD);
Step=50; % number of points in a block

figure;
set(gca, 'xlim', [-40 40]);
set(gca, 'ylim', [0 4]);
title(sprintf('One step = %.2f sec (total %.f sec)',LengthPoint*Step, LengthTotal));

for i=1:Step:length(RD)-Step
    
    if i ~= length(RD)-Step
        
        RDi = RD(:,i:i+Step);
        FMi = FM(:,i:i+Step);
    else
        RDi = RD(:,i:end);
        FMi = FM(:,i:end);
    end
    
    hold on;
    plot(FMi,RDi,'.-r') % plot new block
    xlabel('Temporal modulation (Hz)');
    ylabel('Spectral modulation (cycle/oct)');
    
    fprintf('Push <Enter> to display next segment.\n');
    pause;
    
    plot(FMi,RDi,'.-','Color',[0.5 0.5 0.5]) % change the color as past blocks
end
end
