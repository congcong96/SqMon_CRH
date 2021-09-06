% Makes the dynamic moving ripple using Monty Escabi's toolbox.
% This file is for the cat auditory cortex. For other stations, like
% the IC, the low carrier frequency and the maximum temporal modulation
% rate will need to be changed.
seconds  = 10;
seed = 191105;
% Declare variables.
flo       = 500;        % lower carrier frequency
fhi       = 40000;      % upper carrier frequency 
fRD      = 0.2;        % ripple density bandlimit frequency
fFM      = 0.6;        % temporal modulation bandlimit frequency		  
MinRD = 0;
MaxRD = 0.3; % maximum ripple density (cycles/octave)
MinFM = 2;
MaxFM = 15;% maximum modulation frequency (Hz)

App      = 40;       % peak to peak ripple amplitude (dB)

Fs       = 96000;      % sampling rate - 44100 is the rate for a cd


M        = seconds * Fs;   % number of samples

NS       = 316;        % number of sinusoid carriers
                       % total (10 * MaxRD * #octaves)

NB       = 1;      % number of blocks to divide parameter space
         			  %  into, or number of ripple profiles to add;
         			  %  note that number of ripple components 
         			  %  is NBxNB

Axis     = 'log';	  % carrier frequency axis, 'log' or 'lin',
         			  %  default 'lin'

Block    = 'n';    % breaks up the Fm vs. RD parameter space 
        			  %  into NBxNB discrete blocks, 'y' or 'n' 
          			  %  default 'n'

DF       = 48;    % downsampling factor for spectral profile. Should make this
                  % about 0.5 ms. sampling rate for sound system / sampling
                  % rate for recording system
 
AmpDist  = 'dB';		  % modulation amplitude distribution
                                  %  'dB'  = uniformly distributed, dB scale
                                  %  'lin' = uniformly distributed, lin scale
                                  %  'both' = designs signals with both, and
                                  %           uses the last element in the 
                                  %           App array to designate the 
                                  %           modulation depth for 'lin'

% seed = 2012; % just use the date the sound was made as the seed


filename = sprintf('rn%d-%.0fflo-%.0ffhi-%.0f-%.0fSM-%.0f-%.0fTM-%.0fdb-%.0fkhz-%.0fDF-%.0fsec-seed%d', ...
        NB,flo, fhi, MinRD, MaxRD, MinFM, MaxFM, App, Fs/1000, DF, seconds,seed)


ripnoisemin(filename, flo, fhi, fRD, fFM, MinRD, MaxRD, MinFM, MaxFM, App, M, Fs, NS, NB, Axis, Block, DF, AmpDist,seed);









