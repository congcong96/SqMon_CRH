% Makes the dynamic moving ripple using Monty Escabi's toolbox.

% Declare variables.
flo       = 500;        % lower carrier frequency
fhi       = 40000;      % upper carrier frequency 
fRD      = 0.2;        % ripple density bandlimit frequency
fFM      = 0.6;        % temporal modulation bandlimit frequency
MaxRD    = 4;          % maximum ripple density (cycles/octave)
MaxFM    = 40;			  % maximum modulation frequency (Hz)


% App      = [10 20 30 40 50 30];	  % peak to peak ripple amplitude (dB)
                                          % will make different
                                          % stimulus for each one

App      = [40];       % peak to peak ripple amplitude (dB)

Fs       = 96000;      % sampling rate

minutes  = 1;
M        = minutes * 60 * Fs;   % number of samples

NS       = 316;        % number of sinusoid carriers
                       % total (10 * MaxRD * #octaves)

NB       = 1;       % number of blocks to divide parameter space
         			  %  into, or number of ripple profiles to add;
         			  %  note that number of ripple components 
         			  %  is NBxNB

Axis     = 'log';	  % carrier frequency axis, 'log' or 'lin',
         			  %  default 'lin'

Block    = 'n';    % breaks up the Fm vs. RD parameter space 
         			  %  into NBxNB discrete blocks, 'y' or 'n' 
          			  %  default 'n'

DF       = 48;    % downsampling factor for spectral profile,
 
AmpDist  = 'dB';		  % modulation amplitude distribution
                                  %  'dB'  = uniformly distributed, dB scale
                                  %  'lin' = uniformly distributed, lin scale
                                  %  'both' = designs signals with both, and
                                  %           uses the last element in the 
                                  %           App array to designate the 
                                  %           modulation depth for 'lin'


filename = sprintf('dmr-%.0fflo-%.0ffhi-%.0fSM-%.0fTM-%.0fdb-%.0fDF-%.0fmin', ...
                   flo, fhi, MaxRD, MaxFM, App, DF, minutes)

pause(1);


ripnoise(filename, flo, fhi, fRD, fFM, MaxRD, MaxFM, App, M, Fs, NS, NB, Axis, Block, DF, AmpDist);









