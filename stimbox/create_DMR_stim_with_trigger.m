%% create_DMR_stim_with_trigger
% minutes  = 15;
% seed = 191104;
% the rest is the same as rn1-500flo-40000fhi-0-4SM-0-40TM-40db-96khz-48DF-20min-seed180117
% also chenged only up pulse for triggers
seed = 200707;       %change based on date
cd('E:\Congcong\Documents\stimulus') 
% This to explore the way DMRs are generated. 3/Oct/16
% Added MinRD and MinFM 26/10/16

% the original script is make_ripple_noise.m from Jermyn.
% Makes the ripple noise using Monty Escabi's toolbox.

% Declare variables.
flo       = 500;        % lower carrier frequency
fhi       = 40000;     % upper carrier frequency 
fRD      = 0.2;        % ripple density bandlimit frequency
fFM      = 0.6;        % temporal modulation bandlimit frequency

App      = 40;       % peak to peak ripple amplitude (dB)

Fs       = 96000;      % sampling rate

minutes  = 5;
M        = minutes * 60 * Fs;   % number of samples

NS       = 316;        % number of sinusoid carriers
                       % total (10 * MaxRD * #octaves)

NB       = 1; %16;      % number of blocks to divide parameter space
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

MinRD = 0;
MaxRD = 4;
MinFM = 0;
MaxFM = 40;

filename = sprintf('rn%d-%.0fflo-%.0ffhi-%.0f-%.0fSM-%.0f-%.0fTM-%.0fdb-%.0fkhz-%.0fDF-%.0fmin-seed%d', ...
        NB,flo, fhi, MinRD, MaxRD, MinFM, MaxFM, App, Fs/1000, DF, minutes,seed)


ripnoisemin(filename, flo, fhi, fRD, fFM, MinRD, MaxRD, MinFM, MaxFM, App, M, Fs, NS, NB, Axis, Block, DF, AmpDist,seed);

%% change into wav file

float2sw(sprintf('%s40dB.bin',filename),sprintf('%s.sw',filename));
sw2wav (sprintf('%s.sw',filename), Fs, 16);

%% make stim and trig files
% based on 
% % create_150repeat_ic_dmr_segment_with_trigger(stim, delay)
% % caa 12/14/03
% %
% % Takes a stimulus segment and creates a signed word file with
% % the stimulus, where the segment has 750 msec cosine squared rise/fall ramps.
% filename_in = sprintf('%s.spr',filename);
% fidin = fopen(filename,'r');
% [data, count] = fread(fidin,'float');
% inddot = findstr(filename,'.');
% stim = data(:)';

% load signal
stim = audioread(sprintf('%s.wav',filename));
stim = stim';

trigger = zeros(size(stim));
timeidx = 1:Fs/3:length(stim); % trigger evety one thrid of second

stimlen = length(stim);  % 1 sec stimulus presentation rate

% plot(stim)
% length(stim)
% pause

% Define the envelope
rampUpTime = 0.75; % 750 msec ramps
rampDnTime = 0.75;
rampUpLen = floor(rampUpTime*Fs);
rampDnLen = floor(rampDnTime*Fs);

rampUpEnv = 0.5*(1-cos(2*pi*(0:(rampUpLen-1))/(2*rampUpLen)));
rampDnEnv = 0.5*(1-cos(2*pi*(rampDnLen:-1:1)/(2*rampDnLen)));

% normalize the stimulus and then apply the ramps
mx = max(max(abs(stim)));
signal = round(0.99*32767) ./ mx .* stim; % normalized ramp for int16

% Apply the cosine-squared ramp:
signal(1:rampUpLen) = signal(1:rampUpLen) .* rampUpEnv;
signal((end-rampDnLen+1):end) = signal( (end-rampDnLen+1):end ) .* rampDnEnv;
signal = round(signal);

% % Append zeros so there is dead time to stimulus
% signal = [signal zeros(1, stimlen-length(signal))];

% Now make the trigger: 100 ms pulse.
%single_up_down_pulse = [ones(1,ceil(0.01*Fs)) -ones(1,ceil(0.01*Fs))]; %
% chenged to only up pulse
single_up_down_pulse = [ones(1,ceil(0.01*Fs)) zeros(1,ceil(0.01*Fs))];
triple_up_down_pulse = [single_up_down_pulse single_up_down_pulse single_up_down_pulse];

for i = 1:length(timeidx)
    if i == 1
        trigger(1:length(triple_up_down_pulse)) = triple_up_down_pulse;
    else
        trigger(timeidx(i):timeidx(i)+length(single_up_down_pulse)-1) = single_up_down_pulse;
    end
end

trigger = round(0.5 * 32767) .* trigger;

% length(signal)
% length(trigger)

% plot the signals just to make sure they are correct:
subplot(2,1,1);
plot(signal);
axis([1 timeidx(4) -32767 32767]);
%axis([1 stimlen -32767 32767]);
subplot(2,1,2);
plot(trigger);
axis([1 timeidx(4) -32767 32767]);
%axis([1 stimlen -32767 32767]);
%pause

% Create the output files:
stimtrack = sprintf('%s_stim', filename);
fidstimtrack = fopen([stimtrack '.sw'], 'w');

trigtrack = sprintf('%s_trig', filename);
fidtrigtrack = fopen([trigtrack '.sw'], 'w');


% write out some initial zeros 10s
fwrite(fidstimtrack, zeros(1,ceil(10*Fs)), 'int16');
fwrite(fidtrigtrack, zeros(1,ceil(10*Fs)), 'int16');

% write the data out to file
if ( length(signal) ~= length(trigger) )
    error('The output files must be the same length.');
else
    fwrite(fidstimtrack, signal, 'int16');
    fwrite(fidtrigtrack, trigger, 'int16');
end

% write out some zeros for the end 10s
fwrite(fidstimtrack, zeros(1,ceil(10*Fs)), 'int16');
fwrite(fidtrigtrack, zeros(1,ceil(10*Fs)), 'int16');

fclose('all');

%% change into
sw2wav(sprintf('%s.sw',stimtrack));
sw2wav(sprintf('%s.sw',trigtrack));
