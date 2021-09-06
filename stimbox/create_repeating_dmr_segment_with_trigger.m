function create_repeating_dmr_segment_with_trigger(stim, fs, delay, nreps)
%
% create_repeating_dmr_segment_with_trigger(stim, delay)
%
% Takes a stimulus segment and creates a signed word file with
% nreps repeats of the stimulus, where each segment has 
% 1.0 sec cosine squared rise/fall ramps.
%
% stim : either a vector or a file name specifying where a stimulus
%        composed of float values is stored
%
% fs : sampling rate of the dmr segment.
%
% delay : time, in sec, between end of 1 stimulus
%         presentation and the start of the next presentation
%
% nreps : number of repeats of the dmr segment. The default is 100.
%
% caa 12/14/05
% 100ms trigger, only up
% -- updated by Congcong, 11/04/2019


if ( nargin < 3 || nargin > 4 )
   error('You need 3 or 4 input arguments.');
end

if ( nargin == 3 )
   nreps = 100;
end


if ischar(stim) % stim is a file name
   filename = stim;
   fidin = fopen(stim,'r');
   [data, ~] = fread(fidin,'float');
   inddot = strfind(stim,'.');
   prefix = stim(1:inddot-1);
   clear('stim');
   stim = data(:)';
else
   filename = '';
   prefix = 'signal';
   stim = stim(:)';
end


stimlen = length(stim) + ceil( delay*fs );

% time = (0:length(stim)-1)/fs;
% plot(time, stim)
% length(stim)
% pause


% Define the envelope
rampUpTime = 1.0; % 1.0 sec ramps
rampDnTime = 1.0;
rampUpLen = floor(rampUpTime*fs);
rampDnLen = floor(rampDnTime*fs);

rampUpEnv = 0.5*(1-cos(2*pi*(0:(rampUpLen-1))/(2*rampUpLen)));
rampDnEnv = 0.5*(1-cos(2*pi*(rampDnLen:-1:1)/(2*rampDnLen)));


% normalize the stimulus
mx = max(max(abs(stim)));
signal = round(0.999*32767) ./ mx .* stim; % normalized for int16


% Apply the cosine-squared ramp:
signal(1:rampUpLen) = signal(1:rampUpLen) .* rampUpEnv;
signal((end-rampDnLen+1):end) = signal( (end-rampDnLen+1):end ) .* rampDnEnv;
signal = round(signal);


% Append zeros to make appropriate delay between repeats
signal = [signal zeros(1, stimlen-length(signal))];


% Now make the trigger: 10 ms pulse.
up_dn_pulse = [ones(1,round(0.01*fs)) zeros(1,round(0.01*fs))];
pulselen = length(up_dn_pulse);

trigger = round(0.5 * 32767) * [up_dn_pulse zeros(1, stimlen-pulselen)];
firsttrigger = round(0.5 * 32767) * [up_dn_pulse up_dn_pulse up_dn_pulse zeros(1, stimlen-3*pulselen)];


% plot the signals just to make sure they are correct:

time = (0:length(signal)-1)/fs;

subplot(3,1,1);
plot(time, signal);
axis([min(time)-0.05*max(time) 1.05*max(time) -35000 35000]);

subplot(3,1,2);
plot(time, trigger);
axis([min(time)-0.05*max(time) 1.05*max(time) -35000 35000]);

subplot(3,1,3);
plot(time, firsttrigger);
axis([min(time)-0.05*max(time) 1.05*max(time) -35000 35000]);


% Create the output files:
stimtrack = sprintf('%s_%.0frepeats_%.0fmsdelay_%.0fkhz_stim', prefix, nreps, 1000*delay, round(fs/1000))
trigtrack = sprintf('%s_%.0frepeats_%.0fmsdelay_%.0fkhz_trig', prefix, nreps, 1000*delay, round(fs/1000))

fidstimtrack = fopen([stimtrack '.sw'], 'w');
fidtrigtrack = fopen([trigtrack '.sw'], 'w');


% write out some initial zeros
fwrite(fidstimtrack, zeros(1,ceil(10*fs)), 'int16');
fwrite(fidtrigtrack, zeros(1,ceil(10*fs)), 'int16');


for i = 1:nreps

   fprintf('%.0f of %.0f\n', i, nreps);

   % write the data out to file
   if ( length(signal) ~= length(trigger) )
      error('The output files must be the same length.');
   else
      fwrite(fidstimtrack, signal, 'int16');
      if ( i == 1 )
         fwrite(fidtrigtrack, firsttrigger, 'int16');
      else
         fwrite(fidtrigtrack, trigger, 'int16');
      end
   end

end % (for i)


% Write out one minute of silence, and then close everything
fwrite(fidstimtrack, zeros(1,ceil(60*fs)), 'int16');
fwrite(fidtrigtrack, zeros(1,ceil(60*fs)), 'int16');

fclose(fidstimtrack);
fclose(fidtrigtrack);


% now print out the stimulus parameters for later decoding
txtfile = sprintf('%s_%.0frepeats_%.0fmsdelay_%.0fkhz', prefix, nreps, 1000*delay, round(fs/1000))
fidtxt = fopen([txtfile '.txt'], 'w');

fprintf(fidtxt, '\n');
fprintf(fidtxt, 'File used was: %s\n', filename);
fprintf(fidtxt, 'Ramps were cosine squared\n');
fprintf(fidtxt, 'Stimulus segment length was: %.2f seconds\n', length(stim) / fs);
fprintf(fidtxt, 'onset/offset ramps were %.0f/%.0f ms long\n', 1000*rampUpTime, 1000*rampDnTime);
fprintf(fidtxt, 'Delay between stimulus offset/onset was %.3f ms\n', round(1000*delay) );
fprintf(fidtxt, 'amplitude ranges -- trigger: [0,2^15/2], stimulus: +/- 2^15\n');
fprintf(fidtxt, 'The first trigger was a triple trigger\n\n');
fprintf(fidtxt, '\n');

fclose(fidtxt);



