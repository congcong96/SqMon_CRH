function [stimfile, trigfile, paramfile, frequency, attenuation] = ...
      acute_rat_tms_freq_resp_area_stimulus(center, range, nreps)
% tms_freq_resp_area_stimulus.m
%
% Create tonal stimuli that are used to construct
% frequency response areas. The stimuli created here
% are appropriate for cats, though not for squirrel monkeys.
% The frequency range does not go low enough for monkeys.
%
% [stimfile, trigfile, paramfile, frequency, attenuation] = ...
%       tms_freq_resp_area_stimulus(center, range, nreps)
%
% center : center frequency of FRA. In kHz.
%
% range : range of presented tones. A good choice is 4.
%
% Example: If center = 10 and range = 4 then 
%
%       [min(frequency) max(frequency)] =[2.6589 40.0000] Hz
%
% nreps : how many times each tone at a specific freq/atten will
%       be presented. For single unit data we usually make this 5,
%       which results in a sound file that is about 17 minutes long.
%       Using nreps = 3 may also be adequate for single units too.
%
% frequency and attenuation : a list of the tone frequencies and 
%       attenuations that were used to make the sound file. Each row
%       of [frequency attenuation] represents one presentation.
%       For each repetition there will be 675 combinations presented.
%
% Calling [stimfile, trigfile, paramfile, frequency, attenuation] = tms_freq_resp_area_stimulus;
% makes a stimulus with center = 10, range = 5, and nreps = 1.
%
% stimfile and trigfile are the names of the .sw files into which the tones
% and triggers were saved. These files may then be converted into .wav
% files using the function sw2wav.m
%
% paramfile is a text file holding all the specifications of the tones. The
% list of tone presentations is also listed in this file.
%
% The duration of each tone is 50 ms.
%
% The onset/offset ramps are 5 ms cosine-squared.
%
% The ISI is 300 ms.
%
% The sampling rate of the stimulus is 96 kHz.
%
% 10 s of deadtime precedes the stimuli. 60 s of deadtime occurs after
% all stimuli have been presented. 
%
% caa 11/1/07


if ( nargin == 0 )
   center = 5.6;
   range = 5;
   nreps = 1;
end



% define sampling rate and envelope properties
fs = 96000;  % DVD-Audio sampling rate
dur = 0.05;  % 50 ms tones
stimdur = 0.3;  % 300 ms ISI
stimlen = ceil(stimdur * fs);


% Define/Create the ramps
rampUpTime = 0.0050;  % 5 msec ramps
rampDnTime = 0.0050;
rampUpLen = floor(rampUpTime*fs);
rampDnLen = floor(rampDnTime*fs);

rampUpEnv = 0.5*(1-cos(2*pi*(0:(rampUpLen-1))/(2*rampUpLen)));
rampDnEnv = 0.5*(1-cos(2*pi*(rampDnLen:-1:1)/(2*rampDnLen)));


% Get tms frequencies and amplitudes
[frequency, amplitude] = tms_fra_params(center, range, 1); 
frequency = 1000 .* frequency; % convert to Hz


% Now get the repeats
frequency = repmat(frequency, 1, nreps);
frequency = frequency(:);

amplitude = repmat(amplitude, 1, nreps);
amplitude = amplitude(:);


% Get attenuation values
attenuation = amplitude - max(amplitude);


% Create output files
stimfile = sprintf('tms_freq_resp_area_center%.0fkhz_range%.0foct_nreps%.0f_fs%.0fhz_stim', ...
   center, range, nreps, fs);

trigfile = sprintf('tms_freq_resp_area_center%.0fkhz_range%.0foct_nreps%.0f_fs%.0fhz_trig', ...
   center, range, nreps, fs);

fidstimfile = fopen([stimfile '.sw'], 'w');
fidtrigfile = fopen([trigfile '.sw'], 'w');


% write out 10 seconds of initial deadtime
fwrite(fidstimfile, zeros(1,10*fs), 'int16');
fwrite(fidtrigfile, zeros(1,10*fs), 'int16');


% Now make the trigger: 5 ms biphasic pulse.

up_dn_pulse = [ones(1,round(0.005*fs)) -ones(1,round(0.005*fs))];
pulselen = length(up_dn_pulse);

trigger = round(0.5 * 32767) * [up_dn_pulse zeros(1, stimlen-pulselen)];
firsttrigger = round(0.5 * 32767) * [up_dn_pulse up_dn_pulse up_dn_pulse zeros(1, stimlen-3*pulselen)];


for i = 1:length(frequency)

   fprintf('#%.0f of %.0f\n', i, length(frequency));

   freq = frequency(i);
   amp = amplitude(i);
   atten = amp - max(amplitude);

   time = (0:dur*fs) ./ fs;
   tone = 0.999 .* 2^15 .* 10.^(atten./20) .* sin(2 * pi * freq .* time);

   % now modulate by the cosine squared ramp and add the deadtime
   tone(1:rampUpLen) = tone(1:rampUpLen) .* rampUpEnv;
   tone((end-rampDnLen+1):end) = tone((end-rampDnLen+1):end) .* rampDnEnv;
   tone = [tone zeros(1,stimlen-length(tone))];
   tone = round(tone);

   t = (0:(stimlen-1))/fs;

%    clf;
%    subplot(1,2,1); 
%    plot(t, tone); 
%    axis([min(t) max(t) -35000 35000]);
%    title(sprintf('%.0f - dur=%.0f, freq=%.0f, atten=%.0f', i, 1000*dur, freq, atten));
%    set(gca,'xticklabel','');
%    set(gca,'tickdir', 'out');
%    set(gca,'ticklength', [0.025 0.05]);
% 
%    subplot(2,2,2);
%    nfft = 1024;
%    tk = fft(tone, nfft);
%    tkmag = abs( tk(1:length(tk)/2) );
%    fk = (0:(nfft/2)-1) * fs / nfft;
%    plot(fk, tkmag);
%    set(gca,'tickdir', 'out');
%    set(gca,'ticklength', [0.025 0.05]);
%    xlabel('Frequency in Hz')
%    ylabel('Magnitude')
%    grid; 
% 
%    subplot(2,2,4); 
%    plot(t, firsttrigger); 
%    axis([min(t) max(t) -2^15/1.8 2^15/1.8]);
%    set(gca,'tickdir', 'out');
%    set(gca,'ticklength', [0.025 0.025]);
%    pause;

   % write the data out to file
   if ( length(tone) ~= length(trigger) )
      error('The two output files must be the same length.');
   else
      fwrite(fidstimfile, tone, 'int16');
      if ( i == 1 )
         fwrite(fidtrigfile, firsttrigger, 'int16');
      else
         fwrite(fidtrigfile, trigger, 'int16');
      end
   end

end % (for i)

% Write out 10 seconds of dead time
fwrite(fidstimfile, zeros(1,ceil(10*fs)), 'int16');
fwrite(fidtrigfile, zeros(1,ceil(10*fs)), 'int16');

% Close all files:
fclose(fidstimfile);
fclose(fidtrigfile);



% now print out the stimulus parameters for later decoding
param_matfile = sprintf('acute_rat_tms_freq_resp_area_center%.0fkhz_range%.0foct_nreps%.0f_fs%.0fhz_params', ...
   center, range, nreps, fs);

save(param_matfile, 'fs', 'dur', 'stimdur', 'stimlen', 'rampUpTime', 'rampDnTime', ...
'stimfile', 'trigfile', 'frequency', 'attenuation');



% now print out the stimulus parameters for later decoding
paramfile = sprintf('acute_rat_tms_freq_resp_area_center%.0fkhz_range%.0foct_nreps%.0f_fs%.0fhz_params', ...
   center, range, nreps, fs);

fidtxt = fopen([paramfile '.txt'], 'w');

fprintf(fidtxt, '\n');
fprintf(fidtxt, 'ramps are cos^2\n');
fprintf(fidtxt, 'stimuli are 50 ms tone bursts\n');
fprintf(fidtxt, 'onset/offset ramps are %.0f/%.0f ms\n', 1000*rampUpTime, 1000*rampDnTime);
fprintf(fidtxt, 'interstimulus interval is 300 ms\n');
fprintf(fidtxt, 'amplitude ranges -- trigger: [0,2^15/2], stimulus: +/- 2^15\n');
fprintf(fidtxt, 'The first trigger is a triple trigger\n\n');

fprintf(fidtxt, 'center = %.0f kHz\n', center);
fprintf(fidtxt, 'range = %.0f octaves\n', range);
fprintf(fidtxt, 'nreps = %.0f\n', nreps);

fprintf(fidtxt, '\n\n');
fprintf(fidtxt,'Presentation order of frequency/attenuation tone burst combinations is:\n\n');
fprintf(fidtxt,' Stim        Freq          Amp        Atten\n');
fprintf(fidtxt,' num         (Hz)         (dB)        (dB)\n');
fprintf(fidtxt, '\n');

for i = 1:length(frequency)
   fprintf(fidtxt, '%5.0f  %10.0f  %10.0f  %10.0f\n', i, frequency(i), ...
      amplitude(i), abs(attenuation(i)) );
end

fprintf(fidtxt, '\n');
fclose(fidtxt);


fclose('all');






function [freq, amp] = tms_fra_params(center, range, nreps)
%TMS_FRA_PARAMS Frequency/Amplitude parameters for TMS FRA stimulus
%
% TMS_FRA_PARAMS(CENTER,RANGE,NREPS) returns 45 frequencies, centered 
% at CENTER spanning a frequency extent RANGE defined in terms of number
% of octaves.
%
% There are 45 frequencies as well as 15 amplitude parameters, giving
% a total of 675 frequency/level combinations. The amplitudes are
% equally spaced over a 75 dB range. Each combination represents
% the characteristics of a specific tone stimulus.
%
% NREPS determines how many times each individual frequency/level
% combination is presented. For NREPS > 1, the each tone combination
% is played NREPS times consecutively.
%
% Examples: [f,a] = tms_fra_params(7,5,1) 
%           gives [min(f) max(f)] = [1.3365 39.5980];
%
%           [f,a] = tms_fra_params(10,4,1) 
%           gives [min(f) max(f)] = [2.6589 40.0000];
%
% with each possible tone presented once.
%
%
%TCCALC    tc = TCParam(center, range, spk)
%TCP.i(n,:) = [Fj,Ak]
%where fj and ak are the indices of the frequency and amplitude of the nth
%tone presented
%TCP.v(n,:) = [Fn, An]
%where Fn and An are the frequencies and relative amplitude of the ith tone
%presented.  The frequencies are evenly spaced on a log scale over range
%and centered at center.  The amplitudes are equally spaced spaced over a
%75dB range.


% Check input arguments for errors, etc.
error(nargchk(3,3,nargin));


%----- The following code is from Andrew Tan ------
TCf = [3.0 1.3, 2.6, 22.8, 26.6, 14.3, 3.3, 6.1, 1.0, 24.6, 6.4, 2.4, 1.6, ...
 31.1, 1.2, 3.8, 7.7, 7.1, 5.2, 1.7, 10.5, 15.4, 21.1, 4.4, 3.5, 8.9, 1.1, ...
 4.1, 16.7, 1.5, 4.8, 5.6, 9.7, 8.3, 1.9, 2.0, 11.3, 1.4, 18.0, 12.2, ...
 13.2,  2.2, 19.5, 2.8, 28.8];

TCa = [78, 67, 72, 41, 62, 46, 31, 21, 5, 52, 10, 15, 57, 26, 36];

[TCF, FOrder] = sort(TCf);
FIndex = zeros(size(FOrder));
for i = 1:length(FOrder)
   FIndex(FOrder(i)) = i;
end


[TCA, AOrder] = sort(TCa);
AIndex = zeros(size(AOrder));

for i = 1:length(AOrder)
   AIndex(AOrder(i)) = i;
end


k = 1;

for i = 1:length(AIndex)

   for j = 1:length(FIndex)

      jmod = mod(j+(i-1),length(AIndex));

      if jmod == 0
         jmod = length(AIndex);
      end

      kmod = mod(j, length(FIndex));

      if kmod == 0
         kmod = length(FIndex);
      end

      TCIndex(k,:) = [FIndex(kmod), AIndex(jmod)];
      k = k+1;

   end % (for i)

end % (for j)

TCP.i = TCIndex;

MinFre = center / 2^(range/2);
MaxFre = center * 2^(range/2);

LogFreIncrement = (log10(MaxFre) - log10(MinFre)) / length(FIndex);

%AmpRange in dB
AmpRange = 75;
AmpIncrement = 75 / length(AIndex);

for i = 1:(length(FIndex)*length(AIndex))

   TCP.v(i,1) = 10^(log10(MinFre)+TCIndex(i,1)*LogFreIncrement);
   TCP.v(i,2) = TCIndex(i,2) * AmpIncrement;

end

freq = TCP.v(:,1);
amp = TCP.v(:,2);
freqIndex = TCP.i(:,1);
ampIndex = TCP.i(:,2);

freqnew = repmat(freq(:)', nreps, 1);  % make nreps rows of freq
freq = freqnew(:);

ampnew = repmat(amp(:)', nreps, 1);  % make nreps rows of amp
amp = ampnew(:);




