function stimulus_presentation_order = create_fra_stim_simple(varargin)
% freq_resp_area_stimulus_simple.m
%
% Create tonal stimuli with different frequency/attenuation values.
% 
%
% Default parameters are set within the function.
%
% The parameters that you may set, if you like, are:
%
% fs : sampling rate
% nfreq : number of tone frequencies
% noctaves : number of octaves that frequencies will cover
% flow : lowest frequency
% maxatten : maximum attenuation value
% datten : attenuation step size
% nreps : # of times each freq/atten combination will be repeated
% isi : interstimulus interval
% tonedur : duration of tone
%
% Craig Atencio
% added input parser and save stim_param.mat file
% -- Congcong 11/04/2019

% define parameters for stimuli (default)
% input parser
p = inputParser;
addParameter(p,'fs',96000);% sampling rate
addParameter(p,'nfreq',45);% total number of tones
addParameter(p,'noctaves',6);% number of octaves
addParameter(p,'flow',500);% lowest frequency tone
addParameter(p,'maxatten',70);% maximum attenuation value in dB
addParameter(p,'datten',10);% attenuation value step size in dB
addParameter(p,'nreps',1);% number of repetitions for each freq/atten combination
addParameter(p,'isi',0.300);% interstimulus interval - in seconds
addParameter(p,'tonedur',0.05);% tone length - in seconds

parse(p,varargin{:});

fs = p.Results.fs; 
nfreq = p.Results.nfreq; 
noctaves = p.Results.noctaves; 
flow = p.Results.flow; 
maxatten = p.Results.maxatten; 
datten = p.Results.datten; 
nreps = p.Results.nreps; 
isi = p.Results.isi; 
tonedur = p.Results.tonedur; 


% Define the ramps
rampUpTime = 0.0050;  % 5 msec ramps
rampDnTime = 0.0050;
rampUpLen = floor(rampUpTime*fs);
rampDnLen = floor(rampDnTime*fs);

rampUpEnv = 0.5*(1-cos(2*pi*(0:(rampUpLen-1))/(2*rampUpLen)));
rampDnEnv = 0.5*(1-cos(2*pi*(rampDnLen:-1:1)/(2*rampDnLen)));


% define all frequencies and attenuations that will be used
frequency = round(flow .* 2.^(linspace(0,noctaves,nfreq)));
frequency = frequency(:);
attenuation = 0:datten:maxatten;
attenuation = attenuation(:);


% now create the frequency-attenuation combinations

freq_atten_order = zeros(length(frequency)*length(attenuation), 2);

for i = 1:length(frequency)
   temp = frequency(i) * ones(length(attenuation),1);
   freq_atten_order((i-1)*length(attenuation)+1:i*length(attenuation), :) = [temp attenuation];
end % (for i)


% create duration-frequency-attenuation combinations
nrows = size(freq_atten_order,1);
stimulus_presentation_order = zeros(nreps * nrows, 4);



% Replicate the stimulus combination vector for the desired # of repeats
for i = 1:nreps
   stimulus_presentation_order((i-1)*nrows+1:i*nrows, :) = ...
      [tonedur*ones(nrows,1) isi*ones(nrows,1) freq_atten_order];
end % (for i)


% now randomize the stimulus presentation order
rand('state', 1); % Just set the state of the random number generator
for i = 1:5 % randomize the list 5 times
   r = rand(1,size(stimulus_presentation_order,1));
   [~,ind] = sort(r);
   stimulus_presentation_order = stimulus_presentation_order(ind,:);
end



stimtrack = sprintf('freq_resp_area_stimulus_flo%.0fHz_fhi%.0fHz_nfreq%.0f_natten%.0f_nreps%.0f_fs%.0f_signal', ...
   min(frequency), max(frequency), length(frequency), length(attenuation), nreps, fs)
trigtrack = sprintf('freq_resp_area_stimulus_flo%.0fHz_fhi%.0fHz_nfreq%.0f_natten%.0f_nreps%.0f_fs%.0f_trigger', ...
   min(frequency), max(frequency), length(frequency), length(attenuation), nreps, fs)

fidstimtrack = fopen([stimtrack '.sw'], 'w');
fidtrigtrack = fopen([trigtrack '.sw'], 'w');

% write out 10 seconds of initial deadtime
fwrite(fidstimtrack, zeros(1,10*fs), 'int16');
fwrite(fidtrigtrack, zeros(1,10*fs), 'int16');


nstim = size(stimulus_presentation_order,1);

% Now make the trigger: 5 ms biphasic pulse.
stimlen = ceil(isi * fs);

up_dn_pulse = [ones(1,round(0.01*fs)) zeros(1,round(0.01*fs))];
pulselen = length(up_dn_pulse);

trigger = round(0.5 * 32767) * [up_dn_pulse zeros(1, stimlen-pulselen)];
firsttrigger = round(0.5 * 32767) * [up_dn_pulse up_dn_pulse up_dn_pulse zeros(1, stimlen-3*pulselen)];


for i = 1:nstim

   dur = stimulus_presentation_order(i,1);
   freq = stimulus_presentation_order(i,3);
   atten = stimulus_presentation_order(i,4);

   time = (0:dur*fs) ./ fs;
   tone = 0.999 .* 32767 .* 10.^(-atten./20) .* sin(2 * pi * freq .* time);

   % now modulate by the cosine squared ramp and add the deadtime
   tone(1:rampUpLen) = tone(1:rampUpLen) .* rampUpEnv;
   tone((end-rampDnLen+1):end) = tone((end-rampDnLen+1):end) .* rampDnEnv;
   tone = [tone zeros(1,stimlen-length(tone))];
   tone = round(tone);

   t = (0:stimlen-1)/fs;

%    clf;
%    subplot(2,1,1); 
%    plot(t, tone); 
%    axis([min(t) max(t) -33000 33000]);
%    title(sprintf('%.0f / %.0f  - dur=%.0f, freq=%.0f, atten=%.0f', i, nstim, 1000*dur, freq, atten));
%    set(gca,'xticklabel','');
% 
%    subplot(2,1,2); 
%    plot(t, trigger); 
%    axis([min(t) max(t) 0 33000]);
%    pause(0.25);

   % write the data out to file
   if ( length(tone) ~= length(trigger) )
      error('The two output signals must be the same length.');
   else
      if ( i == 1 )
         fwrite(fidtrigtrack, firsttrigger, 'int16');
      else
         fwrite(fidtrigtrack, trigger, 'int16');
      end
      fwrite(fidstimtrack, tone, 'int16');
   end

end % (for i)

% 5s deadtime after
fwrite(fidstimtrack, zeros(1,5*fs), 'int16');
fwrite(fidtrigtrack, zeros(1,5*fs), 'int16');

fclose(fidstimtrack);
fclose(fidtrigtrack);

% now print out the stimulus parameters for later decoding
txtfile = sprintf('freq_resp_area_stimulus_flo%.0fHz_fhi%.0fHz_nfreq%.0f_natten%.0f_nreps%.0f_fs%.0f_params', ...
   min(frequency), max(frequency), length(frequency), length(attenuation), nreps, fs);
fidtxt = fopen([txtfile '.txt'], 'w');

fprintf(fidtxt, '\n');
fprintf(fidtxt, 'ramps are cos^2\n');
fprintf(fidtxt, 'stimuli are %.0f ms tone bursts\n', 1000*tonedur);
fprintf(fidtxt, 'onset/offset ramps are %.0f/%.0f ms\n', 1000*rampUpTime, 1000*rampDnTime);
fprintf(fidtxt, 'interstimulus interval is %.0f ms\n', 1000*isi);
fprintf(fidtxt, 'amplitude ranges -- trigger: [0,2^23/2], stimulus: +/- 32767\n\n');

fprintf(fidtxt, '\n\n');
fprintf(fidtxt,'Presentation order of duration/frequency/attenuation tone burst combinations is:\n\n');
fprintf(fidtxt,' Stim      Dur       ISI      Freq     Atten\n');
fprintf(fidtxt,' num       (ms)      (ms)     (Hz)     (dB)\n');
fprintf(fidtxt, '\n');

for i = 1:size(stimulus_presentation_order,1)
   fprintf(fidtxt, '%5.0f%10.0f%10.0f%10.0f%10.0f\n', i, 1000*stimulus_presentation_order(i,1), ...
      1000*stimulus_presentation_order(i,2), stimulus_presentation_order(i,3), stimulus_presentation_order(i,4));
end

fprintf(fidtxt, '\n');
fclose(fidtxt);

%save stim_param.mat file
frequency = stimulus_presentation_order(:,3);
attenuation = stimulus_presentation_order(:,4); 
save([txtfile '.mat'], 'attenuation', 'frequency', 'dur', 'fs',...
    'rampDnTime', 'rampUpTime', 'isi', 'stimlen')

