function [signal, trigger, fs, nbits] = filter_ripple_noise(wavfile,es,ms)
% filter_ripple_noise(wavfile) Filter wav file with signal on left
% channel.
%
% Sampling rate must be 96000 Hz. Signal is assumed to be on 
% the left channel.
%
% Procedure for filtering a long .wav file and setting the trigger
% to be monophasic.
% 
% [yfilt, fs, nbits] = filter_ripple_noise(wavfile);
%
% [ytotal,fs,nbits] = wavread(wavfile);
% trigger = ytotal(:,2); % trigger on channel 2, second column
% trigger(trigger<0) = 0;
% new_wave_file Filtered/modified wavfile data
%
% wavwrite([yfilt trigger], fs, nbits, new_wave_file);
% updated to add new magnetic speaker filter ms
% Natsumi 05072019 

narginchk(1,3);

index = findstr(wavfile, '.wav');
if isempty(index)
    error('Need to use file with .wav ending.');
end

if nargin ==1
    wav_outfile = sprintf('%s-filtered.wav', wavfile(1:index-1));
elseif nargin == 2 && es == 1
    wav_outfile = sprintf('%s-filtered-es.wav', wavfile(1:index-1));
elseif nargin == 3 && es == 0 && ms == 1
    wav_outfile = sprintf('%s-filtered-ms.wav', wavfile(1:index-1));
end

fprintf('\nInfile  = %s\n', wavfile);
fprintf('Outfile = %s\n\n', wav_outfile);

%[y, fs, nbits] = wavread(wavfile);
wavfileinfo = audioinfo(wavfile); % updated for MATLAB 2017b (Natsumi 8Dec2017)
nbits = wavfileinfo.BitsPerSample;
[y, fs] = audioread(wavfile);
signal = y(:,1);
trigger = y(:,2);
trigger(trigger<0) = 0;


% trigger = ytotal(:,2); % trigger on channel 2, second column
% trigger(trigger<0) = 0;
% new_wave_file Filtered/modified wavfile data
%
% wavwrite([yfilt trigger], fs, nbits, new_wave_file);



if fs ~= 96000
    error('Sampling rate must be 96000');
end

% Amplitudes for filter. These values produce a fairly flat 
% calibration.
if nargin == 1
    b = speaker_calibration_filter;
elseif nargin == 2 && es == 1
    b = electrostatic_speaker_calibration_filter;
elseif nargin == 3 && es == 0 && ms == 1
    b = speaker_calibration_filter_congcong;
end

signal = filtfilt(b,1,signal);

%wavwrite([signal trigger], fs, nbits, wav_outfile);
audiowrite(wav_outfile,[signal trigger],fs,'BitsPerSample',nbits) %updated for MATLAB2017b (Natsumi 8Dec2017)


return;


  