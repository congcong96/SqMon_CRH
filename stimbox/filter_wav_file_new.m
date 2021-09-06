function filter_wav_file_new(wavfile,es,box,ms)
% filter_ripple_noise(wavfile) Filter wav file with signal on left
% channel.
%
% Sampling rate must be 96000 Hz. Signal is assumed to be on
% the left channel.
%

narginchk(1,4);

index = findstr(wavfile, '.wav');
if isempty(index)
    error('Need to use file with .wav ending.');
end

if nargin ==1
    wav_outfile = sprintf('%s-filtered.wav', wavfile(1:index-1));
elseif nargin == 2 && es == 1
    wav_outfile = sprintf('%s-filtered-es.wav', wavfile(1:index-1));
elseif nargin == 3 && es == 0 && box == 1
    wav_outfile = sprintf('%s-filtered-box.wav', wavfile(1:index-1));
elseif nargin == 4 && es == 0 && box == 0 && ms == 1
    wav_outfile = sprintf('%s-filtered-ms.wav', wavfile(1:index-1));
end

% [y, fs, nbits] = wavread(wavfile);
audiofileinfo = audioinfo(wavfile);
% nsamples = audiofileinfo.TotalSamples;
% nchannels = audiofileinfo.NumChannels;
fs = audiofileinfo.SampleRate;
nbits = audiofileinfo.BitsPerSample;

if fs ~= 96000
    error('Sampling rate must be 96000');
end

if nbits ~= 16
    error('Resolution must be 16 bits.');
end

y = audioread(wavfile);

y = y(:,1);

if fs ~= 96000
    error('Sampling rate must be 96000');
end

if nargin == 1
    b = speaker_calibration_filter;
elseif nargin == 2 && es == 1
    b = electrostatic_speaker_calibration_filter;
elseif nargin == 3 && es == 0 && box == 1
    b = soundbox_speaker_calibration_filter;
elseif nargin == 4 && es == 0 && box == 0 && ms == 1
    b = speaker_calibration_filter_congcong;
end

yfilt = filtfilt(b,1,y);

%wavwrite(yfilt, fs, nbits, wav_outfile);
audiowrite(wav_outfile,yfilt,fs,'BitsPerSample',nbits);

freqz(b,1,1024,fs);

figure;
freqz(yfilt,1,1024,fs);

return;



