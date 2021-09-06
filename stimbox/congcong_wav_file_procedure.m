%% congcong_wav_file_procedure.m
% this will create wav files that have both stim and trig vectors.
% then, a conpensation filtere will be added.
% Craig taught me how to do this on 23 Jan 17
% modified natsumi_wav_file_procedure.m
% for b = soundbox_speaker_calibration_filter
% 8April2019 Natsumi
% updated soundbox_wav_file_procedure for magnetic speaker
% 05/07/2019 Natsumi

work_path = 'E:\Congcong\Documents\stimulus';
%dmr_path = 'C:\Users\Natsumi\Google Drive\MATLAB\DMRs';
%% Part 1: calibration tones
cd(work_path)
folder = 'C:\Users\experiment\Stimuli\Congcong_RecordingStim\calibrationtones';
getfilename = dir(fullfile(folder,'*.wav'));
for idx = 1:length(getfilename)
    getfilename(idx).name
    wavfile = fullfile(folder,getfilename(idx).name);
    filter_wav_file_new(wavfile,0,0,1)
end

%% Part 1.5: puretones burst
cd(work_path)
folder = 'C:\Users\experiment\Stimuli\Congcong_RecordingStim\puretonesburst';
getfilename = dir(fullfile(folder,'*.wav'));
for idx = 1:length(getfilename)
    getfilename(idx).name
    wavfile = fullfile(folder,getfilename(idx).name);
    filter_wav_file_new(wavfile,0,0,1)
end

%% filter FRA stim
folder = 'E:\Congcong\Documents\stimulus\thalamus';
wavfile = fullfile(folder,'tms_freq_resp_area_center4khz_range6oct_nreps1_fs96000hz-stereo.wav');
filter_wav_file_new(wavfile,0,0,1)
%% add trigger for FRA
folder = 'E:\Congcong\Documents\stimulus\thalamus';
stim_file = 'tms_freq_resp_area_center4khz_range6oct_nreps1_fs96000hz-stereo-filtered-ms.wav';
folder2 = 'E:\Congcong\Documents\stimulus\thalamus';
trig_file = 'tms_freq_resp_area_center4khz_range6oct_nreps1_fs96000hz-trig-filtered.wav';

% [y, fs, nbits] = wavread(wavfile);
audiofileinfo = audioinfo(fullfile(folder,stim_file));
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

stim = audioread(fullfile(folder,stim_file));

% [y, fs, nbits] = wavread(wavfile);
audiofileinfo = audioinfo(fullfile(folder2,trig_file));
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

trig = audioread(fullfile(folder2,trig_file));

audiowrite(fullfile(folder,'tms_freq_resp_area_center4khz_range6oct_nreps1_fs96000hz-stim-trig-filtered.wav'),[stim(:) trig(:)],fs,'BitsPerSample',nbits);
%% Part 2: dmr
folder = 'E:\Congcong\Documents\stimulus\thalamus';
stim_file = 'rn1-500flo-40000fhi-0-4SM-0-40TM-40db-96khz-48DF-15min-seed190506_stim.wav';
trig_file = 'rn1-500flo-40000fhi-0-4SM-0-40TM-40db-96khz-48DF-15min-seed190506_trig.wav';

% [y, fs, nbits] = wavread(wavfile);
audiofileinfo = audioinfo(fullfile(folder,stim_file));
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

stim = audioread(fullfile(folder,stim_file));

% [y, fs, nbits] = wavread(wavfile);
audiofileinfo = audioinfo(fullfile(folder,stim_file));
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

trig = audioread(fullfile(folder,trig_file));

outfile = fullfile(folder, 'rn1-500flo-40000fhi-0-4SM-0-40TM-40db-96khz-48DF-15min-seed190506_stim_trig.wav');

audiowrite(outfile,[stim(:) trig(:)],fs,'BitsPerSample',nbits);

filter_ripple_noise(outfile,0,1);

%% Part3; exposure stim calibaration
cd(work_path)
folder = 'C:\Users\experiment\Stimuli\Congcong_RecordingStim\exposurecalibrationtones';
getfilename = dir(fullfile(folder,'*.wav'));
for idx = 1:length(getfilename)
    getfilename(idx).name
    wavfile = fullfile(folder,getfilename(idx).name);
    filter_wav_file_new(wavfile,0,0,1)
end

%% Part 4: exposure stim for recording
folder = 'C:\Users\experiment\Stimuli\Congcong_RecordingStim\originalfiles';
stim_file = 'rat_two_filter_expstim_SAMfr20Hz-SAMdep100-SAMdur150ms-TSdf05oct-TSdur50ms-F1F2df03oct-F1F2dur50ms-05octint-rep40-96khz-stim.wav'
trig_file = 'rat_two_filter_expstim_SAMfr20Hz-SAMdep100-SAMdur150ms-TSdf05oct-TSdur50ms-F1F2df03oct-F1F2dur50ms-05octint-rep40-96khz-trig.wav'

% [y, fs, nbits] = wavread(wavfile);
audiofileinfo = audioinfo(fullfile(folder,stim_file));
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

stim = audioread(fullfile(folder,stim_file));

% [y, fs, nbits] = wavread(wavfile);
audiofileinfo = audioinfo(fullfile(folder,trig_file));
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

trig = audioread(fullfile(folder,trig_file));

outfile = fullfile(folder, 'rat_two_filter_expstim_SAMfr20Hz-SAMdep100-SAMdur150ms-TSdf05oct-TSdur50ms-F1F2df03oct-F1F2dur50ms-05octint-rep40-96khz_stim_trig.wav');

audiowrite(outfile,[stim(:) trig(:)],fs,'BitsPerSample',nbits);

filter_ripple_noise(outfile,0,1);



%% add 10 sec for exposure stim (20May2019)
folder = 'C:\Users\experiment\Stimuli\Congcong_RecordingStim';
infile = 'rat_two_filter_expstim_SAMfr20Hz-SAMdep100-SAMdur150ms-TSdf05oct-TSdur50ms-F1F2df03oct-F1F2dur50ms-05octint-rep40-96khz_stim_trig-filtered-ms.wav';

% [y, fs, nbits] = wavread(wavfile);
audiofileinfo = audioinfo(fullfile(folder,infile));
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

[stim] = audioread(fullfile(folder,infile));

tensecvec = zeros(10*fs,2);

newstim = [tensecvec; stim; tensecvec];

outfile = 'rat_two_filter_expstim_SAMfr20Hz-SAMdep100-SAMdur150ms-TSdf05oct-TSdur50ms-F1F2df03oct-F1F2dur50ms-05octint-rep40-96khz_stim_trig-filtered-ms-2.wav';

audiowrite(outfile,newstim,fs,'BitsPerSample',nbits);
%% filter dmrrep50
folder = 'E:\Congcong\Documents\stimulus\thalamus';
stim_file = 'rn1-500flo-40000fhi-0-0SM-2-15TM-40db-96khz-48DF-10sec-seed19110540dB_50repeats_1000msdelay_96khz_stim.wav'
trig_file = 'rn1-500flo-40000fhi-0-0SM-2-15TM-40db-96khz-48DF-10sec-seed19110540dB_50repeats_1000msdelay_96khz_trig.wav'

% [y, fs, nbits] = wavread(wavfile);
audiofileinfo = audioinfo(fullfile(folder,stim_file));
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

stim = audioread(fullfile(folder,stim_file));

% [y, fs, nbits] = wavread(wavfile);
audiofileinfo = audioinfo(fullfile(folder,trig_file));
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

trig = audioread(fullfile(folder,trig_file));

outfile = fullfile(folder, 'rn1-500flo-40000fhi-0-0SM-2-15TM-40db-96khz-48DF-10sec-seed191105_50repeats_1000msdelay_96khz_stim_trig.wav');

audiowrite(outfile,[stim(:) trig(:)],fs,'BitsPerSample',nbits);

filter_ripple_noise(outfile,0,1);

%% filter FRA10 and add trigger
% filter FRA10
folder = 'E:\Congcong\Documents\stimulus\thalamus';
wavfile = fullfile(folder,'freq_resp_area_stimulus_flo500Hz_fhi32000Hz_nfreq21_natten8_nreps10_fs96000_signal.wav');
filter_wav_file_new(wavfile,0,0,1)
% add trigger
stim_file = 'freq_resp_area_stimulus_flo500Hz_fhi32000Hz_nfreq21_natten8_nreps10_fs96000_signal-filtered-ms.wav';
trig_file = 'freq_resp_area_stimulus_flo500Hz_fhi32000Hz_nfreq21_natten8_nreps10_fs96000_trigger.wav';

% [y, fs, nbits] = wavread(wavfile);
audiofileinfo = audioinfo(fullfile(folder,stim_file));
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

stim = audioread(fullfile(folder,stim_file));

% [y, fs, nbits] = wavread(wavfile);
audiofileinfo = audioinfo(fullfile(folder,trig_file));
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

trig = audioread(fullfile(folder,trig_file));

audiowrite(fullfile(folder,stim_file),[stim(:) trig(:)],fs,'BitsPerSample',nbits);
%% broadband noise
% filter
folder = 'E:\Congcong\Documents\stimulus\thalamus';
wavfile = fullfile(folder,'broadbandnoise_flo500khz_fhi40000khz_natten4_10dBsteps_nreps50_fs96000hz-stim.wav');
filter_wav_file_new(wavfile,0,0,1)
%% add trigger for broadband noise
folder = 'E:\Congcong\Documents\stimulus\thalamus';
stim_file = 'broadbandnoise_flo500khz_fhi40000khz_natten4_10dBsteps_nreps50_fs96000hz-stim-filtered-ms.wav';
folder2 = 'E:\Congcong\Documents\stimulus\thalamus';
trig_file = 'broadbandnoise_flo500khz_fhi40000khz_natten4_10dBsteps_nreps50_fs96000hz-trig.wav';

% [y, fs, nbits] = wavread(wavfile);
audiofileinfo = audioinfo(fullfile(folder,stim_file));
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

stim = audioread(fullfile(folder,stim_file));

% [y, fs, nbits] = wavread(wavfile);
audiofileinfo = audioinfo(fullfile(folder2,trig_file));
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

trig = audioread(fullfile(folder2,trig_file));

audiowrite(fullfile(folder,stim_file),[stim(:) trig(:)],fs,'BitsPerSample',nbits);