function b = speaker_calibration_filter_congcong
% updated speaker_calibration_filter on 05/07/2019
fs = 96000;
fs_half = fs / 2;
fstim = [0 125 250 500 1000 2000 4000 8000 16000 32000 40000 48000];
f = fstim / fs_half;
% for am2t recordings 05/07/2019
% amp = [1 1 1 10^(-5/20) 10^(-10/20)  10^(-15/20) 10^(-12/20) 10^(-10/20) 10^(-5/20) 1 1 0];
% for thalamus recording 11/05/2019
% amp = [1 1 1 10^(-2/20) 10^(-10/20)  10^(-18/20) 10^(-12/20) 10^(-10/20) 10^(-5/20) 1 1 0];
% for thalamus recording 11/18/2019
amp = [1 1 1 10^(-2/20) 10^(-10/20)  10^(-18/20) 10^(-13/20) 10^(-10/20) 10^(-5/20) 1 1 0];

b = fir2(500,f,amp);

end