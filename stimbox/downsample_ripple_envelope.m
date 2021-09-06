%
% [new_spr_file, new_param_file] = downsample_ripple_envelope(specfile, spl, mdb, Sound)
%
% Real Time spectro-temporal receptive field. Uses Lee/Schetzen 
% Aproach via Specto-Temporal Envelope for dB Amplitude Sound distributions 
%
%  specfile       :  Spectral Profile File ( *.spr file )
%  spl            :  Signal RMS Sound Pressure Level
%  mdb            :  Signal Modulation Index in dB
%  Sound          :  Sound Type 
%                    Moving Ripple  :  MR ( Default )
%                    Ripple Noise   :  RN
%
function [new_spr_file, new_param_file] = downsample_ripple_envelope(specfile, spl, mdb, Sound)


%Loading Parameter Data
index = findstr(specfile,'.spr');
suffix = specfile(1:index(1)-1);
paramfile = [suffix '_param.mat'];
f = ['load ' paramfile];
eval(f);

clear App  MaxFM XMax Axis MaxRD RD f phase Block Mn RP f1 f2 
clear Mnfft FM N fFM fRD NB NS LL filename M X fphase Fsn



% Add in the new downsampling calculations for a stimulus sampling rate of
% 44100 Hz, the sampling rate for CDs.
% This will give us about 230/7 = 33 carriers, or 6 carriers per octave, and
% thus a sampling ripple density of 6 cycles per octave, or
% a maximum of 3 cycles per octave for resolution
% The new total downsampling factor in the time domain will be 100, and
% along the spectral axis it will be 7;

dnt = 22; % temporal axis downsampling factor
dnf = 4;  % spectral axis downsampling - gives about 
          % 12 carriers per octave for this case

new_spr_file = [suffix '_12_carriers_per_octave_theo.spr'];
new_param_file = [suffix '_12_carriers_per_octave_param_theo.mat'];





% Finding Mean Spectral Profile and RMS Power
spln = spl - 10*log10(NF);					% Normalized SPL per frequency band

if strcmp(Sound,'RN')
   rmsp = -mdb/2;					% RMS value of normalized Spectral Profile
   pp = mdb^2/12;					% Modulation Depth Variance 
elseif strcmp(Sound,'MR')
   rmsp = -mdb/2;					% RMS value of normalized Spectral Profile
   pp = mdb^2/8;					% Modulation Depth Variance 
end


% Opening Spectral Profile File
fid = fopen(specfile);
frewind(fid);

while ( ~feof(fid) )

   [profile, count] = fread(fid, NT*NF, 'float');

   if ( count )

      % Get original spectral envelope segment,
      % scale it, and then downsample it
      profile = reshape(profile, NF, NT);
      %profile = mdb * profile - rmsp;
      profile = profile(1:dnf:end, 1:dnt:end); % capture downsampled envelope

      % Get size of new envelope segment
      lent = length(profile(1,:));
      lenf = length(profile(:,1));

      % Writing New Spectral Profile File as 'float' file
      tofloat(new_spr_file, reshape(profile, 1 , lent*lenf)); 

   end

   clear profile 

end % ( while ~feof(fid) )

fclose('all');


% Get a sample profiles from the original spectral envelope file
fid = fopen(specfile);
frewind(fid);
profile = fread(fid, NT*NF, 'float');
profile = reshape(profile, NF, NT);
fclose('all');

% Save important parameters for the new spectral envelope file
profile = profile(1:dnf:end, 1:dnt:end);

faxis = faxis(1:dnf:end);
taxis = taxis(1:dnt:end);

% newnt = length(taxis);
% newnf = length(faxis);

newnt = length(profile(1,:));
newnf = length(profile(:,1));

NT = newnt;
NF = newnf;

DF = dnt * DF; % redefine the temporal downsampling factor

clear('profile');

save(new_param_file);



