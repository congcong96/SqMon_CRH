function [stimulus,taxis,faxis] = read_entire_spr_file_2(sprfile,nblocks)
%[stimulus] = read_entire_spr_file(sprfile, ncols)
%
% sprfile : Spectral Profile File. Ends in .spr
% This function assumes the file is saved as 'float'.
%
% stimulus : mxn matrix, which is the entire dmr envelope file.
%
% Each .spr file has an accompanying _param.mat file. The param file must
% reside in the same folder as the sprfile.
%
%
% caa 12/15/06

% added taxis,faxis for output and changed the name from
% read_entire_spr_file to read_entire_spr_file_2
% as I have already one called read_entire_spr_file in infobox
% Natsumi 05June18

narginchk(1,2);

if nargin == 1
    nblocks = inf;
end

Sound = 'MR';
ModType = 'dB';
sprtype='float';


%Loading Parameter Data
index = findstr(sprfile,'.spr');
paramFile = [sprfile(1:index(1)-1) '_param.mat'];
% f = ['load ' paramFile];
% eval(f);
load(paramFile);
clear App  MaxFM XMax Axis MaxRD RD f phase Block Mn RP f1 f2 
clear Mnfft FM N fFM fRD NB NS LL filename M X fphase Fsn

 
%Opening Spectral Profile File
fid = fopen(sprfile);
frewind(fid);
MdB = 40;
RMSP = -MdB/2;
stimulus = [];
blocknum = 1;
while ( ~feof(fid) & (blocknum <= nblocks) )
   [s1, count] = fread(fid,NT*NF,'float');
   if ( count )
      s1 = reshape(s1,NF,NT);
      s1 = MdB * s1 - RMSP;
      stimulus = [stimulus s1];
   end % (if)
    blocknum = blocknum + 1;
end % (while)

[nf, ntrials] = size(stimulus);

fprintf('\n#Freqs = %.0f, #trials = %.0f\n\n', nf, ntrials);


%Closing all opened files
fclose('all');


return;



