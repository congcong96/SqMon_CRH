%function [stimulus] = read_spr_file(sprfile)
%
%       FILE NAME   : read_spr_isk_file
%       DESCRIPTION : Save entire ripple file as a matrix  
%
% sprfile :Spectral Profile File
%
% RETURNED VALUES 
%
% stimulus : The entire sprfile as a matrix.
%
% caa 12/15/06
function [stimulus] = read_spr_file(sprfile)


% Determine the length of the ripple envelope file
[n, NT, NF] = ripple_stim_length(sprfile);
sprlength = n / NF; % number of trials in stimulus file


% if ( length(locator) ~= sprlength )
%    error('Spike train length and envelope file length have different number of trials.');
% end



%Loading Parameter Data
index = findstr(sprfile,'.spr');
paramFile = [sprfile(1:index(1)-1) '_param.mat'];
f = ['load ' paramFile];
eval(f);
clear App  MaxFM XMax Axis MaxRD RD f phase Block Mn RP f1 f2 
clear Mnfft FM N fFM fRD NB NS LL filename M X fphase Fsn

 
%Opening Spectral Profile File
fid = fopen(sprfile,'r');
frewind(fid);

tic
stimulus = [];
n = 1;
while ( ~feof(fid) )
   [s1, count] = fread(fid,NT*NF,'float');
   if ( count )
      s1 = reshape(s1,NF,NT);
      stimulus = [stimulus s1];
%       fprintf('n = %.0f\n', n);
%       n = n + 1;
   end % (if)
end % (while)
toc

stimulus = stimulus - mean ( mean ( stimulus ) );


%Closing all opened files
fclose('all');

return;




