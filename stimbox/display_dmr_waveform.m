function display_dmr_waveform(SpecFile, MdB, ModType, sprtype)
%  FILE NAME    : DISPLAY DYNAMIC MOVING RIPPLE WAVEFORM
%  DESCRIPTION  : Displays a waveform from the spectral-temporal envelope of the stimulus in 'Specfile' as waveform
%
%	SpecFile	: Spectral Profile File (*.spr file)
%
%	MdB		: Signal Modulation Index in dB (usually 40)
%
%	ModType	: Kernel modulation type : 'lin' or 'dB'
%             Default is 'dB'
%             **** You will be using 'dB' ****
%
%	sprtype  : SPR File Type : 'float' or 'int16'
%			     Default == 'float' **** Use This ****
%
% You can run the program, and use the default parameters of 40, 'dB', and
% 'float' by simply issuing the following command:
%
% display_dynamic_moving_ripple_envelope(SpecFile);
%
% The default values will be used and assigned inside the function.
% diaplay 10s at a time
% -- updated by congcong, 11/05/2019
%Parameters
if ( nargin == 1 )
    MdB = 40;
    ModType = 'dB';
    sprtype = 'float';
end

if ( nargin == 2 )
    ModType = 'dB';
    sprtype = 'float';
end

if ( nargin == 3 )
    sprtype = 'float';
end

%Loading Parameter Data
index = findstr(SpecFile,'.spr');
ParamFile = [SpecFile(1:index(1)-1) '_param.mat'];
load(ParamFile, 'NT', 'NF', 'taxis', 'faxis');

%Opening Spectral Profile File
fid = fopen(SpecFile);

if strcmp(ModType,'dB')
    
    RMSP = -MdB/2; % RMS value of normalized Spectral Profile
    
    %Initializing First and Second Spectral Profile Segments
    frewind(fid);
    
    dt = diff(taxis);
    dt = dt(1)*1000;
    
    
    doct = log2(faxis(end)/faxis(1))/length(faxis);
    df = (1)/doct;
    ftick = round(1:df:length(faxis));
    freq = round(faxis(ftick));
    
    
    % Displaying Profile
    while (~feof(fid))
        if strcmp(sprtype,'float')
            S1 = fread(fid,NT*NF*30,'float');% 10s, 0.33s/block
        else
            S1 = fread(fid,NT*NF*30,'int16')/.99/1024/32/2-.5;
        end
        
        S1 = reshape(S1,NF,length(S1)/NF);
    
        profile = S1 - RMSP;
        ncols = size(profile,2);
        time = (0:ncols-1) * dt;
        
        plot(sum(profile)) %plot(sum(profile),'LineWidth',2)
        box off
        set(gca, 'xtick', 1:5000:ncols, 'xticklabel', time(1:5000:ncols)/1000);
        %set(gca, 'ytick', ftick, 'yticklabel', freq);
        %set(gca, 'tickdir', 'out');
        xlabel('Time (s)');
        ylabel('Power');
        fprintf('Push <Enter> to display next segment.\n');
        pause
        % Reading Spectral Profile Data File
    end % (while)
    
elseif ( strcmp(ModType,'lin') )
    error('You shouldn"t be using ripples with linear amplitude scales.');
end % (if)

% Closing all opened files
fclose all;

return;