%
%function []=ripnoisemin(filename,f1,f2,fRD,fFM,MinRD,MaxRD,MinFM,MaxFM,App,M,Fs,NS,NB,Axis,Block,DF,AmpDist,seed)
%
%	
%	FILE NAME 	: RIP NOISE MIN
%	DESCRIPTION 	: Generates Dynamic Moving Ripple and Ripple Noise
%			  Saved as a 'float' file
%			  Saves Spectral Profile as a sequence of MAT files
%			  This File is identical to RIPNOISE except it allows
%			  you to add an arbitrary number of Moving Ripple 
%			  Profiles to generate a Ripple Noise Profile
% changed ripnoiseb to add minumum ranges for FM and RD.
%
%	filename	: Ouput data file name
%       f1              : Lower carrier frequency
%       f2              : Upper carrier Frequency
%	fRD		: Ripple Density Bandlimit Frequency
%	fFM		: Temporal Modulation Bandlimit Frequency
%	MinRD		: Minimum Ripple Density ( Cycles / Octaves )
%	MaxRD		: Maximum Ripple Density ( Cycles / Octaves )
%	MinFM		: Minimum Modulation Frequency ( Hz )
%	MaxFM		: Maximum Modulation Frequency ( Hz )
%       App             : Peak to Peak Riple Amplitude ( dB )
%       M               : Number of Samples
%       Fs              : Sampling Rate
%	NS		: Number of Sinusoids Cariers used for ripple Noise
%	NB		: If 'block'=='n' this parameter coresponds to the 
%			  number of Moving Ripple Profiles which are added to
%			  generate a Ripple Noise Profile
%
%			  If 'block'=='y' the Ripple Density vs. Modulation 
%			  Rate space is divided into NBxNB squares.  When this
%			  is done each individual Moving Ripple Profile probes
%			  a small subspace of the parameters.  Not that for
%			  case a total of NBxNB Moving Ripple Envelopes are 
%			  superimposed to generate a Ripple Noise Envelope
%
%	Axis            : Carrier Freqeuncy Axis Type: 'log' or 'lin'
%			  Default = 'lin'
%	Block		: Breaks up the Fm vs. RD parameter space into
%			  NBxNB discrete blocks, 'y' or 'n', Default : 'n'
%			  See NB for Detailed Description
%	DF		: Temporal Dowsampling Factor For Spectral Profile
%			  Must Be an Integer
%	AmpDist		: Modulation Amplitude Distribution Type
%			  'dB'   = Uniformly Distributed on dB Scale
%			  'lin'  = Uniformly Distributed on linear Scale
%			  'both' = Designs Signals with both
%				   Uses the last element in the App Array to 
%				   Designate the modulation depth for 'lin'
%	seed		: Starting random seed for generating ripple noise or 
%			  dynamic moving ripple parameters (Default=1) - Ignore

function []=ripnoisemin(filename,f1,f2,fRD,fFM,MinRD,MaxRD,MinFM,MaxFM,App,M,Fs,NS,NB,Axis,Block,DF,AmpDist,seed)

%Parameter Conversions
N=32;				%Noise Reconstruction Size  ->> Before changing check ripgensin!!!
LL=1000;			%Noise Upsampling Factor    ->> Before changing check ripgensin!!!
Fsn=Fs/LL;			%Noise Sampling Frequency
Mn=M/LL;			%Noise Signal Length
Mnfft=2^(ceil(log2(Mn)));	%Used for Signal Generation
fphase=5*fRD;			%Temporal Phase Signal Bandlimit Frequency
fphase=2;

%Log Frequency Axis
if Axis=='log'
	XMax=log2(f2/f1);
	X=(0:NS-1)/(NS-1)*XMax;
	faxis=f1*2.^X;
else
	faxis=(0:NS-1)/(NS-1)*(f2-f1)+f1;
	X=log2(faxis/f1);
end

%Generating Ripple Density Signal (Uniformly Distributed Between [MinRD , MaxRD] )
count=1;
if exist('seed')
	seedt=seed;
    disp('ok1')
else
	seedt=1;
end
if Block=='y'
% 	for k=1:NB
% 		for l=1:NB
% 			RD(count,:)=(k-1)*MaxRD/NB ,...
% 			+ MaxRD/NB*noiseunif(fRD,Fsn,Mnfft,seedt);
% 			count = count+1;
%             seedt = seedt+1;
% 		end
% 	end
else
	for k=1:NB
		RD(count,:)=(MaxRD-MinRD)*noiseunif(fRD,Fsn,Mnfft,seedt)+MinRD;
		count = count+1;
        seedt = seedt+1;
	end
end

%Generating the Teporal Modulation Signal 
%Signal Varies Between [-MinFM, -MaxFM] and [MinFM, MaxFM]
%Maximum Bandlimit Freqeuncy is divided by 2
%Beacuse of Clipping of the Parameter FM 
count=1;
if exist('seed')
	seedt=seed;
    disp('ok2')
else
	seedt=1;
end
if Block=='y'
% 	for k=1:NB
% 		for l=1:NB
% 			FM1=noiseunif(fFM/2,Fsn,Mnfft,seedt+NB^2+1);
% 			FM(count,:)=FM1*2*MaxFM/NB + 2*MaxFM/NB*(NB/2-k);
% 			count = count+1;
%             seedt = seedt+1;
% 		end
% 	end
else
	for k=1:NB
        if MinFM == 0
            
            FM1=noiseunif(fFM,Fsn,Mnfft,seedt+NB^2+1);
            FM(count,:)=2*MaxFM*(FM1-.5);
            count=count+1;
            seedt=seedt+1;

        else
%             % only positive
%             FM1=noiseunif(fFM/2,Fsn,Mnfft,seedt+NB^2+1);
%             FM(count,:) = (MaxFM-MinFM)*FM1+MinFM; %the range is from MinFM to MaxFM (no negative).
%             count = count+1;
%             seedt = seedt+1;
            
%             % only negative
%             FM1=noiseunif(fFM/2,Fsn,Mnfft,seedt+NB^2+1);
%             FM(count,:) = (-MaxFM+MinFM)*FM1-MinFM; %the range is from MinFM to MaxFM (no negative).
%             count = count+1;
%             seedt = seedt+1;

% move across positive and negative - cutting the range between -MinFM and MinFM.
            FM1=noiseunif(fFM/2,Fsn,Mnfft*3,seedt+NB^2+1);
            FM1=2*MaxFM*(FM1-.5); % The range is from -MaxFM to MaxFM.
            f=find(FM1>-MinFM & FM1<MinFM); 
            FM1(f) =[]; % eliminate the range from -MinFM to MinFM.
            FM(count,:)=FM1(:,1:Mnfft);
            count = count+1;
            seedt = seedt+1;
            
%             FM1=noiseunif(fFM/2,Fsn,Mnfft,seedt+NB^2+1);
%             FM1 = (MaxFM-MinFM)*FM1+MinFM; % range from MinFM to MaxFM
%             seedt = seedt+1;
%             FM2=noiseunif(fFM/2,Fsn,Mnfft,seedt+NB^2+1);
%             FM2 = (-MaxFM+MinFM)*FM2-MinFM; % range from -MinFM to -MaxFM
%             seedt = seedt+1;
%             FMs=[FM1 FM2];
%             FMs = FMs(randperm(length(FMs))); % randomize the order of FMs
%             FM(count,:)=FMs(1:length(FM1)); %pick up randomly for the half length. Then, the range is from MinFM to MaxFM and from -MinFM to -MaxFM.
%             count = count+1;
            
        end
	end
end

%Genreating Ripple Phase Components
%Note that: RP = 2*pi/Fsn*intfft(FM);
%
if Block=='y'
	for k=1:NB*NB
		RP(k,:)=2*pi/Fsn*intfft(FM(k,:));
	end
else
	for k=1:NB
		RP(k,:)=2*pi/Fsn*intfft(FM(k,:));
	end
end

%Initial Carrier Sinusoid Phases 
phase=2*pi*rand(1,NS);

%For Displaying Statistics
% figure(1),hist(RD,100),title('Ripple Density Distribution')
% figure(2),hist(FM,100),title('Modulation Frequency Distribution')
% figure(3),plot((1:length(FM))/Fsn,FM),title('Modulation Frequency')
% figure(4),plot((1:length(RD))/Fsn,RD),title('Ripple Density')
%save data RD FM RP Fsn

%Generating Ripple Noise
K=0;
flag=0;		%Marks Last segment
for k=2:N:Mn-N-1

	%Extracting RD and RP Noise Segment 
	RDk    = RD(:,k-1:k+N-1);
	RPk    = RP(:,k-1:k+N-1);

	%Interpolating RD and RP
	MM=length(RDk)-1;
	RDint   = interp10(RDk,3);
	RPint   = interp10(RPk,3);
	RDint   = RDint(:,1:N*LL);
	RPint   = RPint(:,1:N*LL);

	%Generating Ripple Noise
	[Y,phase,SpecProf,faxis,taxis]=noisegensinb(f1,f2,RDint,RPint,App,Fs,phase,fphase,K,NB,MaxRD,MaxFM,Axis,DF,AmpDist);
	clear RDint RPint

	%Saving Sound Files
	if strcmp(AmpDist,'lin')

		%Saving Lin AmpDist
		for i=1:length(App)
			tofloat([filename int2str(App(i)) 'Lin.bin'],Y(i,:));
		end

	elseif strcmp(AmpDist,'dB')

		%Saving dB AmpDist
		for i=1:length(App)
			tofloat([filename int2str(App(i)) 'dB.bin'],Y(i,:));
		end

	else
	
		%Saving dB AmpDist
		for i=1:length(App)-1
			tofloat([filename int2str(App(i)) 'dB.bin'],Y(i,:));
		end

		%Saving Lin AmpDist
		tofloat([filename int2str(App(i+1)) 'Lin.bin'],Y(i+1,:));

	end
	clear Y

	%Writing Spectral Profile File as 'float' file
	NT=length(taxis);
	NF=length(faxis);
	tofloat([filename '.spr'],reshape(SpecProf,1,NT*NF));
	clear SpecProf 

	%Updating Display
	K=K+1;
	clc
	disp(['Segment ' num2str(K) ' Done'])

	%Saving Parameters Just in Case it Does Not Reach End
	if K==1
		%Saving Parameters
		f=['save ' filename '_param'];
	end

end

%Saving Parameters
clear RPint RDint RDk RPk k K Y flag f count FM1 filenumber MM ans l;
v=version;
if strcmp(v(1),'5')
	f=['save ' filename '_param -v4'];
	eval(f)
else
	f=['save ' filename '_param'];
	eval(f)
end

%Closing All Opened Files
fclose('all');
