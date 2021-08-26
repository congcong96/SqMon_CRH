function params = strf_parameters(strf, trigger)
%strf_parameters - get strf parameters from significant portions
%   of the strf and the rtf (ripple transfer function.
%
% Calculates temporal, spectral modulation transfer
% functions as well as the separability index for all strfs.
%
% params = strf_parameters(strf, trigger)
% 
% You can use the function plot_strf_parameters.m to
% view the results of the analysis.
%
% The output 'params' is a struct array having length = length(strf)
%
% params has the following form:
%
%    params.exp -> experiment
%    params.site -> recording site
%    params.chan -> channel of neuron
%    params.depth -> probe depth
%    params.position -> neuron depth
%    params.stim -> dmr1, dmr2, ripple noise, etc
%    params.atten -> attenuation of sound
%    params.sm -> max spectral modulation in stimulus
%    params.tm -> max temporal modulation in stimulus
%    params.mdb -> modulation depth of stimulus
%    params.spl -> sound level of stimulus
%    params.n0 -> number of spikes
%    params.w0 -> firing rate
%    params.percent_energy -> different power levels for STRF analysis
%    params.rfenergy -> energy in STRF = root mean square (RMS)
%    params.tmf -> temporal modulation frequency axis
%    params.xmf -> spectral modulation frequency axis
%    params.rtf -> ripple transfer function
%    params.singvals -> singular vals from STRF decomposition
%    params.eigvals -> eigenvalues for decomposition
%    params.tci -> temporal correlation index
%    params.sci -> spectral correlation index
%    params.pli -> phase-locking index
%
%   caa 10/20/02


% parameters to obtain:
% separability index, i.e. the singular and eigen- values
% temporal correlation index
% spectral correlation index
% asymmetry index
% phase-locking index
% best tm
% best sm
% mtfs: tm and sm
% strf energy - use function strf_energy.m
% firing rate
% number of spikes
% duration
%
% and later, feature selectivity index


params = struct(...
    'exp',            [], ...
    'site',           [], ...
    'chan',           [], ...
    'depth',          [], ...
    'position',       [], ...
    'stim',           [], ...
    'atten',          [], ...
    'sm',             [], ...
    'tm',             [], ...
    'mdb',            [], ...
    'spl',            [], ...
    'n0',             [], ...
    'w0',             [], ...
    'percent_energy', [], ...
    'rfenergy',       [], ...
    'tmf',            [], ...
    'xmf',            [], ...
    'rtf',            [], ...
    'singvals',       [], ...
    'eigvals',        [], ...
    'tci',            [], ...
    'sci',            [], ...
    'pli',            []);


percent_energy = [80 85 90 95 97.5 100];
energy = zeros(size(percent_energy));

for i = 1:length(strf)
    
    % Define some initial parameters
    p = 0.001;
    n0 = strf(i).n0contra; %number of spikes
    mdb = strf(i).mdb; % the modulation depth of the ripple stimulus
    w0 = strf(i).w0contra;
    sm = strf(i).sm;
    tm = strf(i).tm;
    
    stim = strf(i).stim; % stimulus type
    fs = strf(i).fs;
    t = strf(i).taxis; %time axis in strf plot
    x = log2(strf(i).faxis ./ min(strf(i).faxis)); %frequency in log2 (500Hz - 40kHz)
    dur = (trigger(end)-trigger(1)) / fs; % duration of the entire stimulus in seconds
    
    % Get the significant strf
    [rfsig] = significant_strf(strf(i).rfcontra, p, n0, mdb, dur);
    
    % take the singular value decomposition to decompose
    % the strf into separable subunits
    [u,s,v] = svd(rfsig); % singluar values arrange in a descending way
    singvals = sum(s); %sindular values (sum of columns as s is a diagnal matrix)
    eigvals = singvals .^ 2;%eigen values of strf, singular value = sqrt(eigenvalue)
    
    % Get phase-locking index
    pli = phase_locking_index(rfsig, mdb, w0, stim);
    
    % Get temporal correlation index
    dt = t(2)-t(1);
    [rt, tshift] = temporal_correlation_index(rfsig, dt);
    
    % Get spectral correlation index
    dx = x(2)-x(1);
    [rx, xshift] = spectral_correlation_index(rfsig, dx);
    
    % get strf energy at different power levels
    cumpower = 100 * cumsum(eigvals) / sum(eigvals+eps);
    
    rfsvd = zeros(size(rfsig,1), size(rfsig,2), length(percent_energy));
    
    for j = 1:length(percent_energy)
        
        % Get rid of noise for plotting appearance
        index_percent_power = find(cumpower <= percent_energy(j)); 
        %idex of eigen values with cumulative power smaller than the enegy level
        num_sing_vals = min([max(index_percent_power) size(u,2) length(singvals)]);
        %max(index_percent_power) shows how many eigenvalues are needed to
        %obtain the required energy level
        
        rftemp = zeros(size(rfsig));
        for k = 1:num_sing_vals
            rftemp = rftemp + singvals(k) * u(:,k) * v(:,k)';
        end % (for k)
        
        % Get STRF energy
        [energy(j)] = strf_energy(rftemp, mdb, stim);
        
        % Save RFs at various SVD power levels
        rfsvd(:,:,j) = rftemp;
       
    end % (for j)
    
    maxtmf = ceil( 1 / dt );
    maxsmf = ceil( 1 / dx );
    
    ntbins = maxtmf; % this will give us 1 Hz tmtf resolution
    nfbins = ceil(maxsmf / 0.15); % make resolution 0.15 cycles / octave
    
    dtfreq = maxtmf / ntbins; % fft temporal frequency resolution, 1 Hz
    dxfreq = maxsmf / nfbins; % fft spectral frequency resolution, 0.15 cycle/octove
   
    for k = 1:size(rfsvd,3)
        % get the ripple transfer function/s
        rfk = fft2(rfsvd(:,:,k), nfbins, ntbins);
        rtftemp = fftshift(abs(rfk));
        
        % Get the tm frequency vector - will be in cycles per second
        if ( mod(size(rtftemp,2),2) )
            tmf = [-(size(rtftemp,2)-1)/2:(size(rtftemp,2)-1)/2]*dtfreq;
        else
            tmf = [-size(rtftemp,2)/2:(size(rtftemp,2)/2-1)]*dtfreq;
        end
        
        itmf0 = find(tmf==0);
        itmfp40 = find(tmf<=40);
        ditmf = max(itmfp40)-itmf0;
        itmf = [itmf0-ditmf:itmf0+ditmf];
        tmf = tmf(itmf);
        
        % Get the sm frequency vector - will be in cycles per octave
        if ( mod(size(rtftemp,1),2) )
            xmf = [-(size(rtftemp,1)-1)/2:(size(rtftemp,1)-1)/2]*dxfreq;
        else
            xmf = [-size(rtftemp,1)/2:(size(rtftemp,1)/2-1)]*dxfreq;
        end
        
        ixmf0 = find(xmf == 0);
        ixmf4 = find(xmf <= 4);
        ixmf = ixmf0:max(ixmf4);
        xmf = xmf(ixmf);
        
        rtf(:,:,k) = rtftemp(ixmf, itmf);
        
        if sum(sum(rtf(:,:,k))) ~= 0
            
            %             % added by Natsumi 2Oct17 to calculate tBMF, sBMF etc
            %             [RTFparams] = rtf_parameters(rtf(:,:,k),xmf,tmf);
            %
            %             % Params --- from rtf_parameters
            %             params(i).tMTF(:,:,k) = RTFparams.tMTF;
            %             params(i).sMTF(:,:,k) = RTFparams.sMTF;
            %             % params(i).Fm = RTFparams.Fm;
            %             % params(i).RD = RTFparams.RD;
            %             % % params(i).tMTF = tMTFu;
            %             % % params(i).sMTF = sMTFu;
            %             % % params(i).Fm = Fmu;
            %             % % params(i).RD = RDu;
            %             params(i).tBMF(k) = RTFparams.tBMF;
            %             params(i).sBMF(k) = RTFparams.sBMF;
            %             params(i).tMTF_BW(k) = RTFparams.tMTF_BW;
            %             params(i).sMTF_BW(k) = RTFparams.sMTF_BW;
            %             params(i).tMTF_BP(k) = RTFparams.tMTF_BP;
            %             params(i).sMTF_BP(k) = RTFparams.sMTF_BP;
            %             params(i).tMTF_Max(k) = RTFparams.tMTF_Max; % added 20Jul2017 Natsumi
            %             params(i).sMTF_Max(k) = RTFparams.sMTF_Max; % added 20Jul2017 Natsumi
            %             params(i).tMTF_Min(k) = RTFparams.tMTF_Min; % added 2Oct2017 Natsumi
            %             params(i).sMTF_Min(k) = RTFparams.sMTF_Min; % added 2Oct2017 Natsumi
            
            % updated to rtf_parameters02 then to  rtf_parameters03
            [RTFparams] = rtf_parameters04(rtf(:,:,k),xmf,tmf);

            % Params --- from rtf_parameters
            params(i).tMTF(:,:,k) = RTFparams.tMTF;
            params(i).sMTF(:,:,k) = RTFparams.sMTF;
            params(i).Fm(:,:,k) = RTFparams.Fm;
            params(i).RD(:,:,k) = RTFparams.RD;
            
            params(i).tBMF(k) = RTFparams.tBMF;
            params(i).sBMF(k) = RTFparams.sBMF;
            params(i).tMTF_BW(k) = RTFparams.tMTF_BW;
            params(i).sMTF_BW(k) = RTFparams.sMTF_BW;
            params(i).tMTF_Min(k) = RTFparams.tMTF_Min; % added 2Oct2017 Natsumi
            params(i).sMTF_Min(k) = RTFparams.sMTF_Min; % added 2Oct2017 Natsumi]
            params(i).tMTF_Max(k) = RTFparams.tMTF_Max; % added 20Jul2017 Natsumi
            params(i).sMTF_Max(k) = RTFparams.sMTF_Max; % added 20Jul2017 Natsumi
            params(i).tMTF_BP(k) = RTFparams.tMTF_BP;
            params(i).sMTF_BP(k) = RTFparams.sMTF_BP;

            % added for rtf_parameters02
            params(i).tMTFu(:,:,k) = RTFparams.tMTFu;
            params(i).sMTFu(:,:,k) = RTFparams.sMTFu;
            params(i).Fmu(:,:,k) = RTFparams.Fmu;
            params(i).RDu(:,:,k) = RTFparams.RDu;
            
            params(i).tMTF_LP(k) = RTFparams.tMTF_LP;
            params(i).sMTF_LP(k) = RTFparams.sMTF_LP;
            params(i).tMTF_HP(k) = RTFparams.tMTF_HP;
            params(i).sMTF_HP(k) = RTFparams.sMTF_HP;
            
            params(i).tBMF1(k) = RTFparams.tBMF1;
            params(i).sBMF1(k) = RTFparams.sBMF1;
            params(i).tMTF_BW1(k) = RTFparams.tMTF_BW1;
            params(i).sMTF_BW1(k) = RTFparams.sMTF_BW1;
            params(i).tMTF_Min1(k) = RTFparams.tMTF_Min1;
            params(i).sMTF_Min1(k) = RTFparams.sMTF_Min1; 
            params(i).tMTF_Max1(k) = RTFparams.tMTF_Max1; 
            params(i).sMTF_Max1(k) = RTFparams.sMTF_Max1;
            
            params(i).tBMF2(k) = RTFparams.tBMF2;
            params(i).sBMF2(k) = RTFparams.sBMF2;
            params(i).tMTF_BW2(k) = RTFparams.tMTF_BW2;
            params(i).sMTF_BW2(k) = RTFparams.sMTF_BW2;
            params(i).tMTF_Min2(k) = RTFparams.tMTF_Min2;
            params(i).sMTF_Min2(k) = RTFparams.sMTF_Min2;
            params(i).tMTF_Max2(k) = RTFparams.tMTF_Max2;
            params(i).sMTF_Max2(k) = RTFparams.sMTF_Max2;
            
            params(i).Max_tMTFu12(k) = RTFparams.Max_tMTFu12;
            params(i).Max_sMTFu12(k) = RTFparams.Max_sMTFu12;
            
        end
        
    end % (for k)
    
    
    % assign parameters to output struct array
    params(i).exp = strf(i).exp;
    params(i).site = strf(i).site;
    params(i).chan = strf(i).chan;
    params(i).depth = strf(i).depth;
    params(i).position = strf(i).position;
    params(i).stim = strf(i).stim;
    params(i).atten = strf(i).atten;
    params(i).sm = strf(i).sm;
    params(i).tm = strf(i).tm;
    params(i).mdb = strf(i).mdb;
    params(i).spl = strf(i).spl;
    params(i).n0 = n0;
    params(i).w0 = w0;
    params(i).percent_energy = percent_energy(:);
    params(i).rfenergy = energy(:);
    params(i).tmf = tmf;
    params(i).xmf = xmf;
    params(i).rtf = rtf;
    params(i).singvals = singvals;
    params(i).eigvals = eigvals;
    params(i).tci = [tshift(:) rt(:)];
    params(i).sci = [xshift(:) rx(:)];
    params(i).pli = pli;
    
    params(i).strfenergy = rfsvd(:);
end % (for i)


