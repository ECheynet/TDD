function [phi,zeta,eta] = TDD(Az,fs,fn,fnMin,fnMax,varargin)
% [phi,zeta,nu] = TDD(Az,fs,fn,fnMin,fnMax,varargin) estimates the mode
% shapes phi, modal displacements nu and modal damping ratios zeta of a 
% structures using the time Domain Decomposition (TDD) mthod [1]
%
%% Inputs
%  * Az: acceleration data. Matrix of size [Nyy x N] where Nyy
% is the number of sensors, and N is the number of time steps
%  * fs: scalar [1 x 1] sampling frequencies
%  * fn: Vector of size [1x M];  "target eigen frequencies"
%  * fnMin: Vector of size [1x M]; Lower boundary for the band-pass
%  filtering of fn.
%  * fnMax: Vector of size [1x M]; Upper boundary for the band-pass
%  filtering of fn.
% Optional inputs as inputParser (varargin):
%  * ModeNormalization: [1 x 1 ]: 1 or 0: option for mode normalization
%  * dataPlot: 1 or 0: option for intermediate data plot (e.g. checking procedure)
%  * Ts: [1 x 1]: float: option for duration of autocorrelation function (for estimation of damping ratio only)
%% Outputs
% * phi:   Vector of size [1x Nyy]; Estimated mode shapes
% * zeta:  Vector of size [1x M]; Estimated modal damping ratios
% * nu: matrix of size [M x N]; Estimated modal displacements
% 
%% Syntax
% [phi,zeta,nu] = TDD(Az,fs,fn,fnMin,fnMax,fn)
% [phi,zeta,nu] = TDD(Az,fs,fn,fnMin,fnMax,'ModeNormalization',0)
% [phi,zeta,nu] = TDD(Az,fs,fn,fnMin,fnMax,'TS',25)
%
%% Author Infos
% E. Cheynet - UiS - last modified 20/19/2018
%
%% Reference
% [1] Byeong Hwa Kim, Norris Stubbs, Taehyo Park, A new method to extract
% modal parameters using output-only responses, Journal of Sound and Vibration,
% Volume 282, Issues 1–2, 6 April 2005, Pages 215-230, ISSN 0022-460X,
% http://dx.doi.org/10.1016/j.jsv.2004.02.026.

%% Preprocessing
Nmodes = numel(fn);
[Nyy,N] = size(Az); %Dimension of displacement matrix.
dt = 1/fs;

p = inputParser();
p.CaseSensitive = false;
p.addOptional('ModeNormalization',1) % option for mode normalization (1 or 0)
p.addOptional('Ts',30) % option for duration of autocorrelation function (for estimation of damping ratio only)
p.parse(varargin{:});
ModeNormalization = p.Results.ModeNormalization;
Ts = p.Results.Ts;

%% Band pass filtering for each modes
Az_filt = zeros(Nmodes,Nyy,N);
for ii=1:Nmodes
    h1=fdesign.bandpass('N,F3dB1,F3dB2',8,fnMin(ii),fnMax(ii),fs);
    d1 = design(h1,'butter');
    for jj=1:Nyy
        Az_filt(ii,jj,:) = filtfilt(d1.sosMatrix,d1.ScaleValues, Az(jj,:));
    end
end

%% Mode shapes (phi), modal displacement (nu) and modal damping ratio (zeta)
phi = zeros(Nmodes,Nyy);
eta = zeros(Nmodes,N);
zeta = zeros(size(fn));
for ii=1:Nmodes,
    G = cov(squeeze(Az_filt(ii,:,:))'); % covariance matrix
    [U,~,~] = svd(G); % svd of covariance matrix
    % mode shapes
    phi(ii,:) = U(:,1)';
    % modal displacement
    eta(ii,:) = (phi(ii,:)*phi(ii,:)')\phi(ii,:)*squeeze(Az_filt(ii,:,:));
    % modal damping ratio
    [Y,t] = NExT(eta(ii,:),dt,Ts,1);% cross-covariance calculated with ifft
    % get the envelop of the curve with the hilbert transform:
    envelop = abs(hilbert(Y));    envelop(1)=Y(1);
    % fit an exponential decay to the envelop
    wn = 2*pi*fn(ii); % -> obtained with peak picking method (fast way)
    [zeta(ii)] = expoFit(envelop,t,wn);
end

% mode normalization
if ModeNormalization==1,
    for ii=1:Nmodes,
        phi(ii,:) = phi(ii,:)./max(abs(phi(ii,:)));% normalisation of the modes
    end
end

    function [IRF,t] = NExT(y,dt,Ts,method)
        % [IRF] = NExT(y,ys,T,dt) implements the Natural Excitation Technique to
        % retrieve the Impulse Response FUnction (IRF) from the cross-correlation
        % of the measured output y.
        %
        % [IRF] = NExT(y,dt,Ts,1) calculate the IRF with cross-correlation
        % calculated by using the inverse fast fourier transform of the
        % cross-spectral power densities  (method = 1).
        %
        % [IRF] = NExT(y,dt,Ts,2) calculate the IRF with cross-correlation
        % calculated by using the unbiased cross-covariance function (method = 2)
        %
        %
        % y: time series of ambient vibrations: vector of size [1xN]
        % dt : Time step
        % method: 1 or 2 for the computation of cross-correlation functions
        % T: Duration of subsegments (T<dt*(numel(y)-1))
        % IRF: impusle response function
        % t: time vector asociated to IRF
        %%
        if nargin<4, method = 2; end % the fastest method is the default method
        if ~ismatrix(y), error('Error: y must be a vector or a matrix'),end
        [Ny,N1]=size(y);
        if Ny>N1,
            y=y';
            [Ny,N1]=size(y);
        end
        % get the maximal segment length fixed by T
        M1 = round(Ts/dt);
        switch method
            case 1
                clear IRF
                for pp=1:Ny,
                    for qq=1:Ny,
                        y1 = fft(y(pp,:));
                        y2 = fft(y(qq,:));
                        h0 = ifft(y1.*conj(y2));
                        IRF(pp,qq,:) = h0(1:M1);
                    end
                end
                % get time vector t associated to the IRF
                t = linspace(0,dt.*(size(IRF,3)-1),size(IRF,3));
                if Ny==1,
                    IRF = squeeze(IRF)'; % if Ny=1
                end
            case 2
                IRF = zeros(Ny,Ny,M1+1);
                for pp=1:Ny,
                    for qq=1:Ny,
                        [dummy,lag]=xcov(y(pp,:),y(qq,:),M1,'unbiased');
                        IRF(pp,qq,:) = dummy(end-round(numel(dummy)/2)+1:end);
                    end
                end
                if Ny==1,
                    IRF = squeeze(IRF)'; % if Ny=1
                end
                % get time vector t associated to the IRF
                t = dt.*lag(end-round(numel(lag)/2)+1:end);
        end
        % normalize the IRF
        if Ny==1,
            IRF = IRF./IRF(1);
        else
        end
    end
    function [zeta] = expoFit(y,t,wn)
        % [zeta] = expoFit(y,t,wn) returns the damping ratio calcualted by fiting
        % an exponential decay to the envelop of the Impulse Response Function.
        %
        % y: envelop of the IRF: vector of size [1 x N]
        % t: time vector [ 1 x N]
        % wn: target eigen frequencies (rad/Hz) :  [1 x 1]
        % zeta: modal damping ratio:  [1 x 1]
        %  optionPlot: 1 to plot the fitted function, and 0 not to plot it.
        %%
        
        % Initialisation
        guess = [1,1e-2];
        % simple exponentiald ecay function
        myFun = @(a,x) a(1).*exp(-a(2).*x);
        % application of nlinfit function
        assert(license('test','Statistics_Toolbox')==1,'The function expoFit requires Matlab Statistics Toolbox.')
        coeff = nlinfit(t,y,myFun,guess);
        % modal damping ratio:
        zeta = abs(coeff(2)./wn);
    end
end



