function [output, U] = mc_create_ppi_regressors(SPM,Y,session,deconvolve);
% Will return task and PPI regressors for each condition in the model
% Can do deconvolution (a la SPM) or convolve and then interact (a la FSL)

%need to do this without an SPM.mat file
%need TR, T, T0
%onsets/durations/names
%still requires an SPM file if using deconvolution

%if you aren't using an SPM file, the format of the SPM argument is:
%SPM{1} = [TR, T, T0];
%SPM{2}.name = {'cond1','cond2',...}
%SPM{2}.ons = {[onsets for 1],[onsets for 2],...}
%SPM{2}.dur = {[dur for 1],[dur for 2],...}

if (deconvolve)
    if (~isstruct(SPM))
        error('PPI with deconvolution requires an SPM structure as the first argument');
    end
end

if (iscell(SPM))
    orig = SPM;
    clear SPM
    SPM.xY.RT = orig{1}(1);
    SPM.xBF.dt = SPM.xY.RT/orig{1}(2);
    SPM.xBF.T0 = orig{1}(3);
    %TR hardcoded to 1 for unit='sec'
    TR = 1;
    dt = SPM.xBF.dt;
    T = orig{1}(2);
    k = numel(Y);
    for i = 1:length(orig{2}.name)
        SPM.Sess(session).U(i).name{1} = orig{2}.name{i};
        ons = orig{2}.ons{i};
        dur = orig{2}.dur{i};
        if (size(ons,1)==1)
            ons = ons';
        end
        if (size(dur,1)==1)
            dur = dur';
        end
        
        u     = ons.^0;
        
        ton       = round(ons*TR/dt) + 33;               % onsets
        tof       = round(dur*TR/dt) + ton + 1;          % offset
        sf        = sparse((k*T + 128),size(u,2));
        ton       = max(ton,1);
        tof       = max(tof,1);
        for j = 1:length(ton)
            if size(sf,1) > ton(j)
                sf(ton(j),:) = sf(ton(j),:) + u(j,:);
            end
            if size(sf,1) > tof(j)
                sf(tof(j),:) = sf(tof(j),:) - u(j,:);
            end
        end
        sf        = cumsum(sf);                         % integrate
        sf        = sf(1:(k*T + 32),:);                 % stimulus
        
        SPM.Sess(session).U(i).u = sf;
    end
    
end

Sess = SPM.Sess(session);

U.name = {};
U.u = [];
U.w = [];

u = length(Sess.U);
for i = 1:u
    for j = 1:length(Sess.U(i).name)
        U.w = [U.w 1];
        U.u = [U.u Sess.U(i).u(33:end,j)];
        U.name{end+1} = Sess.U(i).name{j};
    end
end

% Setup variables
%--------------------------------------------------------------------------
RT      = SPM.xY.RT;
dt      = SPM.xBF.dt;
NT      = round(RT/dt);
fMRI_T0 = SPM.xBF.T0;
%N       = length(xY(1).u);
N       = numel(Y);
k       = 1:NT:N*NT;                       % microtime to scan time indices

output = zeros(N,size(U.w,2)*2 + 1);

% Create basis functions and hrf in scan time and microtime
%--------------------------------------------------------------------------
hrf = spm_hrf(dt);

if (deconvolve)
    
% Create convolved explanatory {Hxb} variables in scan time
%--------------------------------------------------------------------------
    xb  = spm_dctmtx(N*NT + 128,N);
    Hxb = zeros(N,N);
    for i = 1:N
        Hx       = conv(xb(:,i),hrf);
        Hxb(:,i) = Hx(k + 128);
    end
    xb = xb(129:end,:);
    
% Get confounds (in scan time) and constant term
%--------------------------------------------------------------------------
%X0 = xY(1).X0;
    X0 = [];
    M  = size(X0,2);
    
% Remove confounds and save Y in ouput structure
%--------------------------------------------------------------------------
%Yc    = Y - X0*inv(X0'*X0)*X0'*Y;
    
% Specify covariance components; assume neuronal response is white
% treating confounds as fixed effects
%--------------------------------------------------------------------------
    Q = speye(N,N)*N/trace(Hxb'*Hxb);
    Q = blkdiag(Q, speye(M,M)*1e6  );
    
% Get whitening matrix (NB: confounds have already been whitened)
%--------------------------------------------------------------------------
    W = SPM.xX.W(Sess.row,Sess.row);
    
% Create structure for spm_PEB
%--------------------------------------------------------------------------
    clear P
    P{1}.X = [W*Hxb X0];        % Design matrix for lowest level
    P{1}.C = speye(N,N)/4;      % i.i.d assumptions
    P{2}.X = sparse(N + M,1);   % Design matrix for parameters (0's)
    P{2}.C = Q;
    
    C  = spm_PEB(Y,P);
    xn = xb*C{2}.E(1:N);
    xn = spm_detrend(xn);

else
    xn = spm_detrend(interp(Y,NT));
end

% Setup psychological variable from inputs and contrast weights
%----------------------------------------------------------------------
PSY = zeros(N*NT,1);
output(:,1) = Y;
idx = 2;
for i = 1:size(U.u,2)
    %PSY = PSY + full(U.u(:,i) * U.w(i));
    PSY = full(U.u(:,i) * U.w(i));
    
    %end
    %PSY = spm_detrend(PSY);
    if (~deconvolve)
        PSY = conv(PSY,hrf);
        PSY = PSY(1:numel(xn));
    end
    
% Multiply psychological variable by neural signal
%----------------------------------------------------------------------
    PSYxn = PSY.*xn;
    
% Convolve, convert to scan time, and account for slice timing shift
%----------------------------------------------------------------------
    if (deconvolve)
        ppi = conv(PSYxn,hrf);
    else
        ppi = PSYxn;
    end
        
    ppi = ppi((k-1) + fMRI_T0);
    
% Convolve psych effect, convert to scan time, and account for slice
% timing shift
%----------------------------------------------------------------------
    if (deconvolve)
        PSYHRF = conv(PSY,hrf);
    else
        PSYHRF = PSY;
    end
    
    PSYHRF = PSYHRF((k-1) + fMRI_T0);
    
    output(:,idx) = PSYHRF;
    idx = idx + 1;
    output(:,idx) = ppi;
    idx = idx + 1;
end

    
