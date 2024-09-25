function [Y,fdfilt] = mc_dcan_bold_preprocessing(X,TR,Trim,ConfoundsFile,FDthresh,lp_Hz,hp_Hz,bp_order)

if (strcmp(class(X),'single'))
    X = double(X);
end

[r,c] = size(X);
Rr = zeros(r,c);

[path,file,ext] = fileparts(ConfoundsFile);
dat = readtable(fullfile(path,[file ext]),'Delimiter','\t','FileType','text','TreatAsEmpty','n/a');
f = fields(dat);
idxcsf = contains(f,'csf') & ~contains(f,'derivative') & ~contains(f,'power2');
idxwm = contains(f,'white_matter') & ~contains(f,'derivative') & ~contains(f,'power2');
idxgs = contains(f,'global') & ~contains(f,'derivative') & ~contains(f,'power2');

wm = table2array(dat(:,idxwm));
vent = table2array(dat(:,idxcsf));
gs = table2array(dat(:,idxgs));

idxoT = contains(f,'trans_') & ~contains(f,'derivative') & ~contains(f,'power2');
idxoR = contains(f,'rot_') & ~contains(f,'derivative') & ~contains(f,'power2');

mp = table2array([dat(:,idxoT) dat(:,idxoR)]);
tmp = [repmat(0,1,6);diff(mp)];
tmp(:,4:6) = tmp(:,4:6)*50;
fdunfilt = sum(abs(tmp),2);

fmp = mc_notch_filter(mp,0.31,0.43,0.8);
tmp = [repmat(0,1,6);diff(fmp)];
tmp(:,4:6) = tmp(:,4:6)*50;
fdfilt = sum(abs(tmp),2);

fd = fdfilt;

d = [repmat(0,1,6);diff(mp)];
q = [mp.^2 d.^2];
FR = [mp d q];

th = fd<=FDthresh;
th(1:Trim) = 0;
th = th==1;

R = [wm vent gs FR];
mR = mean(R(th,:));
R = R-mR;
n_reg = size(R,2);
bR = zeros(n_reg,2);
x = 1:r;
local_x = x(th);
for i = 1:n_reg
    bR(i,:) = polyfit(local_x,R(th,i)',1);
    R(:,i) = R(:,i) - polyval(bR(i,:),x)';
end

mX = mean(X(th,:));
Xd = X-mX;
bXd = zeros(c,2);
for i = 1:c
    bXd(i,:) = polyfit(local_x,Xd(th,i),1);
    Xd(:,i) = Xd(:,i) - polyval(bXd(i,:),x)';
end

ts1 = Xd(th,:);
ts2 = R(th,:);

%ugh, fix this to just do the regression all at once...
%for i = 1:c
%    b = regress(ts1(:,i),ts2);
%    Rr(:,i) = Xd(:,i) - R*b;
%end
desmtx = [ones(sum(th),1) ts2];
desmtx2 = [ones(numel(th),1) R];

b = pinv(desmtx'*desmtx)*desmtx'*ts1;
Rr = Xd - desmtx2*b;

%interpolation
ix_ol = find(th==0);
ix_in = find(th==1);

int_method = 'linear';
if (numel(ix_in)<=10)
    warning('Too few good frames for interpolation, skipping interpolation');
    int_method = 'none';
end

Rr_int = Rr;

switch int_method
    case 'none'
    case 'linear'
        temp = interp1(ix_in,Rr(ix_in,:),ix_ol);
        Rr_int(ix_ol,:) = temp;
        mean_Rr = mean(Rr(ix_in,:),1);
        if ~th(end)
            Rr_int(end,:) = mean_Rr;
        end
        
        find_nans= mean(temp,2);
        find_nans = find(isnan(find_nans));
        Rr_int(ix_ol(find_nans),:) = repmat(mean_Rr,numel(find_nans),1);
        
    otherwise
end

%bandpass filter
F = 1/TR;
Ny = F/2;
BW_Hz = [lp_Hz hp_Hz];
BW_N = BW_Hz/Ny;
[b,a] = butter(ceil(bp_order/2),BW_N);
padding = zeros(size(Rr_int));
pad_amt = size(padding,1);
temp = cat(1,padding,Rr_int,padding);
Y = filtfilt(b,a,temp);
Y = Y((pad_amt+1):(end-pad_amt),:);





%need
%TR
%trim frames at beginning
%fmriprep path
%fd threshold
%lowpass, highpass, order
%notch filter 

%load WM/CSF/GS regressors from fmriprep confounds file

%load motion paramters from fmriprep confounds file

%calculate filtered motion and filtered fd

%demean, detrend (linear)

%interpolate across filtered fd>0.3 frames least square spectral estimation

%remove nuisance (whole brain, csf, wm, derivs, 24 motion?

%bandpass filter 0.008 - 0.1
%dcan bold 0.009 0.08? filtfilt 2nd order butterworth, no derivs of
%GSR/etc, 24, interpolate fd>0.3, 
%motion filtered fd 0.2 at least 10 minutes data
%concatenate runs, censor at filtered fd>0.08 include if >8 mins (600 frames)


