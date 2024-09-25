function [filteredts]= mc_notch_filter(ts,minf,maxf,TR)

% Sampling frequency (Hz)
Fs = 1 / TR; % sampling frequency

% Center frequency of the notch filter (Hz)
f0 = mean([minf,maxf]);

% Bandwidth of the notch filter (Hz)
BW = abs(maxf-minf);

% Design the notch filter
[b, a] = iirnotch(f0/(Fs/2), BW/(Fs/2));

filteredts = filtfilt(b,a,ts);

