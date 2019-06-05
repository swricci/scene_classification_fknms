function [BGN, SNR, ACT, EVN] = wave_env_metrics(env_db_nr,mode_db,wmax,threshold,chunk_size)

% The function wave_env_metrics calculates summary indices on the wave
% envelope.

% Inputs
    % env_db_nr: noise-reduced wave envelope in decibels
    % mode_db: mode of each frame in wave envelope
    % wmax: max dB in each frame in wave envelope
    % threshold: threshold (in dB) to use in calculaton of EVN, ACT,
    
%Output
    % BGN: Background Noise (set = mode_db)
    % SNR: Signal to noise ratio (difference between max dB and dB BGN
    % ACT: Activity, fraction of values that exceed threshold (ex 3 dB)
    % EVN: Event per second, number of acoustic events per second (where dB
    % env exceeds threshold)
    
% Waveform based summary indices
% The following all are based on wave envelope in dB (some based on nr
% envelope (z_env_db or z_env_db_nr)
%BGN, SNR, EVN, ACT (from Phillips et al 2018)

%% BGN (summary index): mode of the energy distribution in the wave_env
BGN = mode_db;

%% SNR (summary index): difference between max db value in db envelope and
% BGN db
SNR = wmax - mode_db';

%% ACT (summary index): 
%Fraction of values in noise-reduced db env that excees a threshold (3 dB)
for i = 1:size(env_db_nr,1);
    a = find(env_db_nr(i,:) > threshold);
    ACT(i) = length(a)/size(env_db_nr,2);
end

%% EVN (summary index): Events per second, average over noise reduced env
% starts when db envelope crosses threshold from below to above (threshold = 3 dB).
env_gr_thresh = env_db_nr;
env_gr_thresh(env_gr_thresh<threshold) = 0;
env_gr_thresh(env_gr_thresh>=threshold) = threshold;

EVN = [];

for i = 1:size(env_db_nr,1);
    pks = [];
    pks = findpeaks(env_gr_thresh(i,:));
    EVN(i) = length(pks)/chunk_size;
end
