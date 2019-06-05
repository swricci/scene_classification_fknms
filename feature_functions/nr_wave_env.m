function [mode_db,wmax,env_db_nr] = nr_wave_env(wave_env_db,mavg_win)

% The function nr_wave_env reduces the noise in the wave envelope following
% the method use in Towsey et al.

% Inputs
    % wave_env_db: wave envelope in decibels
    % mavg_win: number of points in moving average window (Towsey uses 5)
    
%Output
    % mode_db: mode of energy in wave envelope (for each frame)
    % wmax: maximum dB value in each frame
    % env_db_nr: noise-reduced wave envelope in decibels
    
%Noise removal from waveform (really, wave envelope: uses wave_env_db)
%% 1. find max and min of waveform (each ROW):
wmin = min(wave_env_db,[],2);
wmax = max(wave_env_db,[],2); %wmax is used in SNR calculation

%% 2. Compute 100 bin histogram of dB values
hcount = nan(size(wave_env_db,1),100);
edges = nan(size(wave_env_db,1),101);
for i = 1:size(wave_env_db,1);
    [hcount(i,:),edges(i,:)] = histcounts(wave_env_db(i,:),100); 
end

%% 3. Smooth the histogram using a moving average filter
smooth_hcount = nan(size(hcount,1),100);
for i=1:size(hcount,1);
    smooth_hcount(i,:) = smoothdata(hcount(i,:),mavg_win);
end

%% 4 Mode of histogram (find bin with max count, mid point of bin = mode)
h_mode = nan(1,size(smooth_hcount,1));
h_ind = nan(1,size(smooth_hcount,1));
mode_db = nan(1,size(smooth_hcount,1));
for i = 1:size(smooth_hcount,1);
    [h_mode(i),h_ind(i)] = max(hcount(i,:));
    mode_db(i) = mean(edges(i,h_ind(i):h_ind(i) +1));
end

%% 5. std dev of noise distribution 
% or if N = 0 skip to step 7, since the background noise value 
% will be = mode_db (step 4)
% Note from Towsey: bgn = modal intensity +/- N * stdev.
% if N = 0, bgn = modal intensity

%% 7. Subtract mode from each waveform value, set any values <0 to 0.
env_db_nr = wave_env_db - mode_db';
a = find(env_db_nr < 0);
env_db_nr(a) = 0;
