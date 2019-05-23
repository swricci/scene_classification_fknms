
%% 2. Signal aquisition and pre-processing

%2a. Signal Processing
%split long recording into one-minute segments of audio
%Towsey data ex:
    %fs: 22050
    %seg dur: 60 sec
    %frame size: 512 pts
    %total frames: 2584
    %frame duration: 23.2 ms
    %(22050 * 60)/512 = 2584 frames
    %512 / 22050 = 0.0232 seconds per frame
    
%2b. Production of wave envelope
%Take maximum absolute value in each frame (length of vector = number of
%frames 2584) then convert to dB
%Towsey set minimum value to -90dB


%2c. Spectrograms
%512 or 1024 points
% 0% or 50% overlap
%Hamming window function to each frame prior to FFT
%Spectra smoothed with moving average window width = 3
%Spectral amplitude values:
    %Fourier coefficient (A)
    %Spectral energy/power (A^2)
    %dB (dB = 20 * log10 (A))
    

%% 3. Noise removal
%Three methods: subtraction of mean, median, or mode.
%Median often used - no assumption of noise model
%Towsey uses modal with adaptive level equalization algorithim
    



load STcalibration
dir2process = readtable('dir2process_cont.csv');
%dir = '/Volumes/P3/FKNMS_soundscape/deployment8/LKSPA_Dep8_335851542_ST21/';
chunk_size = 60; %size of data chunk in seconds
frame_size = 1024;
ovlp = 0;
demean = 0;


DirIn =char(dir2process.DirIn(1));
[filelist,ftimes,fend]= mktableSTdir(DirIn);

[y,fstart,fs,metadata]=readST_continuous(char(filelist(4).name),DirIn,STcalib,demean); %no de-noising 

demean = 1;
[y2,fstart2,fs2,metadata2]=readST_continuous(char(filelist(4).name),DirIn,STcalib,demean); %mean

demean = 2;
[y3,fstart3,fs3,metadata3]=readST_continuous(char(filelist(4).name),DirIn,STcalib,demean); %median


%save y, fstart, fs, metadata for file 4 for testing purposes.
load testdata.mat;

%this chunks the data into smaller segments, will need to be futher
%processed...
bwin = chunk_size * fs; %number of points for each column of matrix (30 seconds of data)
x=buffer(y,bwin,ovlp,'nodelay');   
if x(end)==0; x(:,end)=[]; end  % if last colum is zero padded delete. 
% [r,Nwin]=size((x)); 
% msub=repmat(mean(x,1),r,1); % calculate the mean of each column 
% x=x-msub; % remove the mean of each column


z=nan(frame_size,floor((fs*chunk_size)/frame_size),size(x,2));
for i = 1:size(x,2);
    tmp = buffer(x(:,i),frame_size,ovlp,'nodelay');
    if tmp(end)==0; tmp(:,end)=[]; end
    z(:,:,i) = tmp;
end

%wave envelope
for j = 1:size(z,3);
    for k = 1:size(z,2);
        z_env(j,k) = max(abs(z(:,k,j)));
    end
end

z_env_db = 20*log10(z_env);

%% 3a. noise reduction wave envelope (uses z_env)
%1. find max and min of waveform (each ROW):
wmin = min(z_env_db,[],2);
wmax = max(z_env_db,[],2);

%2. Compute 100 bin histogram of dB values
hcount = nan(size(z_env_db,1),100);
edges = nan(size(z_env_db,1),101);
for i = 1:size(z_env_db,1);
    [hcount(i,:),edges(i,:)] = histcounts(z_env_db(i,:),100); %output bins here too
end

figure;histogram('BinEdges', edges(9,:),'BinCounts',hcount(9,:))
figure; plot(hcount');

%3. Smooth the histogram using a moving average filter
windowSize = 5; %can set this!
b = (1/windowSize)*ones(1,windowSize);
a = 1;
smooth_hcount = nan(size(hcount,1),100);
for i=1:size(hcount,1);
    smooth_hcount(i,:) = filter(b,a,hcount(i,:));
end

figure; plot(smooth_hcount');

%4 Mode of histogram
h_mode = nan(1,size(smooth_hcount,1));
h_ind = nan(1,size(smooth_hcount,1));
mode_db = nan(1,size(smooth_hcount,1));
for i = 1:size(smooth_hcount,1);
    [h_mode(i),h_ind(i)] = max(hcount(i,:));
    mode_db(i) = mean(edges(i,h_ind(i):h_ind(i) +1));
end

%5. std dev of noise distribution or if N = 0 skip to step 7, since the
%background noise value will be = mode_db (step 4)

%7. Subtract mode from each waveform value
z_env_db_nr = z_env_db - mode_db';
a = find(z_env_db_nr < 0);
z_env_db_nr(a) = 0;

%SAVE z_env_db z_env_db_nr mode_db

%% Waveform based indices
%BGN, SNR, EVN, ACT (from Phillips et al 2018)

%BGN (summary index): mode of the energy distribution in the waveform
%envelope
BGN = mode_db;

%SNR (summary index): difference between max db value in db envelope and
%BGN db
SNR = wmax - mode_db';

%ACT (summary index): Fraction of values in noise-reduced db env that
%exceed a threshold (Towsey = 3dB).
threshold = 3; %can set this!
for i = 1:size(z_env_db_nr,1);
    a = find(z_env_db_nr(i,:) > threshold);
    ACT(i) = length(a)/size(z_env_db_nr,2);
end

%EVN (summary index): Events per second, average over noise reduced
%one-minute. starts when db envelope crosses threshold from below to above
%(threshold = 3 dB).

%% 2c. Spectrograms
%512 or 1024 points
% 0% or 50% overlap
%Hamming window function to each frame prior to FFT
%Spectra smoothed with moving average window width = 3
%Spectral amplitude values:
    %Fourier coefficient (A)
    %Spectral energy/power (A^2)
    %dB (dB = 20 * log10 (A))
    
% need to create spectrogram from one - minute segments that is NOT
% demeaned/ noise reduced

% x is matrix for each one minute duration segment
% adapt sound_MSPEC where window = hamming, "y" is not demeaned
win = 1024;
spec = nan(win/(2+1,size(x,2));

for j = 4;%:1;%size(x,2));
    %the three lines below essentially create "z" from above...each level
    %of 3rd dimension is for each one -minute segment.
    xspec = buffer(x(:,j), win, ovlp,'nodelay');
    if xspec(end)==0; xspec(:,end) = []; end
    
    %load z here? then loop through third dimension?
    [r, Nwin] = size(xspec);
    
    %calculate spectra
    wo = hamming(r); %hamming window
    zo = xspec.*repmat(wo,1,Nwin); %apply window
    nfft = 2^nextpow2(win);
    Y = fft(zo,nfft,1);
    po=2*abs(Y).^2/(nfft*sum(wo.^2)); % scale for PSD accounting for windowing 
    po=po(1:ceil(nfft/2)+1,:); po(1)=0; % take first 1/2 & zero DC 
    [prows,~] = size(po); % # rows in po. 
    m=0:1:prows-1; f=m*fs/nfft; 
    %smooth with moving average here? win width = 3
    
    %po matrix (no mean removed) is really load in first two frequency bins
    %less than 100 Hz)
    
    t = fstart + (0:1:size(x,2)-1) * (chunk_size/(24*60*60));
    
    %save for this one-minute segment
    
    
    [spec(:,j),f] = sound_MSPEC
    


%% 3b. Noise removal from SPECTROGRAM








