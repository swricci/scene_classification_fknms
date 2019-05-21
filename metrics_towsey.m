
%%2. Signal aquisition and pre-processing

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
    

%%3. Noise removal
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

%noise reduction wave envelope (uses z_env)
%1. find max and min of waveform (each ROW):
wmin = min(z_env_db,[],2);
wmax = max(z_env_db,[],2);

%2. Compute 100 bin histogram of dB values
whist = histogram(z_env_db(1,:),100);
whistcount = whist.BinCounts;

[N,edges,bin] = histcounts(z_env_db(1,:),100);

%3. Smooth the histogram using a moving average filter
windowSize = 5; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
smooth_hist = filter(b,a,whistcount);

%Mode of histogram
[mhist,ind] = max(whistcount);
binedges = whist.BinEdges(ind:ind + 1)




