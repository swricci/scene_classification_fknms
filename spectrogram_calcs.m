% Spectrogram derived metrics


%% Signal processing
%read in data
load STcalibration
dir2process = readtable('dir2process_cont.csv');
%dir = '/Volumes/P3/FKNMS_soundscape/deployment8/LKSPA_Dep8_335851542_ST21/';
chunk_size = 60; %size of data chunk in seconds
frame_size = 4096;
ovlp = 0;

DirIn =char(dir2process.DirIn(1));
[filelist,ftimes,fend]= mktableSTdir(DirIn);

[y,fstart,fs,metadata]=readST_continuous(char(filelist(4).name),DirIn,STcalib);
%Note: y is demeaned and gain adjusted (part of readST function)

%filter data
%high pass filter 50 Hz
[B,A]=butter(1,(50/(fs/2)),'high');
%[H,F] = freqz(B,A,2048,fs);
%figure; plot(F,abs(H));
yfilt = filtfilt(B,A,y);

%reshape y: [fs * chunk_size X # of chunks]
%this chunks the data into smaller segments so each column is the number of
%points in the number of seconds of your chunk size (ex. 60 seconds of data
%* fs = number of points)
bwin = chunk_size * fs; %number of points for each column of matrix (30 seconds of data)
x=buffer(yfilt,bwin,ovlp,'nodelay'); %x=buffer(y,bwin,ovlp,'nodelay');  
if x(end)==0; x(:,end)=[]; end  % if last colum is zero padded delete. 

%steps to demean this matrix?
[r,Nwin]=size((x)); 
msub=repmat(mean(x,1),r,1); % calculate the mean of each column 
x=x-msub; % remove the mean of each column

%Now divide chunks into frames of frame_size points. This creates a 3D
%matrix (z) where:
%dim = [frame size X number of frames in chunk X number of chunks in original recording]
%ex: for 6 hr recording with frame size = 1024, chunk size = 60 sec
% size (z) = 1024 x 2812 x 359 for ONE recording file
z=nan(frame_size,floor((fs*chunk_size)/frame_size),size(x,2));
for i = 1:size(x,2);
    tmp = buffer(x(:,i),frame_size,ovlp,'nodelay');
    if tmp(end)==0; tmp(:,end)=[]; end
    z(:,:,i) = tmp;
end

%Create amp and db spectrograms for a single file (3d matrices = 3d is each
%minute of the recording
for j = 1:15;%size(z,3);    
    
    [r, Nwin] = size(z(:,:,j));
    
    %calculate spectra
    wo = hamming(r); %hamming window
    zo = (z(:,:,j)).*repmat(wo,1,Nwin); %apply window
    nfft = 2^nextpow2(frame_size);
    Y = fft(zo,nfft,1);
    po=2*abs(Y).^2/(nfft*sum(wo.^2)); % scale for PSD accounting for windowing 
    po=po(1:ceil(nfft/2)+1,:); po(1)=0; % take first 1/2 & zero DC 
    %*** should I be zeroing the first row or just the first position?
    
    [prows,~] = size(po); % # rows in po. 
    m=0:1:prows-1; f=m*fs/nfft; 


%     mspec=mean(po,2);
%     tframes = (0:1:size(po,2)-1) * (frame_size/fs);
%     figure; imagesc(tframes,f,10*log10(po)); axis xy; colormap jet; caxis([45, 90]);
    
    %po is the amplitude spectrogram for one minute of the continuous data
    
    po_amp(:,:,j) = po;
    
    po_db(:,:,j) = 10*log10(po); %decibel spectrogram
    
    po_long = reshape(po_amp, size(po_amp,1),[]);
    
end

%Noise removal from spectrogram
    % Step A
    %now we need to create a histogram of intensity values in each
    %frequency bin for the minute
    hcount = nan(size(po_db,1),100);
    edges = nan(size(po_db,1),101);
    for i = 1:size(po_db,1);
        [hcount(i,:),edges(i,:)] = histcounts(po_db(i,:),100); %output bins here too
    end
    
figure;histogram('BinEdges', edges(9,:),'BinCounts',hcount(9,:))
figure; plot(hcount');
    
%smooth with moving average here? win width = 3
    mavg_win_spec = 5; %can set this
    smooth_hcount_po = nan(size(hcount,1),100);
    for i=1:size(hcount,1);
        smooth_hcount_po(i,:) = smoothdata(hcount(i,:),'movmean',mavg_win_spec);
    end
    %figure; plot(smooth_hcount_po');
    h_mode = nan(1,size(smooth_hcount_po,1));
    h_ind = nan(1,size(smooth_hcount_po,1));
    mode_db_freqbin = nan(1,size(smooth_hcount_po,1));
    for i = 1:size(smooth_hcount_po,1);
        [h_mode(i),h_ind(i)] = max(hcount(i,:));
        if h_ind(i) >= 96;
            mode_db_freqbin(i) = mean(edges(i,95:96));
        else
        mode_db_freqbin(i) = mean(edges(i,h_ind(i):h_ind(i) +1));
        end
    end
    
    figure; plot(f,mode_db_freqbin);
    
    %Step B use noise profile from A
    %smooth noise profile
    smooth_mode_db_freq = smoothdata(mode_db_freqbin,'movmean',mavg_win_spec);
    figure; plot(f,smooth_mode_db_freq);
    
    %subtract from po_db matrix to remove noise
    po_db_nr = po_db - smooth_mode_db_freq';
    
    %negative values to 0
    a = find(po_db_nr < 0);
    po_db_nr(a) = 0;

%     figure; imagesc(tframes,f,po_db_nr); axis xy;
%     figure; imagesc(tframes,f,po_db); axis xy;

%spectral indices: result should be a summary vector for each 1 min (so a
%matrix for each recording

%BGNsp = noise profile calculated in step B (smooth_mode_db_freq)

BGNsp = smooth_mode_db_freq;

%PMNsp = power minus noise (max in each freq bin - BGN val)
db_freq_max = max(po_db,[],2);

PMNsp = db_freq_max - BGNsp;

%ACTsp = Activity (fraction of cells in each noise-reduced freq bin that
%exceeds a threshold (3dB).

threshold = 3; %can set this!
for i = 1:size(po_db_nr,1);
    a = find(po_db_nr(i,:) > threshold);
    ACTsp(i) = length(a)/size(po_db_nr,2);
end

%EVNsp (summary index): Events per minute in each noise-reduced frequency
%bin. Event = each time dB value exceeeds 3db threshold from low to high

po_gr_thresh = po_db_nr;
po_gr_thresh(po_gr_thresh<threshold) = 0;
po_gr_thresh(po_gr_thresh>=threshold) = threshold;

EVNsp = [];

for i = 1:size(z_env_db_nr,1);
    pks = [];
    pks = findpeaks(po_gr_thresh(i,:));
    EVNsp(i) = length(pks)/chunk_size;
end


    
 