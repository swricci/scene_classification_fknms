function [po_amp,po_db] = feature_spectrogram(y,fs,chunk_size,frame_size, overlap)
% feature_spectrogram Returns amplitude and decibel spectrograms for an
% audio file.
%   Amplitude and db spectrograms are calculated for each data segment (ex
%   60 sec) of the original file. The result is a 3D matrix with the 3rd
%   dimension representing each data segment, rows = freq bins, columns =
%   time bins.

% Inputs
    % y: demeaned, gain calibrated, filtered waveform.
    % fs: sample frequency
    % chunk_size: duration of the audio "scenes" in seconds
    % frame_size: number of points in each frame (1024, 2048, 4096 etc)
    % overlap: percent overlap of data frames in calculations

%% Reshape y: [fs * chunk_size X # of chunks]
%this chunks the data into smaller segments so each column is the number of
%points in the number of seconds of your chunk size (ex. 60 seconds of data
%* fs = number of points)
bwin = chunk_size * fs; %number of points for each column of matrix (30 seconds of data)
x=buffer(y,bwin,overlap,'nodelay'); %x=buffer(y,bwin,overlap,'nodelay');  
if x(end)==0; x(:,end)=[]; end  % if last colum is zero padded delete. 

%steps to demean this matrix?
[r,Nwin]=size((x)); 
msub=repmat(mean(x,1),r,1); % calculate the mean of each column 
x=x-msub; % remove the mean of each column

%% Divide chunks into frames of frame_size points. 
%This creates a 3D matrix (z) where:
    %dim = [frame size X number of frames in chunk X number of chunks in original recording]
    %ex: for 6 hr recording with frame size = 1024, chunk size = 60 sec
    % size (z) = 1024 x 2812 x 359 for ONE recording file
z=nan(frame_size,floor((fs*chunk_size)/frame_size),size(x,2));
for i = 1:size(x,2);
    tmp = buffer(x(:,i),frame_size,overlap,'nodelay');
    if tmp(end)==0; tmp(:,end)=[]; end
    z(:,:,i) = tmp;
end

%% Create amp and db spectrograms for a single file 
% (3d matrices = 3d is each minute of the recording
for j = 1:size(z,3);    
    
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
    
    po_amp(:,:,j) = po;
    
    po_db(:,:,j) = 10*log10(po); %decibel spectrogram
    
end
end

