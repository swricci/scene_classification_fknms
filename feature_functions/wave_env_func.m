function [wave_env] = wave_env_func(z)
%wave_env_func(y,fs, chunk_size, frame_size, overlap)
% The function wave_env_func returns the wave envelope of an audio file.
% The wave envelope is the maximum absolute value of each frame of a
% wavefile that has been divided into shorter duration "chunks".

% Inputs
    % y: demeaned, gain calibrated waveform. If desired, this can also be
    % the filtered waveform
    % fs: sample frequency
    % chunk_size: duration of the audio "scenes" in seconds
    % frame_size: number of points in each frame (1024, 2048, 4096 etc)
    % overlap: percent overlap of data frames in calculations
    
%Output
    %wave_env: # of chunk_size duration chunks in file x number of frames
    %in chunk
    
% %% Reshape y: [fs * chunk_size X # of chunks]
% %this chunks the data into smaller segments so each column is the number of
% %points in the number of seconds of your chunk size (ex. 60 seconds of data
% %* fs = number of points)
% bwin = chunk_size * fs; %number of points for each column of matrix (30 seconds of data)
% x=buffer(y,bwin,overlap,'nodelay'); %x=buffer(y,bwin,ovlp,'nodelay');  
% if x(end)==0; x(:,end)=[]; end  % if last colum is zero padded delete. 
% 
% %% Divide chunks into frames of frame_size points. 
% %This creates a 3D matrix (z) where:
% %dim = [frame size X number of frames in chunk X number of chunks in original recording]
% %ex: for 6 hr recording with frame size = 1024, chunk size = 60 sec
% % size (z) = 1024 x 2812 x 359 for ONE recording file
% z=nan(frame_size,floor((fs*chunk_size)/frame_size),size(x,2));
% for i = 1:size(x,2);
%     tmp = buffer(x(:,i),frame_size,overlap,'nodelay');
%     if tmp(end)==0; tmp(:,end)=[]; end
%     z(:,:,i) = tmp;
% end

%% Production of wave envelope
%find the maximum absolute value of each column in z. This corresponds to
%the maximum value in each frame of the shorter duration data chunk.
for j = 1:size(z,3);
    for k = 1:size(z,2);
        wave_env(j,k) = max(abs(z(:,k,j)));
    end
end