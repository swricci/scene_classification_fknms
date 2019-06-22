function [z] = chunk_data(y,fs,chunk_size,frame_size,ovlp)

%reshape y: [fs * chunk_size X # of chunks]
%this chunks the data into smaller segments so each column is the number of
%points in the number of seconds of your chunk size (ex. 60 seconds of data
%* fs = number of points)
bwin = chunk_size * fs; %number of points for each column of matrix (30 seconds of data)
x=buffer(y,bwin,ovlp,'nodelay'); %x=buffer(y,bwin,ovlp,'nodelay');  
if x(end)==0; x(:,end)=[]; end  % if last colum is zero padded delete. 

% %steps to demean this matrix?
% [r]=size((x)); 
% msub=repmat(mean(x,1),r,1); % calculate the mean of each column 
% x=x-msub; % remove the mean of each column

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


end

