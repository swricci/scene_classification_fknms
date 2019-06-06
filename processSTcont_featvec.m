%processST_continuous_featurevec

load STcalibration
dir2process = readtable('dir2process_cont.csv');
%dir = '/Volumes/P3/FKNMS_soundscape/deployment8/LKSPA_Dep8_335851542_ST21/';
chunk_size = 60; %size of data chunk in seconds
frame_size = 4096;
ovlp = 0;
BGN_long = [];
SNR_long = [];
ACT_long = [];
EVN_long = [];
ENT_long = [];
tframe_long = [];

DirIn =char(dir2process.DirIn(1));
Site=char(dir2process.Site(1));   
DirOut=char(dir2process.DirOut(1)); 
Deployment=dir2process.Deployment(1);
   
[filelist,ftimes,fend]= mktableSTdir(DirIn);

tic
for file = 2:length(filelist);
    [y,fstart,fs,metadata]=readST_continuous(char(filelist(file).name),DirIn,STcalib);
    
    %filter data
    %high pass filter 50 Hz
    [B,A]=butter(1,(50/(fs/2)),'high');
    %[H,F] = freqz(B,A,2048,fs);
    %figure; plot(F,abs(H));
    yfilt = filtfilt(B,A,y);
  
    %calculate wave envelope
    wave_env = wave_env_func(yfilt,fs,chunk_size,frame_size, ovlp);
    wave_env_db = 20 * log10(wave_env);
    
    %noise removal from wave envelope (use wave_env_db)
    mavg_win = 5;
    [mode_db,wmax,env_db_nr] = nr_wave_env(wave_env_db,mavg_win);
    
    %calculate waveform env based indices
    threshold = 3;
    [BGN(file,:), SNR(file,:), ACT(file,:), EVN(file,:)] = wave_env_metrics(env_db_nr,...
        mode_db,wmax,threshold,chunk_size);
    
    %calculate ENT (temporal entropy)
    for i=1:size(wave_env,1);
        sumXi=sum(wave_env(i,:)); % Sum the envelope 
        At=wave_env(i,:)/sumXi;  % Compute the probability mass function 
        Ht(i)=-(sum(At .* log2(At)))/log2(length(At));
    end
    
    ENT(file,:) = Ht;
    
    %determine start time of each frame
    tframe(file,:) = fstart + (0:1:size(wave_env,1)-1) * (chunk_size/(24*60*60));
    
    % create deployment long index
    BGN_long = cat(2,BGN_long,BGN(file,:));
    SNR_long = cat(2,SNR_long,SNR(file,:));
    ACT_long = cat(2,ACT_long,ACT(file,:));
    EVN_long = cat(2,EVN_long,EVN(file,:));
    ENT_long = cat(2,ENT_long,ENT(file,:));
    tframe_long = cat(2,tframe_long, tframe(file,:));
    
end
    
  toc  
  
eval([Site '_' sprintf('%02.0f',Deployment) '_BGN=BGN;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_SNR=SNR;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_ACT=ACT;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_EVN=EVN;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_ENT=ENT;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_tframe=tframe;']);

  
eval([Site '_' sprintf('%02.0f',Deployment) '_BGN_long=BGN_long;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_SNR_long=SNR_long;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_ACT_long=ACT_long;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_ENT_long=ENT_long;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_EVN_long=EVN_long;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_tframe_long=tframe_long;']);

eval(['save ',DirOut, [Site '_' sprintf('%02.0f',Deployment) '_featvec_results.mat'] ' ' [Site '_' sprintf('%02.0f',Deployment) '_BGN']...
    ' ' [Site '_' sprintf('%02.0f',Deployment) '_SNR'] ' ' [Site '_' sprintf('%02.0f',Deployment) '_ACT']...
    ' ' [Site '_' sprintf('%02.0f',Deployment) '_EVN'] ' ' [Site '_' sprintf('%02.0f',Deployment) '_ENT']...
    ' ' [Site '_' sprintf('%02.0f',Deployment) '_tframe'] ' ' [Site '_' sprintf('%02.0f',Deployment) '_BGN_long']...
    ' ' [Site '_' sprintf('%02.0f',Deployment) '_SNR_long'] ' ' [Site '_' sprintf('%02.0f',Deployment) '_ACT_long']...
    ' ' [Site '_' sprintf('%02.0f',Deployment) '_ENT_long'] ' ' [Site '_' sprintf('%02.0f',Deployment) '_EVN_long']...
    ' ' [Site '_' sprintf('%02.0f',Deployment) '_tframe_long']]); 
 
    
    
