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
LFC_long = [];
HFC_long = [];
EPS_long = [];
EAS_long = [];
ECV_long = [];
ACI_long = [];
pkfreq_long = [];
tframe_long = [];

DirIn =char(dir2process.DirIn(1));
Site=char(dir2process.Site(1));   
DirOut=char(dir2process.DirOut(1)); 
Deployment=dir2process.Deployment(1);
   
[filelist,ftimes,fend]= mktableSTdir(DirIn);

h = waitbar(0,'Calculating feature vectors...');

tic
for file = 2:length(filelist);
    fprintf('Processing %s\n', char(filelist(file).name));
    [y,fstart,fs,metadata]=readST_continuous(char(filelist(file).name),DirIn,STcalib);
    
    %filter data
    %high pass filter 50 Hz
    [B,A]=butter(3,[50, 3000]/(fs/2),'bandpass'); %may need to adjust this?
    [H,F] = freqz(B,A,2048,fs);
    %figure; plot(F,abs(H));
    yfilt = filtfilt(B,A,y);
    %yfilt = bandpass(y,[50 3000],fs); %took a VERY long time to run
    
    z = chunk_data(yfilt,fs,chunk_size,frame_size,ovlp);
  
    %calculate wave envelope
    wave_env = wave_env_func(z);
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
    
    %spectrogram based calculatons
    
    %create amp spectrogram
    [po_amp,f] = feature_spectrogram(yfilt,fs,chunk_size,frame_size,ovlp);
    
    
    %noise removal from spectrogram
    %dB
   po_db = 10*log10(po_amp);
   [po_db_nr,mode_db_freq]=nr_spectrogram(po_db,mavg_win);
    
   
   for k = 1:1:size(po_amp,3);
       tmp = nanmean(po_amp(:,:,k),2);
       [mval,ind] = max(tmp);
       
       pkfreq(k) = f(ind);
   end
   
   peakfreq(file,:) = pkfreq;
   
    %amp
    po_amp2 = po_amp.^2;
    [po_amp2_nr] = nr_spectrogram(po_amp2,mavg_win);
    
    %spectrogram based indices
    bp_lim = [50 3000];
    lf_lim = [50 1500];
    hf_lim = [1500 3000];
    
    %noise removed db spectrogram indices:
    [LFC(file,:),HFC(file,:)]= freq_cover_ind(po_db_nr,lf_lim,hf_lim,threshold,f);
    
    %noise removed amp2 spectrogram indices:
    [EPS(file,:),EAS(file,:),ECV(file,:)]=entropy_metrics(po_amp2_nr,bp_lim,f);
    
    %ACI
    for j=1:1:size(po_amp,3);
        filt = f >=bp_lim(1) & f<= bp_lim(2);
        tmp = po_amp(filt,:,j);
        %tmp = po_amp(:,:,j);

        I = abs(tmp);
        % difference in intensity in each freq bin from one time step to next
        diffI=abs(diff(I,1,2));  
        % sum of intensities in each freq bin
        sumI=sum(I,2);  
        % sum of differences between freq bins  
        sumdiffI=sum(diffI,2);  
        % ACI in each frequency bin 
        ACIfreq(:,j)=sumdiffI./sumI;  

        % sum those frequency bins to calcualte the ACI for that time window 
        ACItmp(j)=sum(ACIfreq(:,j)); 
    end
    
    ACI(file,:) = ACItmp;
    
    %determine start time of each frame
    tframe(file,:) = fstart + (0:1:size(wave_env,1)-1) * (chunk_size/(24*60*60));
    
    % create deployment long index
    BGN_long = cat(2,BGN_long,BGN(file,:));
    SNR_long = cat(2,SNR_long,SNR(file,:));
    ACT_long = cat(2,ACT_long,ACT(file,:));
    EVN_long = cat(2,EVN_long,EVN(file,:));
    ENT_long = cat(2,ENT_long,ENT(file,:));
    LFC_long = cat(2,LFC_long,LFC(file,:));
    HFC_long = cat(2,HFC_long,HFC(file,:));
    EPS_long = cat(2,EPS_long,EPS(file,:));
    EAS_long = cat(2,EAS_long,EAS(file,:));
    ECV_long = cat(2,ECV_long,ECV(file,:));
    ACI_long = cat(2,ACI_long,ACI(file,:));
    pkfreq_long = cat(2,pkfreq_long,peakfreq(file,:));
    tframe_long = cat(2,tframe_long, tframe(file,:));
    
    waitbar(file-1/length(file),h);
    
end
    
  toc  
  
% eval([Site '_' sprintf('%02.0f',Deployment) '_BGN=BGN;']);
% eval([Site '_' sprintf('%02.0f',Deployment) '_SNR=SNR;']);
% eval([Site '_' sprintf('%02.0f',Deployment) '_ACT=ACT;']);
% eval([Site '_' sprintf('%02.0f',Deployment) '_EVN=EVN;']);
% eval([Site '_' sprintf('%02.0f',Deployment) '_ENT=ENT;']);
% eval([Site '_' sprintf('%02.0f',Deployment) '_tframe=tframe;']);

  
eval([Site '_' sprintf('%02.0f',Deployment) '_BGN_long=BGN_long;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_SNR_long=SNR_long;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_ACT_long=ACT_long;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_ENT_long=ENT_long;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_EVN_long=EVN_long;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_LFC_long=LFC_long;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_HFC_long=HFC_long;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_EPS_long=EPS_long;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_EAS_long=EAS_long;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_ECV_long=ECV_long;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_ACI_long=ACI_long;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_pkfreq_long=pkfreq_long;']);
eval([Site '_' sprintf('%02.0f',Deployment) '_tframe_long=tframe_long;']);

eval(['save ',DirOut, [Site '_' sprintf('%02.0f',Deployment) '_featvec_results.mat']...
    ' ' [Site '_' sprintf('%02.0f',Deployment) '_BGN_long'] ' ' [Site '_' sprintf('%02.0f',Deployment) '_ACI_long']...
    ' ' [Site '_' sprintf('%02.0f',Deployment) '_SNR_long'] ' ' [Site '_' sprintf('%02.0f',Deployment) '_ACT_long']...
    ' ' [Site '_' sprintf('%02.0f',Deployment) '_ENT_long'] ' ' [Site '_' sprintf('%02.0f',Deployment) '_EVN_long']...
    ' ' [Site '_' sprintf('%02.0f',Deployment) '_LFC_long'] ' ' [Site '_' sprintf('%02.0f',Deployment) '_HFC_long']...
    ' ' [Site '_' sprintf('%02.0f',Deployment) '_EPS_long'] ' ' [Site '_' sprintf('%02.0f',Deployment) '_EAS_long']...
    ' ' [Site '_' sprintf('%02.0f',Deployment) '_ECV_long'] ' ' [Site '_' sprintf('%02.0f',Deployment) '_tframe_long']]); 