% Process continuous data from FKNMS
% S. Ricci 5/16/19

% Functions used:
%   readST
%   mktableSTdir
%   sound_MSPEC

%% load data, create y vector and reshape to make x matrix
load STcalibration
dir2process = readtable('dir2process_cont.csv');
%dir = '/Volumes/P3/FKNMS_soundscape/deployment8/LKSPA_Dep8_335851542_ST21/';
chunk_size = 30; %size of data chunk in seconds


for k = 1;
    Site = char(dir2process.Site(k));
    DirIn =char(dir2process.DirIn(k));             
    DirOut=char(dir2process.DirOut(k)); 
    Deployment=dir2process.Deployment(k);
    FS=dir2process.FS(k);
    %Sgate=datenum(dir2process.Sgate(k))-0.00416;  % 0.0041 days = 10 minutes 
    %Egate=datenum(dir2process.Egate(k))+0.00416;

% generate file list, file names and diretories 
if exist('DirOut','dir') ~= 1; eval(['system(''mkdir '  DirOut ''')']); end %  make output directory if it does not exist 

[filelist,ftimes,fend]= mktableSTdir(DirIn);

if FS==48000 
win= 2^14;   % number of points for mean spectrum averaging 
ovlp=0;  % overlap in spectral calcs 
end

if FS==96000 
win= 2^15;   % number of points for mean spectrum averaging 
ovlp=0;  % overlap in spectral calcs 
end

if FS==144000 
win= 2^16;   % number of points for mean spectrum averaging 
ovlp=0;  % overlap in spectral calcs 
end

%% readST for each file, sound_MSPEC for each 30 sec in each file
tic 
    for i = 9:15;%length(filelist)  % for each file in the list 

        fprintf('Processing %s\n', char(filelist(i).name)); %display file name of file being processed
        [y,fstart,fs,metadata]=readST_continuous(char(filelist(i).name),DirIn,STcalib);

    %reshape long y vector (~six hours) to matrix where each column is 30
    %seconds of data. Now you can use each column as a "file" for mspec
    %calculations.
        bwin = chunk_size * fs; %number of points for each column of matrix (30 seconds of data)
        x=buffer(y,bwin,ovlp,'nodelay');   
        if x(end)==0; x(:,end)=[]; end  % if last colum is zero padded delete. 
        [r,Nwin]=size((x)); 
        msub=repmat(mean(x,1),r,1); % calculate the mean of each column 
        x=x-msub; % remove the mean of each column

    % create poavg and rms variables
        poavg = nan(win/2+1,size(x,2));
        rms = nan(1,size(x,2));
        % loop through each column in x
        for j = 1:size(x,2);  % for each column in the x matrix 
            [poavg(:,j),f]=sound_MSPEC(x(:,j),fs,win,ovlp); 
            rms(j)= std(x(:,j)); % RMS in the time domain of this file 

        end
    %time axis: start time for each 30 sec data chunk based on fstart
    t = fstart + (0:1:size(x,2)-1) * (chunk_size/(24*60*60));

    %save variables for this file
    eval([Site '_' sprintf('%02.0f',Deployment),'_',datestr(fstart,'yyyymmddHHMMSS'),'_poavg=poavg;']); %poavg matrix (30 sec chunks)
    eval([Site '_' sprintf('%02.0f',Deployment),'_',datestr(fstart,'yyyymmddHHMMSS'),'_rms=rms;']); %rms for each 30 sec chunk
    eval([Site '_' sprintf('%02.0f',Deployment),'_',datestr(fstart,'yyyymmddHHMMSS'),'_t=t;']); %start time for each 30 sec chunk
    eval([Site '_' sprintf('%02.0f',Deployment),'_',datestr(fstart,'yyyymmddHHMMSS'), '_f=f;']); %freq vector for poavg matrix
    
    %save this file's vars in a mat file in DirOut
    eval(['save ',DirOut, [Site '_' sprintf('%02.0f',Deployment),'_',datestr(fstart,'yyyymmddHHMMSS'), '_results.mat']...
        ' ' [Site '_' sprintf('%02.0f',Deployment),'_',datestr(fstart,'yyyymmddHHMMSS'),'_poavg']...
        ' ' [Site '_' sprintf('%02.0f',Deployment),'_',datestr(fstart,'yyyymmddHHMMSS'),'_rms']...
        ' ' [Site '_' sprintf('%02.0f',Deployment),'_',datestr(fstart,'yyyymmddHHMMSS'),'_t']...
        ' ' [Site '_' sprintf('%02.0f',Deployment),'_',datestr(fstart,'yyyymmddHHMMSS'), '_f']]);
    %save filelist for whole deployment
    eval([Site '_' sprintf('%02.0f',Deployment) '_filelist=filelist;']);
    eval([Site '_' sprintf('%02.0f',Deployment) '_metadata=metadata;']);
    eval(['save ',DirOut, [Site '_' sprintf('%02.0f',Deployment) '_results.mat'] ' ' [Site '_' sprintf('%02.0f',Deployment) '_filelist']...
        ' ' [Site '_' sprintf('%02.0f',Deployment) '_metadata']]);

    clear poavg rms tchunk metadata x y fstart fs bwin
    end
    

end
toc    
%%

% files = 'lkspa_dep8_continuous_wavfiles.txt';
% fileID = fopen(files);
% lkspa_filelist = textscan(fileID,'%s');
% fclose(fileID);
% lkspa_filelist = lkspa_filelist{1};

% 
% 
% [~,F,T,Pxx]=spectrogram(y(1:fs*3600),2^12,2^11,[],fs); %0.05sec windows, 0% overlap
% t=(0:1:(length(y)-1))*(1/fs);
% t = t(1:(fs*3600));
% figure; imagesc(T,F,10*log10(Pxx));set(gca,'FontSize',12);
% caxis([40,90]);
% axis xy; xlabel('Time (sec)');
% ylim([0 25000]);
% set(gca,'Ytick',0:5000:25000,'YtickLabel',[0,5,10,15,20,25],'ygrid','on');
% ylabel('Frequency (kHz)'); xlim([0,max(t)]); set(gca,'xtick',0:5:max(t));
% title('25000 Hz Spectrogram (colourbar = 45-90 dB)'); colormap jet;
% 
% 
% 
% dirIn = '/Volumes/P3/FKNMS_soundscape/deployment8/LKSPA_Dep8_335851542_ST21/'; %needs ending / to work with mktableSTdir
% 
% if fs==48000 
% win= 2^14;   % number of points for mean spectrum averaging 
% ovlp=0;  % overlap in spectral calcs 
% poavg=nan(win/2+1,length(filelist));    % preallocate space  - Power in each freq bin and file 
% rms=nan(1,length(filelist));            % preallocate space - RMS pressure in each file  
% UTC=nan(1,length(filelist)); 
% end