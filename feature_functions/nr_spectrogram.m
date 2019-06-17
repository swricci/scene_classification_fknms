function [po_db_nr,mode_db_freq] = nr_spectrogram(po_db,mavg_win)
%nr_spectrogram Removes background noise from spectrogram.
% For each minute segment, a background noise profile is constructed for
% each frequency bin. This is subtracted from the decibel spectrogram,
% resulting in a noise-removed spectrogram

%pre-allocate space
mode_db_freq = nan(size(po_db,3),size(po_db,1));

%% Create histogram of intensity values in each freq bin

for c = 1:size(po_db,3);
    tmp = po_db(:,:,c);
    hcount = nan(size(tmp,1),100);
    edges = nan(size(tmp,1),101);
    
    for i = 1:size(tmp,1);
        [hcount(i,:),edges(i,:)] = histcounts(tmp(i,:),100); %output bins here too
    end
    
%% Smooth with moving avg filter
    smooth_hcount_po = nan(size(hcount,1),100);
    for i=1:size(hcount,1);
        smooth_hcount_po(i,:) = smoothdata(hcount(i,:),'movmean',mavg_win);
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
        mode_db_freqbin(c,i) = mean(edges(i,h_ind(i):h_ind(i) +1));
        end
    end
    
%% Subtract noise profile from db spectrogram to remove noise
    %smooth noise profile
    mode_db_freq(c,:) = smoothdata(mode_db_freqbin(c,:),'movmean',mavg_win);
    %figure; plot(f,smooth_mode_db_freq);
    
    %subtract from po_db matrix to remove noise
    po_db_nr(:,:,c) = tmp - mode_db_freq(c,:)';
    
    %negative values to 0
    a = find(po_db_nr < 0);
    po_db_nr(a) = 0;
end