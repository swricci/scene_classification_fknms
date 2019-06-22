function [LFC,HFC] = freq_cover_ind(po_db_nr,lf_lim,hf_lim,threshold,f)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
for d =1:1:size(po_db_nr,3);
    tmp_nr = po_db_nr(:,:,d);
    
    lf = f >=lf_lim(1) & f<= lf_lim(2);
    hf = f >=hf_lim(1) & f<= hf_lim(2);
    
    tmp_nr_lf =tmp_nr(lf,:);
    tmp_nr_hf = tmp_nr(hf,:);
    
    ncell_lf = length(find(tmp_nr_lf > threshold));
    ncell_hf = length(find(tmp_nr_hf > threshold));
    
    LFC(d) = ncell_lf/numel(tmp_nr_lf);
    HFC(d) = ncell_hf/numel(tmp_nr_hf);
end
end

