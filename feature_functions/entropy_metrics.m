function [EPS,EAS,ECV] = entropy_metrics(po_amp2_nr,lf_lim,f)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

lf = f >=lf_lim(1) & f<= lf_lim(2);

for i = 1:1:size(po_amp2_nr,3);
    tmp = po_amp2_nr(lf,:,i);
    
    %EPS (Entropy of Spectral Peaks)
    [maxval,maxind] = max(tmp,[],1);
    [maxcount,maxbins]=histcounts(maxind,size(tmp,1));
    
    sumXi=sum(maxcount); % Sum the envelope 
    At=maxcount/sumXi;  % Compute the probability mass function 
    Hp=-(nansum(At .* log2(At)))/log2(length(At));
    EPS(i) = 1 - Hp;
    
    
    %EAS (Entropy of Average Spectrum)
    poavg = nanmean(tmp,2);
    sumXi_poavg=sum(poavg); % Sum the envelope 
    At_poavg=poavg/sumXi_poavg;  % Compute the probability mass function 
    Ha=-(nansum(At_poavg .* log2(At_poavg)))/log2(length(At_poavg));
    EAS(i) = 1 - Ha;

    %ECV (Entropy of the Spectrum of Coeffecients of Variation)
    povar = var(tmp,[],2);
    pocv = povar./poavg;
    
    %Xi=abs(hilbert(poavg));
    Xi=pocv;
    sumXi_pocv=sum(pocv); % Sum the envelope 
    At_pocv=pocv/sumXi_pocv;  % Compute the probability mass function 
    Hc=-(nansum(At_pocv .* log2(At_pocv)))/log2(length(At_pocv));
    ECV(i) = 1 - Hc;

end
end

