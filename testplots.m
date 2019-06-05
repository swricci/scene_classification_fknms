ts = cat(2, LKP_08_20190130055936_rms,LKP_08_20190130115933_rms,LKP_08_20190130175930_rms,LKP_08_20190130235928_rms,...
    LKP_08_20190131055925_rms,LKP_08_20190131115923_rms,LKP_08_20190131175920_rms,LKP_08_20190131235917_rms,...
    LKP_08_20190201055914_rms,LKP_08_20190201115912_rms,LKP_08_20190201175910_rms,LKP_08_20190201235907_rms);

time = cat(2, LKP_08_20190130055936_t,LKP_08_20190130115933_t,LKP_08_20190130175930_t,LKP_08_20190130235928_t,...
    LKP_08_20190131055925_t,LKP_08_20190131115923_t,LKP_08_20190131175920_t,LKP_08_20190131235917_t,...
    LKP_08_20190201055914_t,LKP_08_20190201115912_t,LKP_08_20190201175910_t,LKP_08_20190201235907_t);

figure; plot(time,20*log10(ts));

power = cat(2, LKP_08_20190130055936_poavg,LKP_08_20190130115933_poavg,LKP_08_20190130175930_poavg,LKP_08_20190130235928_poavg,...
    LKP_08_20190131055925_poavg,LKP_08_20190131115923_poavg,LKP_08_20190131175920_poavg,LKP_08_20190131235917_poavg,...
    LKP_08_20190201055914_poavg,LKP_08_20190201115912_poavg,LKP_08_20190201175910_poavg,LKP_08_20190201235907_poavg);


freq = cat(2, LKP_08_20190130055936_f,LKP_08_20190130115933_f,LKP_08_20190130175930_f,LKP_08_20190130235928_f,...
    LKP_08_20190131055925_f,LKP_08_20190131115923_f,LKP_08_20190131175920_f,LKP_08_20190131235917_f,...
    LKP_08_20190201055914_f,LKP_08_20190201115912_f,LKP_08_20190201175910_f,LKP_08_20190201235907_f);


figure; imagesc(time,freq,10*log10(power));set(gca,'FontSize',12); 
caxis([40,90]);  
axis xy; colormap jet;
datetick('x');

figure; imagesc(time,freq,10*log10(power));set(gca,'FontSize',12); 
caxis([40,90]);  ylim([0 3000]);
axis xy; colormap jet;
datetick('x');


ts = cat(2, LKP_08_20190130055936_rms, LKP_08_20190131055925_rms);

time = cat(2, LKP_08_20190130055936_t,LKP_08_20190131055925_t);

figure; plot(time,20*log10(ts));

power = cat(2, LKP_08_20190130055936_poavg,LKP_08_20190131055925_poavg);


freq = cat(2, LKP_08_20190130055936_f, LKP_08_20190131055925_f);

figure; imagesc(time,freq,10*log10(power));set(gca,'FontSize',12); 
caxis([40,90]);  
axis xy; colormap jet;
datetick('x');


figure; imagesc(LKP_08_20190130055936_t,LKP_08_20190130055936_f,10*log10(LKP_08_20190130055936_poavg));set(gca,'FontSize',12); 
axis xy; colormap jet;
datetick('x');

figure; imagesc(LKP_08_20190131055925_t,LKP_08_20190131055925_f,10*log10(LKP_08_20190131055925_poavg));set(gca,'FontSize',12); 
axis xy; colormap jet;
datetick('x');