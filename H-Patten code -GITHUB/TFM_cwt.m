%% using cwtft in spectrog with moarletParam = 100;
mp = 100;
CLIMs_EZ = [5, 40;  % 1
     5,25;          %2
    10,40;          %3
    10,40;          %4
     5,40;          %5
    25,45;          %6
    15,25;          %7
    10,30;          %8
    10,30;          %9
     0,40;          %10
     5,20;          %11
    10,70;          %12
    10,30;
     0,20;
    10,30;          %15
     0,30;
    10,35;
     0,30;
    10,35;
    10,40;          %20
    10,30;
     0,30;
     5,30;
     5,30;
     5,30];

CLIMs_PZ = [5, 40;
     5,25;
    10,40;
   -10,20;
     0,35;          %5
   -10,30;
     5,25;
    10,30;
    10,30;
     5,40;          %10
     5,20;
    10,70;          %12
    10,30;
     0,20;
    10,30;          %15
     0,30;
     0,20;
     0,15;
    10,20;
    10,30;          %20
    10,30;
     0,30;
     5,20;
     5,30;
     5,20];

 %% plot tfm
 % for time-frequency map
 
fmax = 200;
for i = 1:N
    fs      = FS(i);
    XLIM    = XLIMs(i,:);
    range   = fs*(XLIM(1)-lag):fs*(XLIM(2)+lag);
    t_range = range/fs;
    
    cfs = Spectrog(data_EZ{i}(1:100*fs),1/fs,1:fmax,mp);
    TFM = pow2db(abs(cfs).^2);
    base= mean(TFM(:,1:50*fs),2);
    cfsLVFA_EZ = Spectrog(LVFAm_EZ{i},1/fs,1:fmax,mp);
    TFM_LVFA_EZ= pow2db(abs(cfsLVFA_EZ).^2)-repmat(base,[1,size(cfsLVFA_EZ,2)]);
    
    cfs = Spectrog(data_PZ{i}(1:100*fs),1/fs,1:fmax,mp);
    TFM = pow2db(abs(cfs).^2);
    base= mean(TFM(:,1:50*fs),2);
    cfsLVFA_PZ = Spectrog(LVFA_PZ{i},1/fs,1:fmax,mp);
    TFM_LVFA_PZ= pow2db(abs(cfsLVFA_PZ).^2)-repmat(base,[1,size(cfsLVFA_PZ,2)]);
    
    figure()
    subplot(2,1,1)
    imagesc('XData',t_range,'YData',1:fmax,'CData',TFM_LVFA_EZ); hold on
    plot([XLIM(1),XLIM(1)],[0,fmax],'r--');
    plot([XLIM(2),XLIM(2)],[0,fmax],'r--');
    xlim(XLIM+[-lag,lag]);
    ylim([0,fmax]);
    ylabel('Frequency $f$ (Hz)','interpreter','latex');
    title(['Patient ',int2str(i),': TFM of ',types{i},' in EZ'],'interpreter','latex');
    caxis(CLIMs_EZ(i,:));
    colormap jet;
    c=colorbar;
    ct=get(c,'Title');
    set(ct,'String','$\Delta$db','interpreter','latex')
    shading flat
    set(gca,'YDir','normal','fontsize',20);
    subplot(2,1,2)
    imagesc('XData',t_range,'YData',1:fmax,'CData',TFM_LVFA_PZ); hold on
    plot([XLIM(1),XLIM(1)],[0,fmax],'r--');
    plot([XLIM(2),XLIM(2)],[0,fmax],'r--');
    xlim(XLIM+[-lag,lag]);
    ylim([0,fmax]);
    ylabel('Frequency $f$ (Hz)','interpreter','latex');
    xlabel('Time (second)','interpreter','latex');
    title(['Patient ',int2str(i),': TFM of ',types{i},' in PZ'],'interpreter','latex');
    caxis(CLIMs_PZ(i,:));
    colormap jet;
    c=colorbar;
    ct=get(c,'Title');
    set(ct,'String','$\Delta$db','interpreter','latex')
    shading flat
    set(gca,'YDir','normal','fontsize',20);
    savefig(['./figs/tfm_Pat_',int2str(i),'.fig']);
    
    close all
end
