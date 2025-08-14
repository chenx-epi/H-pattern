% The example of EZ (info_EZ in file bic_info.mat)
waxis = info_EZ.waxis;
bic   = info_EZ.bic;

fmax  = 200;  % max frequency for bic plot

%% plot amplitude
figure('visible','on','position',[600,100,600,250])
%contourf(waxis,waxis,abs(bic),[0.85,0.9],50,'LineStyle','none'), grid on
contourf(waxis,waxis,abs(bic),'LineStyle','none'), grid on
%Con.LineWidth =1;
hold on
plot([0,fmax/2],[0,fmax/2],'r-','linewidth',2);
plot([0,fmax],[0,0],'r-','linewidth',2);
plot([fmax/2,fmax],[fmax/2,0],'r-','linewidth',2);
%scatter(waxis(row),waxis(col),50,'*');
title('Bicoherence amplitude')
xlabel('f_1 (Hz)'), ylabel('f_2 (Hz)')
set(gcf,'Name','Hosa BICOHER amplitude');
set(gca,'clim',[0,1],'fontsize',20);
xlim([0,fmax]);
ylim([0,fmax/2]);
C  = parula(10); C(1:5,:)=ones(5,3);  % you can choose other colormap, like jet, hot etc.
Cs = colormap(C);
colorbar('Ticks',0:0.2:1);

%% plot skewness
figure('visible','on','position',[600,100,600,250])
contourf(waxis,waxis,real(bic),'LineStyle','none'), grid on
hold on
plot([0,fmax/2],[0,fmax/2],'r-','linewidth',2);
plot([0,fmax],[0,0],'r-','linewidth',2);
plot([fmax/2,fmax],[fmax/2,0],'r-','linewidth',2);
title('Skewness')
xlabel('f_1 (Hz)'), ylabel('f_2 (Hz)')
set(gcf,'Name','Hosa BICOHER skewness');
set(gca,'clim',[-1,1],'fontsize',20)
xlim([0,fmax]);
ylim([0,fmax/2]);
C  = parula(20); C(7:14,:)=ones(8,3);
Cs = colormap(C);
colorbar('Ticks',-1:0.2:1);

%% plot asymmetry
figure('visible','on','position',[600,100,600,250])
contourf(waxis,waxis,imag(bic),'LineStyle','none'), grid on
hold on
plot([0,fmax/2],[0,fmax/2],'r-','linewidth',2);
plot([0,fmax],[0,0],'r-','linewidth',2);
plot([fmax/2,fmax],[fmax/2,0],'r-','linewidth',2);
title('Asymmetry')
xlabel('f_1 (Hz)'), ylabel('f_2 (Hz)')
set(gcf,'Name','Hosa BICOHER asymmetry');
set(gca,'clim',[-1,1],'fontsize',20)
xlim([0,fmax]);
ylim([0,fmax/2]);
C  = parula(20); C(7:14,:)=ones(8,3);
Cs = colormap(C);
colorbar('Ticks',-1:0.2:1);
