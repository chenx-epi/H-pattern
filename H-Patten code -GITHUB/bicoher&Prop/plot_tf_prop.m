%% plot propagation (in prop.mat) 
%    from EZ to fit PZ (EZtoPZ)
[PZ1, t]= lsim(EZtoPZ.sys_tf,EZtoPZ.EZ,EZtoPZ.t);
figure('position',[600 100 600 200],'Name','EZ to PZ'); 
plot(t, PZ,'linewidth',2); hold on;
plot(t, PZ1,'r-','linewidth',2); xlim(XLIM); box off;
xlabel('Time (second)','interpreter','latex');
title('Estimating non-{\itd}H from {\itd}H');
set(gca,'fontsize',20);

%% plot propagation (in prop.mat) 
%    from PZ to fit EZ (PZtoEZ)
[EZ1, t]= lsim(PZtoEZ.sys_tf,PZtoEZ.PZ,EZtoPZ.t);
figure('position',[600 100 600 200],'Name','PZ to EZ'); 
plot(t, EZ, 'linewidth',2); hold on
plot(t, EZ1,'r-','linewidth',2); xlim(XLIM); box off;
title('Estimating {\itd}H from non-{\itd}H');
xlabel('Time (second)','interpreter','latex');
set(gca,'fontsize',20);
