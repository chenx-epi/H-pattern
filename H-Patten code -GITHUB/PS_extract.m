ind = find(PS>1.9e-4);  % find the index where PS>1.9e-4
di  = diff(ind);        % find segment jump
% find the stard and end index of each segment
pstart = ind(find(di>1)+1); pstart = [ind(1),pstart];
pend   = ind(di>1);         pend   = [pend,ind(end)];

segs_ind = [pstart',pend'];
PS_nopulse = PS;

% for each segment, extract the pulse
for i = 1:size(segs_ind,1)
    ind_i = segs_ind(i,1):segs_ind(i,2);
    [max_i,mi] = max(PS(ind_i));
    max_i_ind  = ind_i(mi);  % find the index for the pulse peak
    
    ind_i1 = max(1,segs_ind(i,1)-20):max_i_ind;
    seg_i1 = PS(ind_i1);
    sd_i1  = sign(diff(seg_i1));
    sd_i1_minus = find(sd_i1<0);
    ind_i1_left = ind_i1(sd_i1_minus(end))+1;   % index of left side of pulse
    
    ind_i2 = max_i_ind:segs_ind(i,2)+20;
    if ind_i2(end)>numel(PS)
        ind_i2 = max_i_ind:numel(PS);
    end
    seg_i2 = PS(ind_i2);
    sd_i2  = sign(diff(seg_i2));
    sd_i2_plus  = find(sd_i2>0);
    if isempty(sd_i2_plus)
        ind_i2_right = numel(PS);
    else
        ind_i2_right= ind_i2(sd_i2_plus(1));      % index of right side of pulse
    end
    
    ind_i_large = max(1,segs_ind(i,1)-20):segs_ind(i,2)+20;
    if ind_i_large(end)>numel(PS)
        ind_i_large = max_i_ind:numel(PS);
    end
    %plot(t_range(ind_i_large),PS(ind_i_large),'k-o'); hold on
    %plot(ind_i_large,PS(ind_i_large),'k-o'); hold on
    ind_i_pulse = ind_i1_left:ind_i2_right;
    %plot(t_range(ind_i_pulse),PS(ind_i_pulse),'r-*');
    %plot(ind_i_pulse,PS(ind_i_pulse),'r-*');
    
    x = [segs_ind(i,1)-20:ind_i1_left,ind_i2_right:segs_ind(i,2)+20];
    if segs_ind(i,2)+10>numel(PS)
        x=[max(1,segs_ind(i,1)-10):ind_i1_left,ind_i2_right:numel(PS)];
    end
    v = PS(x);
    xq= ind_i1_left+1:ind_i2_right-1;
    %vq=interp1(x,v,xq,'linear','extrap');
    vq=interp1(x,v,xq,'spline','extrap');
    if max(vq)>PS(max_i_ind)
        vq=interp1(x,v,xq,'linear','extrap');
    end
    if min(vq)<PS(ind_i1_left)
        vq=interp1(x,v,xq,'linear','extrap');
        disp(i);
    end
    
    PS_nopulse(xq) = vq;
%     
%     figure()
%     plot(t_range(ind_i_large),PS(ind_i_large)); hold on
%     plot(t_range(ind_i_large),PS_nopulse(ind_i_large));
end

figure()
plot(t_range,PS,'linewidth',2); hold on
plot(t_range,PS_nopulse); hold on
modu= smooth(t_range,PS_nopulse,0.08,'loess'); modu = modu';
plot(t_range,modu,'r-');

%% using cwtft in spectrog with moarletParam = 100;
mp =100;

cfs = Spectrog(data133,1/fs,1:300,mp);
TFM = pow2db(abs(cfs).^2);
base= mean(TFM(:,1:50*fs),2);
TFM_base=TFM; %-repmat(base,[1,size(TFM,2)]);

cfsPS = Spectrog(PS,1/fs,1:300,mp);
TFM_PS= pow2db(abs(cfsPS).^2);%-repmat(base,[1,size(cfsPS,2)]);

cfsPSn = Spectrog(PS_nopulse,1/fs,1:300,mp);
TFM_PSn= pow2db(abs(cfsPSn).^2);%-repmat(base,[1,size(cfsPS,2)]);


figure()
imagesc('XData',t_range,'YData',1:300,'CData',TFM_PS);
xlim([t_range(1),t_range(end)]);
ylim([0,200]);
ylabel('Frequency $f$ (Hz)','interpreter','latex');
xlabel('Time (second)','interpreter','latex');
title('TFM of PS','interpreter','latex');
caxis([-90,-50]);
c=colorbar;
ct=get(c,'Title');
set(ct,'String','db','interpreter','latex')
shading flat
set(gca,'YDir','normal','fontsize',20);

figure()
imagesc('XData',t_range,'YData',1:300,'CData',TFM_PSn);
xlim([t_range(1),t_range(end)]);
ylim([0,200]);
ylabel('Frequency $f$ (Hz)','interpreter','latex');
xlabel('Time (second)','interpreter','latex');
title('TFM of PS with no pulse','interpreter','latex');
caxis([-90,-50]);
c=colorbar;
ct=get(c,'Title');
set(ct,'String','db','interpreter','latex')
shading flat
set(gca,'YDir','normal','fontsize',20);
