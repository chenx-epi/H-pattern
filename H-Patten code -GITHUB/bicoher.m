function info = bicoher(y,  nfft, wind, nsamp, overlap, Fs, fmax,folder,pre,i)
%BICOHER - Direct (FD) method for estimating bicoherence
%	[bic,waxis] = bicoher (y,  nfft, wind, segsamp, overlap)
%	y     - data vector or time-series
%	nfft - fft length [default = power of two > segsamp]
%	       actual size used is power of two greater than 'nsamp'
%	wind - specifies the time-domain window to be applied to each
%	       data segment; should be of length 'segsamp' (see below);
%		otherwise, the default Hanning window is used.
%	segsamp - samples per segment [default: such that we have 8 segments]
%	        - if x is a matrix, segsamp is set to the number of rows
%	overlap - percentage overlap, allowed range [0,99]. [default = 50];
%	        - if x is a matrix, overlap is set to 0.
%	bic     - estimated bicoherence: an nfft x nfft array, with origin
%	          at the center, and axes pointing down and to the right.
%	waxis   - vector of frequencies associated with the rows and columns
%	          of bic;  sampling frequency is assumed to be 1.
%   Fs      - sampling frequency
%   fmax    - maximal frequency in Bspec

%   info.phi- nonlinearity of the signal

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc.
%       $Revision: 1.7 $
%  A. Swami   January 20, 1995

%     RESTRICTED RIGHTS LEGEND
% Use, duplication, or disclosure by the Government is subject to
% restrictions as set forth in subparagraph (c) (1) (ii) of the
% Rights in Technical Data and Computer Software clause of DFARS
% 252.227-7013.
% Manufacturer: United Signals & Systems, Inc., P.O. Box 2374,
% Culver City, California 90231.
%
%  This material may be reproduced by or for the U.S. Government pursuant
%  to the copyright license under the clause at DFARS 252.227-7013.


% --------------------- parameter checks -----------------------------

[ly, nrecs] = size(y);
if (ly == 1) y = y(:);  ly = nrecs; nrecs = 1; end

if (exist('nfft') ~= 1)            nfft = 128; end
if (exist('overlap') ~= 1)      overlap = 50;  end
overlap = max(0,min(overlap,99));
if (nrecs > 1)                  overlap = 0;   end
if (exist('nsamp') ~= 1)          nsamp = 0;  end
if (nrecs > 1)                    nsamp = ly;  end

if (nrecs == 1 & nsamp <= 0)
    nsamp = fix(ly/ (8 - 7 * overlap/100));
end
if (nfft  < nsamp)   nfft = 2^nextpow2(nsamp); end

overlap  = fix( nsamp * overlap/100);
nadvance = nsamp - overlap;
nrecs    = fix ( (ly*nrecs - overlap) / nadvance);

% ----------------------------------------------------------------------
if (exist('wind') ~= 1) wind = hanning(nsamp); end
[rw,cw] = size(wind);
if (min(rw,cw) ~= 1 | max(rw,cw) ~= nsamp)
    disp(['Segment size  is ',int2str(nsamp)])
    disp(['"wind" array  is ',int2str(rw),' by ',int2str(cw)])
    disp(['Using default Hanning window'])
    wind = hanning(nsamp);
end
wind = wind(:);

% ---------------- accumulate triple products ----------------------
bic  = zeros(nfft,nfft);
%Pyy  = zeros(nfft,1);
Py1y2= zeros(nfft,nfft);
Py12 = zeros(nfft,nfft);

mask = hankel([1:nfft],[nfft,1:nfft-1]);   % the hankel mask (faster)
ind  = [1:nsamp];
yvar = 0;

for k = 1:nrecs
    ys = y(ind);
    ys = (ys(:) - mean(ys)) .* wind;
    Yf = fft(ys,nfft)  / nsamp;
    CYf     = conj(Yf);
    %Pyy     = Pyy + Yf .* CYf;
    yvar=yvar+var(ys);
    
    Yf1f2= Yf * Yf.';
    Yf12 = CYf(mask);
    bic  = bic + Yf1f2 .* Yf12;
    Py1y2   = Py1y2 + abs(Yf1f2).^2;
    Py12    = Py12 + abs(Yf12).^2;
    
    ind = ind + nadvance;
end

yvar    = yvar / nrecs;
bic1     = bic / nrecs;
%Pyy     = Pyy / nrecs;
%Pmask   = Pyy(mask);
%bic = abs(bic).^2 ./ (Pyy * Pyy.' .* Pmask);

Py1y2= Py1y2 / nrecs;
Py12 = Py12 / nrecs;
bic  = bic1 ./ sqrt(Py1y2 .* Py12);
bic0 = fftshift(bic);

% ------------ contout plot of magnitude bispectum --------------------

fs=Fs/nsamp;
m=round(fmax/fs)+1;
%waxis = (-m:m)'*fs;
if (rem(nfft,2) == 0)
    waxis0 = (-nfft/2:(nfft/2-1))'.*fs;
else
    waxis0 = (-(nfft-1)/2:(nfft-1)/2)'.*fs;
end
range = nfft/2+(1:m);
waxis = waxis0(range);
bic   = bic0(range,range);

% ------------- total nonlinear interaction --------------
bicsum= sum(sum(tril(bic1(1:round(m/2),:))))+ sum(sum(triu(fliplr(bic1(round(m/2)+(1:round(m/2)),:)))));
bicsum= bicsum./sqrt(yvar).^3;

[colmax,row] = max(abs(bic));
[maxval,col] = max(colmax);
row = row(col);
disp([' Sum: bic = ',num2str(bicsum),newline,' Max: bic(',num2str(waxis(row)),',',num2str(waxis(col)),') = ', ...
    num2str(maxval)])

info.bicsum = bicsum;
info.frow   = waxis(row);
info.fcol   = waxis(col);
info.maxval = maxval;
info.bic    = bic;
info.waxis  = waxis;

% set domain
mask = triu(ones(size(bic)));
nrow = size(mask,1);
mask(floor(nrow/2)+1:end,:)=0;
mask = fliplr(triu(fliplr(mask)));   

thres = sqrt(6/nrecs);
mask_t = abs(bic)>thres;

bicabs = abs(bic).*mask.*mask_t;
Is  = ones(size(bic)); Is = Is-eye(nrow)./2;
phi = 8.*bicabs^2.*fs^2.*Is;
info.phi = sqrt(mean(mean(phi)));
disp([' Nonlinearity: phi = ', num2str(info.phi)]);

mask(mask==0)=nan;
% -------------- plot bicoherence ---------------------
figure('visible','on','position',[600,100,600,250])
contourf(waxis,waxis,abs(bic).*mask.*mask_t,50,'LineStyle','none'), grid on
hold on
plot([0,fmax/2],[0,fmax/2],'r-','linewidth',2);
plot([0,fmax],[0,0],'r-','linewidth',2);
plot([fmax/2,fmax],[fmax/2,0],'r-','linewidth',2);
%scatter(waxis(row),waxis(col),50,'*');
title('Bicoherence amplitude')
xlabel('f_1 (Hz)'), ylabel('f_2 (Hz)')
set(gcf,'Name','Hosa BICOHER amplitude');
set(gca,'clim',[0,1],'fontsize',20)
xlim([0,fmax]);
ylim([0,fmax/2]);
C = parula(50); cnum = floor(thres.*50); C(1:cnum,:)=ones(cnum,3);
colormap(C)
colorbar
savefig(gcf,[folder,pre,'_bic_amp_Pat_',num2str(i),'.fig'],'compact');
saveas(gcf,[folder,pre,'_bic_amp_Pat_',num2str(i),'.pdf']);
close all;

figure('visible','on','position',[600,100,600,250])
contourf(waxis,waxis,real(bic).*mask,50,'LineStyle','none'), grid on
hold on
plot([0,fmax/2],[0,fmax/2],'r-','linewidth',2);
plot([0,fmax],[0,0],'r-','linewidth',2);
plot([fmax/2,fmax],[fmax/2,0],'r-','linewidth',2);
title('Skewness')
xlabel('f_1 (Hz)'), ylabel('f_2 (Hz)')
set(gcf,'Name','Hosa BICOHER real');
set(gca,'clim',[-1,1]/2,'fontsize',20)
xlim([0,fmax]);
ylim([0,fmax/2]);
%colormap jet
colorbar
savefig(gcf,[folder,pre,'_bic_skew_Pat_',num2str(i),'.fig'],'compact');
saveas(gcf,[folder,pre,'_bic_skew_Pat_',num2str(i),'.pdf']);
close all;

figure('visible','on','position',[600,100,600,250])
contourf(waxis,waxis,imag(bic).*mask,50,'LineStyle','none'), grid on
hold on
plot([0,fmax/2],[0,fmax/2],'r-','linewidth',2);
plot([0,fmax],[0,0],'r-','linewidth',2);
plot([fmax/2,fmax],[fmax/2,0],'r-','linewidth',2);
title('Asymmetry')
xlabel('f_1 (Hz)'), ylabel('f_2 (Hz)')
set(gcf,'Name','Hosa BICOHER imag');
set(gca,'clim',[-1,1]/2,'fontsize',20)
xlim([0,fmax]);
ylim([0,fmax/2]);
%colormap jet
colorbar
savefig(gcf,[folder,pre,'_bic_asym_Pat_',num2str(i),'.fig'],'compact');
saveas(gcf,[folder,pre,'_bic_asym_Pat_',num2str(i),'.pdf']);
close all;

return
