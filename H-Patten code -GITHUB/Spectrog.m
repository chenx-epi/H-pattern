function cfs=Spectrog(phi_e,dt,freqs,mp)

% Set up paramters for morlet wavelet transform
if isempty(mp) || nargin<4, mp=7; end   %morletParam = 7;
morletFourierFactor = 4*pi/(mp+sqrt(2+mp^2));

% Set up structure defining scales between min and max pseudo-frequencies
scales = 1./(morletFourierFactor * freqs);

cfstruct = cwtft({phi_e,dt},'scales',scales,'wavelet',{'morl',mp});
cfs = cfstruct.cfs;
