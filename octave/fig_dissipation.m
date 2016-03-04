%-------------------------------------%
% Figure of the Laplacian Dissipation %  
%-------------------------------------%
nx = 64; 


%--- small-scale dissipation
ksig   =  20;
psig   =  8;
rtsig  =  0.01;

%--- large-scale dissipation
klam   =  2;
plam   = -2;
rtlam  =  0.005;

%--- wave number grid
kmax   = fix(nx/3);

%============== internal part ============

%--- spectral grid
k =  1:kmax;

sig = rtsig*(1/ksig)^(2*psig);
lam = rtlam*(1/klam)^(2*plam);

r   = k;
r2  = k.^2;

D_sig   = -sig*r2.^psig;
D_lam   = -lam*r2.^plam;

D       = D_sig + D_lam;

ff = figure
ax = axes
set(ax,'XLim',[1,kmax]);
hold
pp = plot(r,D_sig + D_lam,'r');
set(pp,'LineWidth',[2]);


