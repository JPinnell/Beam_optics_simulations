%% Mode converter HG_mn <-> LG_pl
% v1. Jonathan Pinnell (2020)
% This script simulates an optical mode converter that can convert between
% LG_pl and HG_mn modes such that p = min(m,n), l = m-n.
% Based on  Beijersbergen et al. "Astigmatic laser mode converters and 
% transfer of orbital angular momentum" (1992).
% To alternate between an initial LG and HG mode, uncomment lines 15-16 
% or 17-18.

% Make coordinates
H = 1000; dx = 8e-3; x = dx.*(-H/2:(H/2-1)); [X,Y] = meshgrid(x,-x);
[Phi,R] = cart2pol(X,Y);

lambda = 633e-6; % wavelength of light

w0 = 0.6; % waist radius
M = 5; N = 3; % HG mode indices
Mode0 = HG((X+Y)/sqrt(2),(X-Y)/sqrt(2),M,N,1,w0); % diagonally orientated HG mode
% P = 3; L = 0; % LG mode indices
% Mode0 = LG(R,Phi,P,L,1,w0); % LG mode

f = 200; % focal length of normal lens before mode converter
T = exp(-1i.*pi./(lambda.*f).*(X.^2+Y.^2)); % tranmission function of normal lens

% Calculate focal length of cylindrical lens
wf = lambda*f/pi/w0; % beam size at focal plane
zR = pi*wf^2/lambda; % rayleigh range of focused beam
fc = zR/(1 + 1/sqrt(2)); % focal length of cylindrical lens
Tc = exp(-1i.*pi./(lambda.*fc).*(X.^2)); % tranmission function of cylindrical lens

d = sqrt(2)*fc; % distance between cylindrical lenses

% Mode converter (see original paper for experimental setup arrangement)
Mode = FresnelProp(Mode0,dx,lambda,f);
Mode = Mode.*T;
Mode = FresnelProp(Mode,dx,lambda,f-d/2);
Mode = Mode.*Tc;
Mode = FresnelProp(Mode,dx,lambda,d);
Mode = Mode.*Tc;
Mode = FresnelProp(Mode,dx,lambda,f-d/2);
Mode = Mode.*T;
Mode = FresnelProp(Mode,dx,lambda,f);

% plots
Q = 400;
figure('color','w','units','pixels','position',[100 100 2*Q Q]);

subplot(1,2,1);
imagesc(abs(Mode0));
text(H/20,H/20,'Before','FontSize',20);
set(gca,'units','pixels','position',[0 0 Q Q],'visible','off');

subplot(1,2,2);
imagesc(abs(Mode));
text(H/20,H/20,'After','FontSize',20);
set(gca,'units','pixels','position',[Q 0 Q Q],'visible','off');