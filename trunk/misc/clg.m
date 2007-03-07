% constructionbyaddition lowg version
% DFT the old fashioned way
% msb, dec 2006
%% constants
clear
coex=-2i*pi;
sq2=sqrt(2);
twopi=2*pi; 
gam=4258.0; % gamma, physical constant for proton imaging
dcoffset=20+1i*32;
% datafile='pinecone_19';
%datafile='four_tube_phantom_lowg/data/pinecone_01';
datafile='kiwi_lowg_pinecone_s_20070105_02/data/pinecone_02';
ktrajfile='lowg_circ_ktraj.mat';
%% initialize
fov=12.8; % field of view in cm
nres=64; % nominal resolution (x,y) in pixels, 64 is good
[XX,YY]=meshgrid(((1:nres)-floor(nres/2)-1)*fov/nres); % meshgrid for evaluation
load(ktrajfile); % load the measured k trajectory
% gmax=4.6833; % high g version maximum gradient used, not the varian system max
gmax=3.643;  % low gradient version
Ta=.065;                % duration of acquisition in ms
kf=nres/(2*fov);
delt=1.0/198511.166253; % defined from sw in procpar
NNN=floor(Ta/delt);trimoff=230;  
kss=ktraj(trimoff+1:NNN+trimoff)+.09+.013-((1:NNN)'*2.22e-06)-.05i;
phi=real(phi(trimoff+1:NNN+trimoff))+imag(phi(trimoff+1:NNN+trimoff));

kss=kss*2.0128;
kx=real(kss); 
ky=imag(kss);
delt=1./(gam*gmax*fov); % dwell time
[the_data_RE,the_data_IM]=LOAD_FID(datafile,2);
the_data=the_data_RE+1i*the_data_IM + dcoffset;
freqoffset=240;
fofactor=exp(-2i*pi*freqoffset*(1:NNN)*delt);
thesignal=the_data(trimoff+1:NNN+trimoff).'.*exp(-1i*phi).'.*fofactor;
%thesignal=ones(1,NNN);
%% density compensation
kr=abs(kss); % k radial
kfe=max(kr); % actual radius of k-space disc
%kfe=kf*sq2; % efective kf
% denscomp=(kr+eps).*sqrt(1-(kr./kfe).^2); % make the filter for a rossette
%  denscomp=kr.*sqrt(1-(kr./kfe).^2); % make the filter for a rossette
%  denscomp=denscomp./mean(denscomp); % normalize
%  denscomp=0.3+0.7*denscomp;
%  thesignal=thesignal.*denscomp'; % apply the filter
% plot3(kx,ky,abs(thesignal))
% drawnow
% pause(.2)
%% main
NNH=floor(NNN/2);NL=128;spec=zeros(64,64,128);
aten=exp(-1.5*(1:NNN)'/NNN);
for nx=1:64
    for ny=1:64
        basefid=aten.*exp(coex*(kx*XX(nx,ny)+ky*YY(nx,ny))).*thesignal';
        spec0=fft([basefid;zeros(NNN,1)]);
        spec(nx,ny,:)=[spec0(2*NNN-63:2*NNN);spec0(1:64)];
    end
end
figure(1);imagesc(abs(squeeze(spec(:,33,:))));colormap gray
    