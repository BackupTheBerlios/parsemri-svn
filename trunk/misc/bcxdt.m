% bcXdt
% process pinecone calibration data
% developed with remarkable Lamarkian genetic algorithms
% 11/27/06 msb v_v
% 11/28/06 msb -_-
%====================================================
% bcXdt.m
%====================================================
%% initialization
% clear
% figure;
% there are 8 sequence parameter combos at 8 locations
% start with numbers below, add multiples of 8 to generate other location
% numbers:
%      spin-echo: 5-8
% high-gradients: odd numbers
%        fat-sat: 3,4 7,8
%  circumscribed: all of them
combos=8; % see above 
% combo=2; % low-gradient
% combo=7; % the full monty
combo=1; % high-gradient
locationnums=[1 2 3 4 5 6 7 8]; % which spatial locations
%numlocations=size(locationnums,2);
alllocations=combo:combos:(max(locationnums)+1)*combos; % first number is the starting point
locations=alllocations(locationnums); % select the proper file numbers
% allxoffsets=[0.5 0.4 0.00 -0.2 -0.50  0.6  0.3].'; % spatial location of point phantoms
% allyoffsets=[0.5 0.7 0.75  0.8  0.45 -0.1 -0.35].';
allxoffsets=[-0.29; -1.65; -2.00;  -1.45;  1.00;  1.55;  1.15;  0.10];
allyoffsets=[ 1.85;  1.10;  0.10;  -1.15;  1.25;  0.30; -1.10; -1.65];
xoffset=allxoffsets(locationnums); % select the spatial locations
yoffset=allyoffsets(locationnums);
count=0;

%%  load and process signals, end up with unwrapped phase
for location=locations
    count=count+1;
    disp(['loading file ' num2str(location)]);
    if location < 10  % read in the data from disk
        [RE,IM]=LOAD_FID(['pinecone_0' num2str(location)]);
    else
        if location < 100
            [RE,IM]=LOAD_FID(['pinecone_' num2str(location)]);
        else
            [RE,IM]=LOAD_FID(['pinecone_' num2str(location)]); % 3 digits, the same for now
        end            
    end
    rawsig=RE+i*IM+(20+1i*32);   % compose signal and add offset
    sumsig=sum(rawsig,2);           % signal averaging
    filtsig(:,count)=sumsig;        % use raw signal average
    phase(:,count)=unwrap(angle(filtsig(:,count))); % unwrap the signal
end
clear count

%% take the difference of the unwrapped angle, take the spectrum
% of that, and simply suppress all frequencies beyond the 
% higher frequency of the signal, and use
% the inverse transform of this. --DT
phasediff=diff(phase); % find differences of phase

%% generate k-space trajectory
nfitlength=3500;
X=[xoffset yoffset]; % spatial locations
Y=phasediff(401:nfitlength+400,:);
nperp=200;nperm=800;

y=Y(:,1);
for n=2:size(Y,2);
y=[y;Y(:,n)];
end
%% the nlin fit, fit the traj model to the data
themodel=2;
switch themodel
    case 1
        %% model 1
        params0=[0.3688 -0.0800 0.3814 -0.0842 0.0334 0.0074 -3.9265 -0.8617 -0.7142 -.8554];
        % params0=[0.1 0.1 0.1 0.1 2*pi/nperp 2*pi/nperm 0 0 0 0];
        [params,r,J]=nlinfit(X,y,@ptcalmodel,params0);
        pa=params;a1=pa(1);a2=pa(2);a3=pa(3);a4=pa(4);omp=pa(5);omm=pa(6);pha1=pa(7);pha2=pa(8);pha3=pa(9);pha4=pa(10);
        for n=1:nfitlength
            g(n,1)=a1*sin(omp*n+pha1)+a2*sin(omm*n+pha2);
            g(n,2)=a3*cos(omp*n+pha3)+a4*cos(omm*n+pha4);
        end
    case 2
        %% model 2
        % this model includes the global phi value, assumes it is proportional to
        % gradient waveforms
        % the starting parameters: a parameter set from a previous fit:
        params0=[0.3688 -0.0800 0.3814 -0.0842 0.0334 0.0074 -3.9265 -0.8617 -0.7142 -.8554 0 0];
        % params0=[0.1 0.1 0.1 0.1 2*pi/nperp 2*pi/nperm 0 0 0 0 .001 0.001];
        [params,r,J]=nlinfit(X,y,@ptcalmodel2,params0);
        pa=params;
        a1=pa(1);
        a2=pa(2);
        a3=pa(3);
        a4=pa(4);
        omp=pa(5);
        omm=pa(6);
        pha1=pa(7);
        pha2=pa(8);
        pha3=pa(9);
        pha4=pa(10);
        nlength=17890-400;
        for n=1:nlength
            g(n,1)=a1*sin(omp*n+pha1)+a2*sin(omm*n+pha2);
            g(n,2)=a3*cos(omp*n+pha3)+a4*cos(omm*n+pha4);
            phi(n)=pa(11)*g(n,1)+pa(12)*g(n,2);
        end
    case 3
        %% model 3
        % this model includes the global phi value, assumes it is proportional to
        % gradient waveforms, plus an exponentially decaying component
        params0=[0.3688 -0.0800 0.3814 -0.0842 0.0334 0.0074 -3.9265 -0.8617 -0.7142 -.8554 0 0 0 0];
        % params0=[0.1 0.1 0.1 0.1 2*pi/nperp 2*pi/nperm 0 0 0 0 .001 0.001];
        [params,r,J]=nlinfit(X,y,@ptcalmodel3,params0);
        pa=params;
        a1=pa(1);
        a2=pa(2);
        a3=pa(3);
        a4=pa(4);
        omp=pa(5);
        omm=pa(6);
        pha1=pa(7);
        pha2=pa(8);
        pha3=pa(9);
        pha4=pa(10);
        for n=1:nfitlength
            g(n,1)=a1*sin(omp*n+pha1)+a2*sin(omm*n+pha2);
            g(n,2)=a3*cos(omp*n+pha3)+a4*cos(omm*n+pha4);
            phi(n)=pa(11)*g(n,1)+pa(12)*g(n,2)+pa(13)*exp(-pa(14)*n);
        end
end
%% separate the results
gam=4258;
twopi=2*pi;
fov=12.8; % field of view in cm
if mod(combo,2)==1 % odd numbers are high gradients
    gmax=4.683; % for high-gradient version
else
    gmax=3.643; % low gradient version
end
delt=1./(gam*gmax*fov);
krd0=1.644*(twopi*gam*cumsum(g(:,1))*delt+1.6); % scale the results
kpe0=1.644*(twopi*gam*cumsum(g(:,2))*delt+.115);
phi0=phi; 

%% show results
ktraj=krd0+i*kpe0;
figure;plot(ktraj)
title(num2str(combo));
drawnow;
axis image

%% save results
save( ['ktraj_pt_phan_model2_combo_' num2str(combo)] , 'ktraj','phi0')
