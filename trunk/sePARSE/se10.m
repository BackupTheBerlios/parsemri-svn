function se10(parameter_name,parameter_values,varargin)
%% run PARSE with several different values of a parameter
% examples: 
% se10; % tests the code
% se10('freq_off',-52:2:100,'dataset','fourtube_se');
% se10('echo_time',1000:1000:11000,'dataset','fourtube_se');
% se10('dataset','fourtube_se');
% This is the currently the top level script for sePARSE code.
% add your datasets and specific parameters in dataset_load.m
% add new optional parameter try..catch blocks in se_opt.m
% Varargin is a matlab keyword that is used to pull in optional name, value
% pairs which are then parsed by the local function varparser
% separse4 does most of the work
% uses NIfTI toolbox from
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=8797&objectType=File
%
% msb 05 Feb '07
% msb 19 Feb '07
% msb 22 Feb '07
% 
% dependencies:
% se10
% |- plotter
% |- make_nii, from mathworks user site.
% |- save_nii, from mathworks user site.
% |- dataset_load
% |- se_opts
% |- varparser
% |- load_procpar
% '- query_procpar
seversion=10;
%% init some things...
tic % start the process timer
if exist('parameter_name','var') 
    disp([parameter_name ' will be evaluated at ']);
    disp(parameter_values)
else % if the code is called without parameters then execute test case
    disp('running sanity check.')
    parameter_name='test';
    parameter_values='0';
    s=['mv psis.mat psis' num2str(seversion) '.bak'];
    system(s);
end
try % flag to show progress or not with frequently updated figs
    plot_a.show=varparser('show',varargin);
    if plot_a.show
        disp('progress display is turned on.')
    else
        disp('progress display is turned off.')
    end
catch
    plot_a.show=1;
end
if strcmp(parameter_name,'dataset')
    dataset=parameter_values;
    disp(['analyzing ' dataset])
    parameter_name=parameter_values;
    parameter_values=1;
else
    try % choose the dataset using input value or default value...
        dataset=varparser('dataset',varargin);
        disp(['using dataset ' dataset])
    catch
        dataset='fourtube'; % an real example dataset to see if the code is working.
        disp(['using a default dataset of ' dataset])
    end
end
%% select dataset 
% also contains some dataset specific defaults
settings={}; % init with empty cell array that dataset_load fills in
dataset_load % call the dataset and opt parameter loading script
%% make a place to save the data
savedir=['sweepresults_' num2str(seversion) '_' parameter_name]; % where should the results go
try
    mkdirresult=mkdir( savedir );
    if mkdirresult==1
        disp(['running will write to directory ' savedir '. continuing in 5 sec.'])
        pause(5)
    end
catch
    error('couldn not make dir for results')
end 
%% run the sweep
sweep_count=0; % parameter value trial counter
figure(2);
figx=size(parameter_values,2); % one little figure for each trial
figy=4; % number of maps displayed per parameter value
for parameter_value=parameter_values % try each parameter value
    sweep_count=sweep_count+1;
    coefficientmaps=separse4(mr_data,plot_a,parameter_name,parameter_value,varargin{:},settings{:});
    yres=size(coefficientmaps,1)/figy;
    yrng=1:yres; 
    paramvolume(:,:,sweep_count)=coefficientmaps;
    save([savedir '/' parameter_name num2str(parameter_value)],'coefficientmaps') % save current result
    set(0,'CurrentFigure',2);  % select without bringing window forward.
    for idx=1:figy % display the images as we sweep the parameter
        plot_a.position=[(sweep_count-1)/figx (idx-1)/figy (1/figx) (1/figy)];
        subplot('Position',plot_a.position) % plot them all in a row.
        imagesc(coefficientmaps(yrng+((idx-1)*yres),:)) % display the current result
        axis off
        colormap gray;
    end
    drawnow    
    elapsed_time=toc;
    disp([ 'elapsed time: ' num2str(elapsed_time)])
end
%% save all the results as a single 3D nifti format volume for inspection
% paramvolume(65:128,:,:)=paramvolume(65:128,:,:)*100000; % scaling factor because Mag image is small valued.
save([savedir '/paramvolume'], 'paramvolume' );
try
    nii=make_nii(paramvolume,[2 2 2]);
    save_nii(nii,[savedir '/sigvolume.nii']);
catch
   disp('could not make the nifti volumes. Do you have the nii tools installed?')
end
for beeep=1:4
beep;
pause(0.6);
end
% done with the sweep

%% main functions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

function coefficientmaps = separse4(mr_data,plot_a,varargin)
% example: cfm = separse2('fourtube',1)
% mark bolding mbolding@uab.edu 
% 29 Nov '06 msb
% 03 Jan '07 msb new trajectories
% 19 Jan '07 msb converting to SE
% 05 Feb '07 msb converting to SE second try, merging pre- and main code
% 19 Feb '07 msb is "try catch" crazy
% 22 Feb '07 msb is adding more parameters... SE third try rolling in
% functions and making p into a structure
% 27 Feb '07 msb SE works I think, now trying to search from echo and not
% first data point.
% 02 Mar '07 msb refining the SE search. 
% 
% original:
% D Twieg    UAB Biomedical Engineering
% uses hexagonal evaluation grid.   

%% tweakers, infrequently changed options 
trimoff=220; %  how many initial data points to discard
% trimoff=8280 % echo time for current spin echo data
shift_fft=0; % calibration v. signal time offset 
dc_offset=21+1i*32;
hammfilter=1; % hamming filter in k-space 1=on, 0=off
maxIter=200; % the maximum number of iterations
plot_a.fig=1; % figure number for plot
% sq2=sqrt(2);
twopi=2*pi;
% gam=4258.0;
% gmax=3.643;           % low-g
fov = 12.8;             % field of view in cm, calibration traj dependent
nres = 64;              % nominal resolution for reconstruction
disp_res=72;            % display grid resolution
delxnom=fov/nres;       % the nominal x resolution for the hexagonal grid

%% look for varargs, optional inputs to the function 
disp(['=== loading ' num2str(size(varargin,2)/2) ' optional settings ==='])
se_opts % call the option loading script to parse the varargin
disp('=== end of optional settings ===')

%% try to grab some values from "procpar" 
% Varian stores all of the acquisition parameters in procpar.
% it contains offsets, FOV, thk, powers, etc.
try
    % using resto the resonance offset frequency as an example
    procpar_vars=load_procpar(mr_data);
    showprocvar('thk',procpar_vars);  % slice thickness
    showprocvar('pss0',procpar_vars); % slice offset
    sw=showprocvar('sw',procpar_vars);   % spectral width
    showprocvar('at',procpar_vars);   % spectral width
catch
    disp('could not read values from the procpar file...')
end

%% bring all of the plots to the front at the start
if(plot_a.show==1)         
    figure(plot_a.fig)
end

%% load the calibration trajectory, kss 
phi=[];
ktraj=[];
disp(['loading calibration data '  calibrationdata])
load(calibrationdata);  % brings in phi and ktraj
disp(['processing map ' num2str(whichfid)]);
% gmax=4.6833;          % hi-g maximum read gradient used. Not the Varian
% system maximum gradient!
Ta=.068;                % duration of acquisition in ms
% kf=nres/(2*fov);
% delt=1.0/(gam*gmax*fov); % theoretical time step 
% delt=1.0/198511.166253; % defined from sw in procpar
delt=1.0/double(sw); % defined from sw in procpar
numTimePts=floor(Ta/delt);     % number of acquisition points
kss=ktraj(shift_fft+trimoff+1:numTimePts+trimoff+shift_fft);  % k trajectory with lead in and lead out trimmed off
phi=real(phi)+imag(phi);
phi=phi(shift_fft+trimoff+1:numTimePts+trimoff+shift_fft);  % trim phi to the correct length and offset
kss=kss*1.61;           % empirical scaling 
if hammfilter==1        % hamming filter to correct for sampling density
    kss=hammit(kss);
end

%% lengths and stopping tolerances for progressive conjugate 
% gradient search.  For different rosettes, the lengths should be
% integer multiples of a single "echo time" -- i.e., the number of
% samples between successive crossings of the k-origin.
% You can determine this length by plotting the magnitude of the signal
% and counting the number of points between peaks.  The first number, N1, in
% NLIST should be an integer multiple of the echo time, and it should be
% sufficiently long to give a complete circle of rosette
% lobes, covering k-space more or less symmetrically.  Check this by
% plotting kss(1:N1). 
swoopsize=121.0;                        % one pass through the center of k space
N1=floor(5*swoopsize);                  % one "slice"
NLIST=[(N1*cumsum(ones(1,5)*4))';numTimePts];  % analyze a little data then a little more
NLENGTHS = size(NLIST,1);               % how many sections is it divided into
% FLIST=[0.3 ; 0.12 ; 0.1 ; 0.1 ; 0.1];   % tighten the tolerance as more data is added
FLIST=[0.2 ; 0.1 ; 0.05 ; 0.05 ; 0.05 ; 0.05]; % tighten the tolerance as more data is added
whichplot='slice';plotter;fprintf('%s',whichplot)
%% makehex1; % builds numVox, xx, yy, xa, ya, XX, YY hexagonal evaluation grid
% makehex1.m
% generate hexagonal sampling grid
% D. Twieg

%% make hexagonal grid
[xx,yy]=hexgrid80(delxnom);

%% count hex grid points and define arrays of locations
[nl,xa,ya,ia,ja]=locations80(xx,yy,fov); 
[XX,YY]=meshgrid(((1:disp_res)-floor(disp_res/2))*fov/disp_res); % standard rectilinear grid for data display

%% initialize parameter arrays and structures
numVox      = nl-1;             % numVox is the number of grid points
Voxels      = 1:numVox;         % handy array of voxel indicies
p.Mperp     = zeros(1,numVox);  % parameter estimate array. initial mag
p.fRa       = zeros(1,numVox);  % parameter estimate array. freq and decay
p.Rb        = zeros(1,numVox);  % parameter to account for 2*R2 prime
ExponentA   = zeros(size(xx));  % storage for hex grid representation
ExponentB   = zeros(size(xx));
Magnitude   = zeros(size(xx));  
CC          = zeros(1,maxIter); % array to store current cost result

%% decay model
p.echodecay.pre = @(p) exp(p.fRa - real(p.Rb));
p.echodecay.post= @(p) exp(p.fRa + real(p.Rb));

%% load the data.
disp(['loading ' mr_data]);
try
    for signal_index=whichfid % load the set of wanted signals
    [the_data_RE,the_data_IM]=LOAD_FID(mr_data,signal_index);
    the_data(:,signal_index)=the_data_RE+1i*the_data_IM + dc_offset; % combine re and im parts and dc offset
    end
catch
    error('could not load the data.')
end
theSignal=sum(the_data(trimoff+1:numTimePts+trimoff,:),2)'; % trim the signal to the proper size
theSignal=theSignal.*exp(-1i*phi)'; % phi correction 
theSignal=theSignal.*exp(-1i*twopi*freq_off*(1:numTimePts)*delt); % frequency offset 
theSignal=theSignal*sig_weight*numVox*delt/max(abs(theSignal)); % regularization
clear the_data_RE the_data_IM

%% Load psis, the basis functions representing the effect of k-traj on
% signal. only need calculate psis once for a given traj and grid
psis=psisGen(kss,numVox,numTimePts,xa,ya);


%% prep done. main code.
disp('setup is finished. beginning search...');
old_timeLength = 0; % the time length from previous iteration
iteration = 1; % record number of steps and estimations in parameter space
CCold = 1e+15; % last "current cost"
CCnew = 1e+14; % current "current cost"
for nlen = 1:NLENGTHS  % for each data length, this is the main code loop
    timeLength = NLIST(nlen); 
    timeVec=1:timeLength;
    while (abs(CCold-CCnew) > FLIST(nlen)*CCold || timeLength ~= old_timeLength ) % while not converged or new data
        % parameter gradient evaluation  xxxxxxxxxxxx      
        % define direction of 1-D search
        sigEstimate = estimateSignal(p,psis,timeVec);   % estimate the signal
        sigdiff = (theSignal(timeVec)-sigEstimate);     % estimate the error 
        if (timeLength ~= old_timeLength)               %  more data, reset search direction
            g = calcGrad(psis,timeVec,sigdiff,p);
            h = g;
            whichplot='gradients';plotter;fprintf('%s',whichplot)
        else
            gold = g; 
            g = calcGrad(psis,timeVec,sigdiff,p);
            ga = max(0,((g-gold)*g')/(gold*gold'));
            h = g + ga*h;
            whichplot='gradients';plotter;fprintf('%s',whichplot)
        end
        delpf   = NewtonMeth(p,psis,timeVec,sigEstimate,theSignal,h); % how far do we step?
        p.Mperp = p.Mperp + delpf*(h(Voxels));                  % take a step magnitude
        p.fRa   = p.fRa   + delpf*(h(Voxels+numVox));           % take a step freq and Ra
        p.Rb    = p.Rb    + real(delpf*(h(Voxels+numVox*2)));   % take a step Rb
        
        % scale and regrid
        for nl = Voxels                               % read parameter values into the hex grid locations
            Magnitude(ia(nl),ja(nl)) = p.Mperp(nl);   % amplitude, complex magnitude image
            ExponentA(ia(nl),ja(nl)) = p.fRa(nl);     % freq and part of R that may reverse after echo?
            ExponentB(ia(nl),ja(nl)) = p.Rb(nl);      % R correction after echo, now we can find R2 and R2'
        end
        M0Rect  = griddata(yy,xx,abs(Magnitude),XX,YY,'cubic');                  % resample M on a rect grid
        frRect  = griddata(yy,xx,imag(ExponentA)/(twopi*delt),XX,YY,'cubic');    % resample f on a rect grid
        R2Rect  = griddata(yy,xx,real(-1*ExponentA)/delt,XX,YY,'cubic');         % resample R2 on a rect grid
        R2primeRect = griddata(yy,xx,real(-1*ExponentB)/delt,XX,YY,'cubic');     % resample R2prime
         
        % measure, show, and record progress
        sigEstimate = estimateSignal(p,psis,timeLength);
        sigdif = theSignal(timeVec)-sigEstimate;
        CC(iteration)= sigdif*sigdif'; 
        whichplot='progress';plotter;fprintf('%s',whichplot)
        old_timeLength = timeLength;
        CCold = CCnew;
        CCnew = CC(iteration);
        iteration = iteration + 1;
    end
end

%% all done, show results 
coefficientmaps=[frRect;M0Rect;R2Rect;R2primeRect];
beep;
return % from separse4

%% sub-functions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%% Newtons method
function delpf = NewtonMeth(p,psis,timeVec,signalEstimate,theSignal,h)
numVox=size(p.Mperp,2);
Voxels=1:numVox;
% delta = 1e-6;  % step size for gradient estimation. <<<
% delta = 1e-7;
delta = 1e-8;
sigdiff = signalEstimate-theSignal(timeVec);
f_0 = sigdiff*sigdiff';
p0=p;

p0.Mperp = p.Mperp + delta*(h(Voxels));         % M perp
p0.fRa   = p.fRa   + delta*(h(Voxels+numVox));  % f and R
p0.Rb    = p.Rb    + delta*(h(Voxels+numVox*2));% R prime
sigdiff  = estimateSignal(p0,psis,timeVec) - theSignal(timeVec);
f_a_p = sigdiff*sigdiff';

p0.Mperp = p.Mperp - delta*(h(Voxels));         % M perp
p0.fRa   = p.fRa   - delta*(h(Voxels+numVox));  % f and R
p0.Rb    = p.Rb    - delta*(h(Voxels+numVox*2));% R prime
sigdiff  = estimateSignal(p0,psis,timeVec) - theSignal(timeVec);
f_a_n = sigdiff*sigdiff';

f_1_0 = (f_a_p-f_a_n)/(2*delta);    % First derivative @ alpha = 0
f_1_a_p = (f_a_p-f_0)/delta;        % First derivative @ alpha = delta/2
f_1_a_n = (f_0-f_a_n)/delta;        % First derivative @ alpha = -delta/2
f_2_0 = (f_1_a_p-f_1_a_n)/delta;    % Second derivative @ alpha = 0

delpf = max(0,-f_1_0/f_2_0);        % min point
% delpf = -(f_1_0/f_2_0);

return

%% signal estimation
function sigEst = estimateSignal(p,psis,timeVec)

sigEst=zeros(size(timeVec));    % array to hold the signal
b=p.Mperp;                      % initial magnetization
preechodecay = p.echodecay.pre(p); % imag part p.fRa is freq, real parts are relaxation
postechodecay= p.echodecay.post(p);

sigEst(1)=b*psis(:,1);          % first time point
for t=timeVec(2:end)            % iterate over time points
    if t < p.echo_time
        b=b.*(preechodecay);    % phase evolution and decay due to model
    else
        b=b.*(postechodecay);   % after the echo some of the phase changes reverse
    end
    sigEst(t)=b*psis(:,t);      % effect of motion through k-space
end

return

%% gradient calculation
function g = calcGrad(psis,timeVec,sigdiff,p)
numVox=size(p.fRa,2);
Voxels=1:numVox;

CC_conj = conj(-sigdiff);
fR_grad = zeros(numVox,1);          % freq and decay
Rp_grad = zeros(numVox,1);          % R prime
M       = transpose(p.Mperp);       % M perp
wfpre   = transpose(p.echodecay.pre(p)); % before spin echo
wfpost  = transpose(p.echodecay.post(p)); % after spin echo
wf_n    = ones(numVox,1);
M_grad  = psis(:,1)*CC_conj(1);
for t = timeVec
    if t < p.echo_time
        wf_n = wf_n.*wfpre;
    else
        wf_n = wf_n.*wfpost;
    end
    temp    = wf_n.*psis(:,t)*CC_conj(t);
    M_grad  = M_grad + temp;
    fR_grad = fR_grad + (t-1)*temp;
    Rp_grad = Rp_grad + (t-p.echo_time)*temp;
end
fR_grad = M.*fR_grad;
Rp_grad = M.*Rp_grad;
g(Voxels)            = -conj(2*M_grad);
g(Voxels+numVox)     = -conj(2*fR_grad);
g(Voxels+numVox*2)   = -conj(2*Rp_grad);
return

%% generate psis
function psis=psisGen(kss,numVox,numTimePts,xa,ya)
try
    disp('looking for cached psis, delete the psis.mat file to recalculate');
    disp('psis depends on the k traj and evaluation grid');
    load psis
catch
    disp('psis could not be loaded. regenerating psis...');
    kx=real(kss);
    ky=imag(kss);
    psis = zeros(numVox,numTimePts);
    % here psis is the array of ideal local responses to the encoding gradients for unit-magnitude on-resonance
    % magnetization vectors at each pixel (xa,ya) location.  These are "basis"
    % functions of a sort.
    figure(3)
    for nloc = 1:numVox % for each location...
        psis(nloc,:) = exp(-2i*pi*(kx(1:numTimePts)*xa(nloc)+ky(1:numTimePts)*ya(nloc))); % effect at each Voxel due to k-traj
        if mod(nloc,500) == 0  % show calculation progress once every 500 steps
            imagesc(real(psis));
            drawnow;
        end
    end
    fprintf('%s','saving psis ')
    save('psis','psis')
    disp('done')
end
return

%% hexgrid
function [xx,yy]=hexgrid80(delxnom)
w=1.128/(2*delxnom);
delx=1.0/(sqrt(3)*w); % circumscribed hexagon; 
% gives minimum bias error
dely=delx*sqrt(3)/2;
xx=zeros(80,80);
yy=xx;
for nc=1:80
    for nr=1:2:79
        yy(nr,nc)=(nr-41)*dely;  % y locations
        yy(nr+1,nc)=(nr-40)*dely;
        xx(nr,nc)=(nc-33)*delx; % x locations
        xx(nr+1,nc)=((nc-33)+0.5)*delx;
    end
end
return

%% locations
function [nl,xa,ya,ia,ja]=locations80(xx,yy,fov)
RA=1.0*fov/2; % radius
nl=1;
xa=zeros(80*80,1);
ya=xa;
ia=xa;
ja=xa;
for nc=1:80 % number of columns
    for nr=1:80 % number of rows
        if sqrt(xx(nr,nc).^2+yy(nr,nc).^2) < RA % if it is in the radius
            xa(nl)=xx(nr,nc);   % x coordinates
            ya(nl)=yy(nr,nc);   % y coordinates
            ia(nl)=nr; % row index
            ja(nl)=nc; % col index
            nl=nl+1;
        end
    end
end
return

%% hamming filter
function kss=hammit(kss)
kr=abs(kss);                            % find the radius of each point
kfe=max(kr);                            % kfe set as max radius
kch=(0.54+0.46*cos(pi*kr/(sqrt(2)*kfe)));   % circular Hamming, circumscribed
kss=kss.*kch;
return

%% display a value from procpar
function value=showprocvar(procvar,procpar_vars)
value=query_procpar(procvar,procpar_vars);
disp([procvar ' was '])
disp(value)
return

%% eq. 8.15, 8.18  from Haacke

%                                   { exp(-t*R2)*exp(-t*R2')      | 0    < t < tau
% Mtransverse(t) = Mtransverse(0) * { exp(-t*R2)*exp((t-TE)*R2')  | tau  < t < 2tau
%                                   { exp(-t*R2)*exp((TE-t)*R2')  | 2tau < t


