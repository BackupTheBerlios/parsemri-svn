function se9(parameter_name,parameter_values,varargin)
%% run PARSE with several different values of a parameter
% examples: 
% se9;
% se9('freq_off',-52:2:100,'dataset','fourtube_se');
% se9('echo_time',1000:1000:11000,'dataset','fourtube_se');
% se9('dataset','fourtube_se');
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
% se9
% |- plotter
% |- make_nii, from mathworks user site.
% |- save_nii, from mathworks user site.
% |- dataset_load
% |- se_opts
% |- varparser
% |- load_procpar
% '- query_procpar

%% init some things...
tic % start the process timer
if exist('parameter_name','var') % a default value for when no sweep is run
    disp([parameter_name ' will be evaluated at ']);
    disp(parameter_values)
else
    parameter_name='foo';
    parameter_values=0;
end
try % flag to show progress or not with frequently updated figs
    show=varparser('show',varargin);
    if show
        disp('progress display is turned on.')
    else
        disp('progress display is turned off.')
    end
catch
    show=1;
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
%% select dataset =======================================================
% also contains some dataset specific defaults
settings={}; % init with empty cell array that dataset_load fills in
dataset_load % call the dataset and opt parameter loading script
%% make a place to save the data
savedir=['sweepresults_' parameter_name]; % where should the results go
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
figy=4;
yres=64;
yrng=1:yres;
for parameter_value=parameter_values % try each parameter value
    sweep_count=sweep_count+1;
    coefficientmaps=separse4(mr_data,show,parameter_name,parameter_value,varargin{:},settings{:});
    paramvolume(:,:,sweep_count)=coefficientmaps;
    elapsed_time=toc;
    save([savedir '/' parameter_name num2str(parameter_value)],'coefficientmaps') % save current result
    set(0,'CurrentFigure',2);
    for idx=1:figy
        subplot('Position',[(sweep_count-1)/figx (idx-1)/figy (1/figx) (1/figy)]) % plot them all in a row.
        imagesc(coefficientmaps(yrng+((idx-1)*yres),:)) % display the current result
        axis off
        colormap gray;
    end
    drawnow
end
disp([ 'elapsed time: ' num2str(elapsed_time)]);
%% save all the results as a single 3D nifti format volume for inspection
paramvolume(65:128,:,:)=paramvolume(65:128,:,:)*100000; % scaling factor because Mag image is small valued.
save paramvolume paramvolume
try
    nii=make_nii(paramvolume,[2 2 2]);
    save_nii(nii,[savedir '/sigvolume']);
catch
   disp('could not make the nifti volumes. Do you have the nii tools installed?')
end

%% functions ==================================================

function coefficientmaps = separse4(mr_data,show,varargin)
% coefficientmaps = separse2(dataset,show,varargin)
% example: cfm = separse2('fourtube',1)
% mark bolding mbolding@uab.edu 
% 29 Nov '06 msb
% 03 Jan '07 msb new trajectories
% 19 Jan '07 msb converting to SE
% 05 Feb '07 msb converting to SE second try, merging pre- and main code
% 19 Feb '07 msb is "try catch" crazy
% 22 Feb '07 msb is adding more parameters... SE third try rolling in
% functions and making p into a structure
% 
% original:
% D Twieg    UAB Biomedical Engineering
% uses hexagonal evaluation grid.   

%% constants ===========================================================
% also see the tweakers

plot_a=1; % figure number for plot
maxIter=200; % the maximum number of iterations
coex=-2i*pi;
sq2=sqrt(2);
twopi=2*pi;
% gam=4258.0;
% gmax=3.643;           % low-g
% fov = 12.8;           % field of view in cm calibration traj dependent
% nres = 64;            % nominal resolution for reconstruction

%% tweakers <<<<<<<<<<<<<<<<<<<<
trimoff=220; %  how many initial data points to discard
% trimoff=8280 % echo time for current spin echo data
shift_fft=0; % calibration v. signal time offset 
dc_offset=21+1i*32;
hammfilter=0; % hamming filter in k-space 1=on, 0=off

%% look for varargs, optional inputs to the function =====================
disp(['=== loading ' num2str(size(varargin,2)/2) ' optional settings ==='])
se_opts % call the option loading script to parse the varargin
disp('=== end of optional settings ===')

%% try to grab some values from "procpar" ===============================
% Varian stores all of the acquisition parameters in procpar.
% it contains offsets, FOV, thk, powers, etc.
try
    % using resto the resonance offset frequency as an example
    procpar_vars=load_procpar(mr_data);
    showprocvar('thk',procpar_vars)
    showprocvar('pss0',procpar_vars)
catch
    disp('could not read some values from the procpar file...')
end

%% bring all of the plots to the front at the start
if(show==1)         
    figure(plot_a)
end

%% load the calibration trajectory, kss =================================
phi=[];
ktraj=[];
try
    load(calibrationdata);  % brings in phi and ktraj
catch
    error('Could not load the calibration data. Booo.');
end
disp(['processing map ' num2str(whichfid)]);
% gmax=4.6833;          % hi-g maximum read gradient used. Not the Varian
% system maximum gradient!
Ta=.068;                % duration of acquisition in ms
% Ta=.025;              % post echo
% Ta=0.045;
% kf=nres/(2*fov);
% delt=1.0/(gam*gmax*fov); % theoretical time step 
delt=1.0/198511.166253; % defined from sw in procpar
numTimePts=floor(Ta/delt);     % number of acquisition points
kss=ktraj(shift_fft+trimoff+1:numTimePts+trimoff+shift_fft);  % k trajectory with lead in and lead out trimmed off
phi=real(phi)+imag(phi);
phi=phi(shift_fft+trimoff+1:numTimePts+trimoff+shift_fft);  % trim phi to the correct length and offset
kss=kss*1.61;           % empirical scaling

%% build Hamming filter to correct for sample spacing ====================
if hammfilter==1
    kr=abs(kss);                            % find the radius of each point
    kfe=max(kr);                            % kfe set as max radius
    kch=(0.54+0.46*cos(pi*kr/(sq2*kfe)));   % circular Hamming, circumscribed
    kss=kss.*kch;
end

%% lengths and stopping tolerances for progressive conjugate =============
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
%% makehex1; % builds numVox, xx, yy, XX, YY hexagonal evaluation grid
% makehex1.m
% generate hexagonal sampling grid
% D. Twieg
fov = 12.8;             % field of view in cm calibration traj dependent
nres = 64;              % nominal resolution for reconstruction
delxnom=fov/nres;
w=1.128/(2*delxnom);
delx=1.0/(sqrt(3)*w); % circumscribed hexagon; 
% gives minimum bias error
dely=delx*sqrt(3)/2;
xx=zeros(80,80);
yy=xx;
%% make hexagonal grid
for nc=1:80
    for nr=1:2:79
        yy(nr,nc)=(nr-41)*dely;  % y locations
        yy(nr+1,nc)=(nr-40)*dely;
        xx(nr,nc)=(nc-33)*delx; % x locations
        xx(nr+1,nc)=((nc-33)+0.5)*delx;
    end
end
%% count hex grid points and define arrays
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
[XX,YY]=meshgrid(((1:64)-33)*fov/64); % standard rectilinear grid for data display

numVox=nl-1; % numVox is the number of grid points
%% initialize parameter arrays and structures
Voxels      = 1:numVox;         % handy array of voxel indicies
p.Mperp     = zeros(1,numVox);  % parameter estimate array. initial mag
p.fRa       = zeros(1,numVox);  % parameter estimate array. freq and decay
p.Rb        = zeros(1,numVox);  % parameter to account for 2*R2 prime
p.echo_time = echo_time;
ExponentA   = zeros(size(xx));  % storage for hex grid representation
ExponentB   = zeros(size(xx));
Magnitude   = zeros(size(xx));  
CC          = zeros(1,maxIter); % array to store current cost result

%% load the data.=========================================================
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

%% set up psis ========================================================
% only need calculate psis once for a given traj and grid

try
    disp('using cached psis, delete the psis.mat file to recalculate');
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
    for nloc = Voxels % for each location...
        psis(nloc,:) = exp(coex*(kx(1:numTimePts)*xa(nl)+ky(1:numTimePts)*ya(nl))); % effect at each Voxel due to k-traj
        whichplot='psis';plotter;fprintf('%s',whichplot)
    end
    save('psis','psis')
end

disp('setup is finished. beginning search...');
%% prep done. main code.==================================================
old_time_length = 0; % the time length from previous iteration
iteration = 1; % record number of steps and estimations in parameter space
CCold = 1e+15; % last "current cost"
CCnew = 1e+14; % current "current cost"
for nlen = 1:NLENGTHS  % for each data length, this is the main code loop
    time_length = NLIST(nlen); 
    while (abs(CCold-CCnew) > FLIST(nlen)*CCold || time_length ~= old_time_length ) % while not converged or new data
        sigEstimate = estimateSignal(p,psis,time_length);    % estimate the signal
        sigdiff = (theSignal(1:time_length)-sigEstimate); % estimate the error 
        
        %---------------------parameter gradient evaluation----------------       
        % define direction of 1-D search
        if (time_length ~= old_time_length) %  more data, reset search direction
            g = calcGrad(psis,time_length,sigdiff,p);
            h = g;
            whichplot='gradients';plotter;fprintf('%s',whichplot)
        else
            gold = g; 
            g = calcGrad(psis,time_length,sigdiff,p);
            ga = max(0,((g-gold)*g')/(gold*gold'));
            h = g + ga*h;
            whichplot='gradients';plotter;fprintf('%s',whichplot)
        end
        delpf   = NewtonMeth(p,psis,time_length,sigEstimate,theSignal,h); % how far do we step?
        p.Mperp = p.Mperp + delpf*(h(Voxels));                  % take a step magnitude
        p.fRa   = p.fRa   + delpf*(h(Voxels+numVox));           % take a step freq and Ra
        p.Rb    = p.Rb    + real(delpf*(h(Voxels+numVox*2)));   % take a step Rb
        for nl = Voxels                               % read parameter values into the hex grid locations
            Magnitude(ia(nl),ja(nl)) = p.Mperp(nl);   % amplitude, complex magnitude image
            ExponentA(ia(nl),ja(nl)) = p.fRa(nl);     % freq and part of R that may reverse after echo?
            ExponentB(ia(nl),ja(nl)) = p.Rb(nl);      % R correction after echo, now we can find R2 and R2'
        end
        M0Rect  =griddata(yy,xx,abs(Magnitude),XX,YY,'cubic');   % resample M on a rect grid
        frRect  =griddata(yy,xx,imag(ExponentA)/(twopi*delt),XX,YY,'cubic');   % resample f on a rect grid
        R2Rect =griddata(yy,xx,real(-1*ExponentA)/delt,XX,YY,'cubic');  % resample R2 on a rect grid
        R2primeRect =griddata(yy,xx,real(-1*ExponentB)/delt,XX,YY,'cubic');
         
        %----------------------------------------
        sigEstimate = estimateSignal(p,psis,time_length);
        sigdif = theSignal(1:time_length)-sigEstimate;
        CC(iteration)= sigdif*sigdif'; 
        whichplot='progress';plotter;fprintf('%s',whichplot)
        old_time_length = time_length;
        CCold = CCnew;
        CCnew = CC(iteration);
        iteration = iteration + 1;
    end
end

%% all done, show results -----------------
coefficientmaps=[frRect;M0Rect;R2Rect;R2primeRect];
return % FIN! 

%% sub-functions ====================================
%% Newtons method
function delpf = NewtonMeth(p,psis,numTimePts,sesto,theSignal,h)
numVox=size(p.Mperp,2);
Voxels=1:numVox;
% delta = 1e-6;  % step size for gradient estimation. <<<
% delta = 1e-7;
delta = 1e-8;

sigdiff = sesto-theSignal(1:numTimePts);
f_0 = sigdiff*sigdiff';
p0.echo_time=p.echo_time;

p0.Mperp = p.Mperp + delta*(h(Voxels));         % M perp
p0.fRa   = p.fRa   + delta*(h(Voxels+numVox));  % f and R
p0.Rb    = p.Rb    + delta*(h(Voxels+numVox*2));% R prime
sigdiff  = estimateSignal(p0,psis,numTimePts) - theSignal(1:numTimePts);
f_a_p = sigdiff*sigdiff';

p0.Mperp = p.Mperp - delta*(h(Voxels));         % M perp
p0.fRa   = p.fRa   - delta*(h(Voxels+numVox));  % f and R
p0.Rb    = p.Rb    - delta*(h(Voxels+numVox*2));% R prime
sigdiff  = estimateSignal(p0,psis,numTimePts) - theSignal(1:numTimePts);
f_a_n = sigdiff*sigdiff';

f_1_0 = (f_a_p-f_a_n)/(2*delta);    % First derivative @ alpha = 0
f_1_a_p = (f_a_p-f_0)/delta;        % First derivative @ alpha = delta/2
f_1_a_n = (f_0-f_a_n)/delta;        % First derivative @ alpha = -delta/2
f_2_0 = (f_1_a_p-f_1_a_n)/delta;    % Second derivative @ alpha = 0

delpf = max(0,-f_1_0/f_2_0);        % min point

return

%% signal estimation
function sigEst = estimateSignal(p,psis,timeLength)

sigEst=zeros(1,timeLength);     % array to hold the signal
b=p.Mperp;                      % initial magnetization
preechodecay =exp(p.fRa - real(p.Rb));
postechodecay=exp(p.fRa + real(p.Rb)); % imag part p.fRa is freq, real parts are relaxation

sigEst(1)=b*psis(:,1);          % first time point
for t=2:timeLength              % iterate over time points
    if t < p.echo_time
        b=b.*(preechodecay);    % phase evolution and decay due to model
    else
        b=b.*(postechodecay);   % after the echo some of the phase changes reverse
    end
    sigEst(t)=b*psis(:,t);      % effect of motion through k-space
end

return

%% gradient calculation
function g = calcGrad(psis,timeLength,sigdiff,p)
numVox=size(p.fRa,2);
Voxels=1:numVox;

CC_conj = conj(-sigdiff);
fR_grad = zeros(numVox,1);          % freq and decay
Rp_grad = zeros(numVox,1);          % R prime
M       = transpose(p.Mperp);       % M perp
wf      = transpose(exp(p.fRa - real(p.Rb))); % before spin echo
wfprime = transpose(exp(p.fRa + real(p.Rb))); % after spin echo
wf_n    = ones(numVox,1);
M_grad  = psis(:,1)*CC_conj(1);
for t = 1:timeLength
    if t < p.echo_time
        wf_n = wf_n.*wf;
    else
        wf_n = wf_n.*wfprime;
    end
    temp    = wf_n.*psis(:,t)*CC_conj(t);
    M_grad  = M_grad + temp;
    fR_grad = fR_grad + (t-1)*temp;
    if t > p.echo_time
        Rp_grad = Rp_grad + (t-p.echo_time)*temp;
    end
end
fR_grad = M.*fR_grad;
Rp_grad = M.*Rp_grad;
g(Voxels)            = -conj(2*M_grad);
g(Voxels+numVox)     = -conj(2*fR_grad);
g(Voxels+numVox*2)   = -conj(2*Rp_grad);
return

%% display a value from procpar
function showprocvar(procvar,procpar_vars)
value=query_procpar(procvar,procpar_vars);
disp([procvar ' was '])
disp(value)
return
%% eq. 8.15, 8.18  from Haacke

%                                   { exp(-t*R2)*exp(-t*R2')      | 0    < t < tau
% Mtransverse(t) = Mtransverse(0) * { exp(-t*R2)*exp((t-TE)*R2')  | tau  < t < 2tau
%                                   { exp(-t*R2)*exp((TE-t)*R2')  | 2tau < t


