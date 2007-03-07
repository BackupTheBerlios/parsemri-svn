%% options for separse
% place a reasonable default for every option or setting in this file
% use varparser to pull the variables out of varargin
%
try % frequency offset, this is used to untwist the signal.
    freq_off=varparser('freq_off',varargin);
    disp(['using frequency offset of ' num2str(freq_off)])
catch
    freq_off=0; % frequency offset generic frequency
    disp(['using a default frequency offset of ' num2str(freq_off)])
end

try % signal weighting for regularization of parameter space. It is multiplied times the signal.
    sig_weight=varparser('sig_weight',varargin);
    disp(['using sig_weight of ' num2str(sig_weight)])
catch
    sig_weight=25; % signal normalization factor
    disp(['using a default sig_weight of ' num2str(sig_weight)])
end

try % which traces do we want from the data file
    whichfid=varparser('whichfid',varargin);
    disp(['using whichfid of ' num2str(whichfid)])
catch
    whichfid=1;    %default signal to use from the data
    disp(['using the first signal trace in the data'])
end

try % which calibration trajectory should be used
    calibrationdata=varparser('caldata',varargin);
    disp(['using optional calibrationdata ' calibrationdata])
catch
    calibrationdata='circ_lowg_calibration_09jan2007/circ_lowg_calibration.mat'; % default 
    disp(['using default calibrationdata ' calibrationdata])
end

try % parse mode,  is this fid or spin echo data
    p.parse_mode=varparser('parse_mode',varargin);
catch
    p.parse_mode='fid'; % default 
end
disp(['using ' p.parse_mode ' signal mode'])

try % echo time,  is this fid or spin echo data?
    p.echo_time=varparser('echo_time',varargin);
catch
    if strcmp(p.parse_mode,'se') 
        p.echo_time=8182; % default for se determined from point phantom data.
    else
        p.echo_time=Inf; % default for FID data 
    end
end
disp(['using ' num2str(p.echo_time) ' echo time'])