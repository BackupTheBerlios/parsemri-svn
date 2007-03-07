function value=varparser(name,pv_pairs)
% varparser.m
% msb Feb 2007
% Searches through varargins for a variable name and returns
% the first instance of the value. 
% Throws an error if it cannot find the var.
% example:
% value = varparser('foo',varargin);
%%
n=length(pv_pairs)/2;
if n ~= floor(n)
        error 'Property/value pairs must come in pairs.'
end
for i=1:n
    p_i = pv_pairs{2*i-1};
    v_i = pv_pairs{2*i};
    if strcmp(name,p_i)
        value=v_i;
        return
    end
end

error(['varparser could not find the parameter ' name])
