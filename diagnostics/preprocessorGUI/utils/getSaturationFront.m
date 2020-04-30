function [tau, sw, fsw, cumsw] = getSaturationFront(fluid, varargin)
% Compute simple Buckley-Leverett profile at 1 volume unit injected
%
% SYNOPSIS:
%    [tau, s, fs, cums]  = getSaturationFront(fluid)
%    [tau, s, fs, cums]  = getSaturationFront(fluid, 'pn1', pv1, ...)
% 
% REQUIRED PARAMETERS:
%  fluid - fluid structure containing functions krW(sw), krO(so), muW(p) and
%          muO(p).
%
%OPTIONAL PARAMETERS
%  nPoints - number of descrete points in output-vectors. Distributed
%            equally between tau = 0 and tau(end) which is the last
%            tau-value for which s>0 (deafult 250).
%  pRef    - reference pressure for computing viscosity values 
%            (default 200*barsa)
%  sw0     - initial (constant) water satuaration (deafult is connate water 
%            saturation) 
% 
% RETURNS:
%  tau     - time (of flight) vector
%  sw      - water satuation values corr to tau
%  fsw     - water fractional flow values corr to tau
%  cumsw   - cummulative water volume corr to tau (cumsw(end)=1)

opt = struct('nPoints', 250, 'pRef', 200*barsa, 'sw0', []);
opt = merge_options(opt, varargin{:});
% we do computations on a fine resolution, resample at end
np = 10000;

sw  = initVariablesAD_diagonal((np:-1:0)'/np); % negative direction

% compute frac-flow
mobw  = fluid.krW(sw)/fluid.muW(opt.pRef);
if isfield(fluid, 'krO')
    krO = fluid.krO;
else
    krO = fluid.krOW;
end
mobo  = krO(1-sw)/fluid.muO(opt.pRef);
fsw = mobw./(mobw + mobo);
% tof is dfsw/dsw
tau = fsw.jac{1}.diagonal;
% set sw to double
[sw, fsw] = deal(value(sw), value(fsw));

% get connate water saturation (remember sw from 1 down to 0)
conix = find(value(mobw)==0, 1, 'first');
swcon = sw(conix);
% work with saturation increase
sd   = sw-swcon;
% add zero unless present
if tau(1) > 0
    [tau, sd] = deal([0; tau], [sd(1); sd]);
end
% integrate water saturation increase (sd*dtau)
cumsw = [0; cumsum(sd(1:end-1).*diff(tau))];
% find index where integral equals one
ix   = find(cumsw>=1-sqrt(eps), 1, 'first');
if isempty(ix)
    lv = cumsw(end);
    ix = find(cumsw>=lv-1e-3, 1, 'first');
    warning('May be inaccurate, error %2.2e', 1-cumsw(ix));
end
% remove multiple values where tau == 0
ix1 = find(tau>0, 1, 'first');
ix1 = min(ix, max(1, ix1-1));
[tau, sw, fsw, cumsw] = deal(tau(ix1:ix), sw(ix1:ix), fsw(ix1:ix), cumsw(ix1:ix));

% if s0>swcon we need to rescale time and cums
if ~isempty(opt.sw0) && opt.sw0>swcon
    ii = find(sw<=opt.sw0, 1, 'first');
    if ~isempty(ii)
        [tau, sw, fsw, cumsw] = deal(tau(1:ii), sw(1:ii), fsw(1:ii), cumsw(1:ii));
    end
    truecums = cumsw(end) - (opt.sw0-swcon)*tau(end);
    % rescale
    cumsw = cumsw/cumsw(end);
    tau  = tau/truecums;
end

% check for repeated tau-values
[~, ixu] = uniquetol(tau);
% for linear fw we may get only one point at random position
if numel(ixu) ==1
    [tau, sw, fsw, cumsw] = deal([0; tau(end)], sw([1, end]), fsw([1, end]), [0; cumsw(end)]);
else
    [tau, sw, fsw, cumsw] = deal(tau(ixu), sw(ixu), fsw(ixu), cumsw(ixu));
end
assert(all(diff(tau)>0), 'Something went wrong :(');

% adjust last point such that cums(end) = 1 exactly
a = (1-cumsw(end-1))/(cumsw(end)-cumsw(end-1));
cumsw(end) = 1;
tau(end) = tau(end-1) + a*(tau(end)-tau(end-1));
sw(end) = sw(end-1) + a*(sw(end)-sw(end-1));
fsw(end) = fsw(end-1) + a*(fsw(end)-fsw(end-1));

% resample:
Fs  = griddedInterpolant(tau, sw);
Ffs = griddedInterpolant(tau, fsw);
Fcs = griddedInterpolant(tau, cumsw);

np  = opt.nPoints;
tau  = (0:np)'*tau(end)/np;
sw    = Fs(tau);
fsw   = Ffs(tau);
cumsw = Fcs(tau);
end

