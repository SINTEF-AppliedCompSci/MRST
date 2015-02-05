function [u, v, g, info] = lineSearch(u0, v0, g0, d, f, opt)
% lineSearch -- helper function which performs line search based on
% Wolfe-conditions. 
% 
% SYNOPSIS:
% [u, v, g, info] = lineSearch(u0, v0, g0, d, f, opt)
%
% PARAMETERS:
% u0 :  current control vector
% v0 :  current objective value
% g0 :  current gradient
% d  :  search direction, returns an error unless d'*g0 > 0
% f  :  objective function s.t. [v,g] = f(u)
% opt:  struct of parameters for line search
%
% RETURNS:
% u    : control vector found along d
% v    : objective value for u
% g    : gradient for u
% info :  struct containing fields:
%   flag : 1 for success (Wolfe conditions satisfied)
%          -1 max step-length reached without satisfying Wolfe-conds
%          -2 line search failed in maxIt iterations (returns most recent 
%          control without further considerations)
%   step : step-length
%   nits : number of iterations used
c1      = opt.wolfe1;
c2      = opt.wolfe2;
sgf     = opt.safeguardFac;
incTol  = opt.stepIncreaseTol;
maxIt   = opt.lineSearchMaxIt;

assert(d'*g0 >=0)
% convenience-function to gather info for a "point" in a struct, where
% a is the step-length, v is the value and dv is the dirctional derivative
% along d
assignPoint = @(a,v,dv)struct('a', a, 'v', v, 'dv', dv);
p0          = assignPoint(0, v0, d'*g0);  % 0-point

% function for wolfe conditions wrt p0
w1 = @(p) p.v >= p0.v + c1*p.v*p0.dv;
w2 = @(p) abs(p.dv) <= c2*abs(p0.dv);

% maximal step-length s.t. u = u0+aMax*d is feasible, i.e. s.t. 0 <= u <= 1
aMax = getAlphaMax(u0, d);
assert(aMax>=1-sqrt(eps));

% end-points of initial interval
p1 = p0;
p2 = assignPoint(max(aMax, 1), -inf, -inf);
% initial step-length:
a  = 1;

lineSearchDone   = false;
it = 0;
% XXX figure(2), clf
best = struct('u', u0, 'v', v0);
while ~lineSearchDone && it < maxIt
    it = it +1;
    u = u0 + a*d;
    [v, g]  = f(u); % function evaluation / simulation
    if v > best.v, best = struct('u', u, 'v', v); end
    p = assignPoint(a, v, d'*g); % new point
    if w1(p) && w2(p)
        lineSearchDone = true;
        flag = 1;
    else
        if (p.a > aMax*(1-sqrt(eps))) && (p.dv > 0) % max step-length reached
            if p.v > p1.v
                lineSearchDone  = true;
                flag = -1;
                fprintf('Wolfe conditions not satisfied for this step ...\n');
            else % oscillating derivative, half interval
                a = (p1.a+p2.a)/2;
                fprintf('Oscillating derivative, cutting interval in half ...\n')
            end
        else % logic for refining/expanding interval of interest [p1 p2]
            if p.a > p2.a
                if w1(p2)
                    p1 = p2;
                end
                p2 = p;
            else
                if ~w1(p) || p2.dv < 0
                    p2 = p;
                else
                    p1 = p;
                end
            end
            % find next candidate-step by interpolation
            a = argmaxCubic(p1, p2);
            % safe-guarding and thresholding:
            if a > p2.a
                a = max(a, (1+sgf)*p2.a);
                a = min(a, min(incTol*p2.a, aMax));
            elseif a > p1.a
                a = max(a, p1.a + sgf*(p2.a-p1.a));
                a = min(a, p2.a - sgf*(p2.a-p1.a));
            else
                a = (p1.a+p2.a)/2;
                fprintf('Cubic interpolation failed, cutting interval in half ...')
            end
        end
    end
end
% check if line search succeeded
if ~lineSearchDone
    flag = -2;
    fprintf('Line search unable to succeed in %d iterations ...\n', maxIt);
    [v, u] = deal(best.v, best.u);
end
info = struct('flag', flag, 'step', a, 'nits', it);
end

function alphaMax = getAlphaMax(u, d)
% find maximal a s.t. 0 <= u+a*d <= 1
step0 = inf(numel(u),1);
ix0   = and(abs(d)>eps, d<0);
step0(ix0) = -u(ix0)./d(ix0);

step1 = inf(numel(u),1);
ix1   = and(abs(d)>eps, d>0);
step1(ix1) = (1-u(ix1))./d(ix1);

alphaMax = min(min(step0), min(step1));
end

   
    

    