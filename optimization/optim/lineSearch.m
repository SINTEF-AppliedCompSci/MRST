function [u, v, g, info] = lineSearch(u0, v0, g0, d, f, c, opt)
% lineSearch -- helper function which performs line search based on
% Wolfe-conditions. 
% 
% SYNOPSIS:
% [u, v, g, info] = lineSearch(u0, v0, g0, d, f, c, opt)
%
% PARAMETERS:
% u0 :  current control vector
% v0 :  current objective value
% g0 :  current gradient
% d  :  search direction, returns an error unless d'*g0 > 0
% f  :  objective function s.t. [v,g] = f(u)
% c  :  struct with linear inequality and equality constraints
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

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
c1      = opt.wolfe1;
c2      = opt.wolfe2;
sgf     = opt.safeguardFac;
incTol  = opt.stepIncreaseTol;
maxIt   = opt.lineSearchMaxIt;

% Assert search direction is increasing
assert(d'*g0 >=0)

% Define convenience-function to gather info for a "point", where a is the 
% step-length, v is the value and dv is the dirctional derivative:
assignPoint = @(a,v,dv)struct('a', a, 'v', v, 'dv', dv);
p0          = assignPoint(0, v0, d'*g0);  % 0-point

% Function handles for wolfe conditions wrt p0
%w1 = @(p) p.v >= p0.v + c1*p.v*p0.dv;
w1 = @(p) p.v >= p0.v + c1*p.a*p0.dv;
w2 = @(p) abs(p.dv) <= c2*abs(p0.dv);

% Maximal step-length s.t. u = u0+aMax*d is feasible (should always be >= 1)
aMax = getAlphaMax(u0, d, c);
assert(aMax>1-sqrt(eps)); 
aMax = max(1, aMax);

% End-points of initial interval
p1 = p0;
p2 = assignPoint(aMax, -inf, -inf);
% Initial step-length:
%a  = 1;
a = min(1, aMax/2);

lineSearchDone   = false;
it = 0;
while ~lineSearchDone && it < maxIt
    it = it +1;
    u = u0 + a*d;
    [v, g]  = f(u); % function evaluation / simulation
    p = assignPoint(a, v, d'*g); % new point
    if w1(p) && w2(p)
        lineSearchDone = true;
        flag = 1;
    else
        if (p.a > aMax*(1-sqrt(eps))) && ... % max step-length reached
               (p.v > p0.v) && ... % the step yielded an improvement
               (p.dv > 0) % continuing further would improve (but not allowed)
            lineSearchDone  = true;
            flag = -1;
            fprintf('Wolfe conditions not satisfied for this step ...\n');
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
            % Find next candidate-step by interpolation
            a = argmaxCubic(p1, p2);
            % Safe-guarding and thresholding:
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
% Check if line search succeeded
if ~lineSearchDone
    flag = -2;
    fprintf('Line search unable to succeed in %d iterations ...\n', maxIt);
    % Although line search did not succeed in maxIt iterations, we ensure
    % to return the greater of p1 and p2's objective value none the less.
    if p1.v > p2.v
        v = p1.v;
        u = u0 + p1.a*d;
        % and we will re-compute the gradients for these controls
        [v, g]  = f(u);
    end
end
info = struct('flag', flag, 'step', a, 'nits', it);
end

function alphaMax = getAlphaMax(u, d, c)
% Find maximal a s.t. Ai*(u+a*d) <= bi
[A, b] = deal(c.i.A, c.i.b);

ignore_ix = false(numel(b),1);
isSmall = (abs(d)<sqrt(eps));
ignore_ix(1:2*numel(u)) = [isSmall; isSmall];

s = (b-A*u)./(A*d);
s = s(~ignore_ix);
alphaMax = min(s(s>eps));
end

   
    

    