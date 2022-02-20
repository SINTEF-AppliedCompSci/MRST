function [p, s] = samplePiecewiseLinear(f, sInit, tolAbs)
% Utility to linearize smooth parametric curve f(s)
%
% SYNOPSIS:
%   [p, s] = samplePiecewiseLinear(f, sInit, tolAbs)
%
% DESCRIPTION:
%   This function picks points {s_i} along the curve segment 
%       f(s), sInit(1) <= s <= sInit(end)
%   such that for each output parameter interval [s_i, s_{i+1}] the distance 
%   from the line segment between f(s_i) and f(s_{i+1})) to the
%   corresponding curve segment is smaller that tolAbs. 
%   
%   The above requirement is ignored whenever norm( f(s_{i+1})-f(s_i) ) < tolAbs
%
% PARAMETERS:
%   f       - function handle for f:R -> R^n.
%   sInit   - Initial set of parameter values, must at least contain two
%             elements such that sInit(1) is minimal and sInit(end) is maximal 
%             paramter value. Any other values (1<i<end) is used as initial 
%             linearisation, i.e., s -> interp(sInit, f(sInit), s)
%   tolAbs  - Absolute tollerance, see description above 
%
% RETURNS:
%   p - points along the curve
%   s - parameters sich that p = f(s)
%
% EXAMPLE:
%   [p, s] = samplePiecewiseLinear(@(t)[t.*cos(t), t.*sin(t)], [0 4*pi], 0.1)
%   plot(p(:,1), p(:,2), 'o-')

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

% compute initial points
assert(numel(sInit)>1, 'Input must contain at least two parameter values (sInit)')
s = sInit(:);
p = f(s);

% use squared tollerance
tolSq = tolAbs^2;

% if p is one-dimensional, f is not parametric
if size(p, 2) == 1
    warning('Requires parametric function, setting f -> [t, f(t)]')
    f = @(t)[t, f(t)];
end

done = false;
ki   = 1;
opts = optimset('Display','off', 'TolFun', tolAbs/10);
warn = false;
while ~done
    % current points
    pcur = f(s([ki, ki+1],:));
    
    % if distance between points are less than tol, dont refine 
    if  norm(diff(pcur)) <= tolAbs
        dMaxSq = 0;
    else
        % find largest distance squared along segement
        [sOpt, dMaxSq, flag] = fminbnd(@(s)-distSq(f(s), pcur), s(ki), s(ki+1), opts);
         dMaxSq = -dMaxSq;
        if ~(flag == 1) % if something fails, just choose midpoint
           warn = true;
           sOpt   = (s(ki) + s(ki+1))/2;
           dMaxSq = distSq(f(sOpt), pcur);
        end
    end
    
    if dMaxSq >= tolSq
        % Add new parameter
        s = [s(1:ki); sOpt; s(ki+1:end)];
    else
        % move to next interval
        ki = ki +1;
        if ki >= numel(s)
            done = true;
        end
    end
end
% recalculate points
p = f(s);
if warn
    warning('fminbnd had problems finding a solution, examination of resulting trajectory is recommended');
end
end

% -------------------------------------------------------------------------
function d2 = distSq(p, l)
    % squared distance from point p to line segment between l(1,:) and l(2,:)
    p1 = l(1,:);
    v  = diff(l);  % should be non-zero here
    v2 = v*v';
    m  = p-p1;
    % t-value for p1 + t*v:
    t = (m*v')/v2;
    % distance to line
    d2 = m*m' -(t^2)*v2;
    % add squared of distance to segment along line
    d2 = d2 + (min(0, t)^2 + max(0, t-1)^2)*v2;
end
