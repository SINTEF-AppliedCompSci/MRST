function x_final = criticalPointChop(x0, x_final, xc, e, type, dir)
%Perform one or two-sided stability chop for an updated value
%
% SYNOPSIS:
%   x = criticalPointChop(x, x+dx, xc, e)
%
% REQUIRED PARAMETERS:
%   x0       - Quantity before update.
%
%   x        - Quantity after update
%
%   xc      - Critical point
%
%   e       - Epsilon used for stability
%
%   type    - If provided, indicates the type of chop: Above will chop so
%             that the passed value is on the other side of the critical
%             point, below the opposite and both will chop in both
%             directions. Default: Both
%
%   dir     - Direction that will be chopped. Either any, increasing or
%             decreasing. Default: Any, chop both increasing and decreasing
%             quantities.
% RETURNS:
%   x_final - Chopped value.
%

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

    assert(e > 0);
    if nargin < 6
        dir = 'any';
        if nargin < 5
            type = 'both';
            if nargin < 4
                e = 1e-8;
            end
        end
    end
    if numel(xc) == 1
        xc = repmat(xc, size(x0));
    end
    % x0: before
    % x: after
    % xc: critical point
    % e: Epsilon around critical point
    % type: below (chop point just below), above (chop just above) and both
    % dir: increasing, decreasing, both
    inc = x0 < xc & x_final > xc; % Increasing, crossing critical point
    dec = x0 > xc & x_final < xc; % Decreasing, crossing critical point
    switch lower(dir)
        case 'any'
            active = inc | dec;
        case 'increasing'
            active = inc;
        case 'decreasing'
            active = dec;
        otherwise
            error('Direction ''%s'' is not known. Supported values: ''any'', ''increasing'', ''decreasing''.', dir)
    end
    x = x_final(active);
    x0 = x0(active);
    xc = xc(active);
    dec = dec(active);
    inc = inc(active);
    % x_0       x_c - e | x_c | x_c + e       x_next
    % If both: first chop to x_c - e/2, then to x_c + 3*e/2
    %          below critical - into critical - above critical
    % If below: chop to x_c - e/2
    % if above: chop to x_c + e/2
    % Reverse definitions when coming from above
    switch lower(type)
        case 'both'
            % Define middle of interval
            middle = x0 < xc + e & x0 > xc - e;
            % Deal with increasing
            lo = inc & ~middle;
            hi = inc & middle;
            % If already inside critical region  -> chop to just outside interval
            x(hi) = min(x(hi), xc(hi) + 3*e/2);
            % If below critical region -> chop to critical region, lower end
            x(lo) = xc(lo) - e/2;

            % Deal with decreasing case (most defs. are flipped)
            lo = dec & middle;
            hi = dec & ~middle;
            x(lo) = max(x(lo), xc(lo) - 3*e/2);
            x(hi) = xc(hi) + e/2;
        case 'above'
            % If we are crossing from below, we chop to just above
            x(inc) = min(x(inc), xc(inc) + e);
            % Descending variant
            x(dec) = max(x(dec), xc(dec) - e);
        case 'below'
            % If we are crossing ascending, chop to just below the critical
            % point unless we are already sufficiently closed
            inc_low = inc & x0 < xc - e;
            x(inc_low) = min(x(inc_low), xc(inc_low) - e/2);
            inc_hi = dec & x0 > xc + e;
            x(inc_hi) = max(x(inc_hi), xc(inc_hi) + e/2);
        otherwise
            error('Type ''%s'' is not known. Supported values: ''both'', ''above'', ''below''.', type)
    end
    x_final(active) = x;
end