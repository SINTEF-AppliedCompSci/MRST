function [F,Phi] = computeFandPhi(arg1, varargin)
%Compute flow-capacity/storage-capacity diagram (F,Phi)
%
% SYNOPSIS:
%   [F,Phi] = computeFandPhi(pv, tof)
%   [F,Phi] = computeFandPhi(rtd)
%   [F,Phi] = computeFandPhi(rtd, 'sum', true/false);
%
% DESCRIPTION:
%   Compute the flow-capacity/storage-capacity based upon either:
%    1. Pore volumes and time-of-flight values per cell. In this case,
%       Phi is defined as the cumulative pore volume and F is the
%       cumulative flux as function of residence time. Assuming
%       incompressible flow, the instant flux is computed as the ratio of
%       pore volume and the residence time.
%    2. Residence-time distribution. In this case, F is the integral of the
%       RTD values, whereas Phi is the first-order moment (i.e., the
%       integral of the product of time and RTD values).
%   Making an analogue to 1D displacement theory, the F-Phi curve is the
%   equivalent to a plot of the fractional flow versus saturation.
%   Technical description: see Section 13.2-13.3 in the MRST book, Shavali
%   et al. (SPE 146446), and Shook and Mitchell (SPE 124625)
%
% PARAMETERS:
%   vol - pore volume, one value per cell
%   tof - two-component vector with time-of-flight from injector and
%         time-of-flight from producer, one value per cell
%   rtd - structure with residence-time distribution from computeRTD
%
% RETURNS:
%   F   - flow capacity = cumulative flux for increasing time-of-flight
%         values, where the flux per cell is defined from the relation
%                volume = flux * total_travel_time
%   Phi - storage capacity = cumulative pore volume for increasing
%         time-of-flight values
%
% SEE ALSO:
%   `computeTOFandTracer`, `computeLorenz`, `computeSweep`

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

if isstruct(arg1)
    [F,Phi] = computeFandPhiFromDist(arg1, varargin{:});
    return
else
    pv = arg1; tof = varargin{1};
end

assert(size(tof,2) == 2, 'Tof input must have two columns.')

t      = sum(tof,2);     % total travel time
[ts,order] = sort(t);    % sort cells base on travel time
v      = pv(order);      % put volumes in order
Phi    = cumsum(v);      % cumulative sum
vt     = full(Phi(end)); % total volume of region
Phi    = [0; Phi/vt];    % normalize to units of pore volumes
flux   = v./ts;          % back out flux based on incompressible flow
ff     = cumsum(flux);   % cumulative sum
ft     = full(ff(end));  % total flux computed
F      = [0; ff/ft];     % normalize and store flux
end

