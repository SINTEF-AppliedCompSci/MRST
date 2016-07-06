function v = volumesOfTraps(Gt, res, varargin)
% Compute volumes of (a subset of) precomputed structural traps
% 
% SYNOPSIS:
%   v = volumesOfTraps(Gt, res, varargin)
%
% PARAMETERS:
%   Gt       - top surface grid
%   res      - trap structure, as computed by trapAnalysis
%   varargin - may contain a vector with indices of traps for which to
%              compute volumes, or can be left empty, if volumes are to be
%              computed for _all_ traps
%
% RETURNS:
%   v - vector containing volumes for each specified trap
%       NB: The volumes computes are total volumes, not pore volumes!
%
% SEE ALSO:
%   trapAnalysis    
%
%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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

    if nargin == 2
        trap = 1:max(res.traps);
    else
        trap = varargin{1};
    end
    v = zeros(1,numel(trap));
    for i = 1:numel(trap)
        ind     = res.traps == trap(i);
        z       = Gt.cells.z(ind);
        fill_z  = res.trap_z(trap(i));

        assert(all((fill_z - z)>=0));

        H       = Gt.cells.H(ind);
        h_plume = min((fill_z - z),H);    
        v(i)    = sum(max(eps, Gt.cells.volumes(ind).*h_plume));
    end
end
