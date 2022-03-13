function resSol = initResSolVE_s(G, p0, s0, varargin)
% Wrapper for initResSol which adds any extra properties needed by the
% vertical-equil module solvers. See resSol for details of valid
% arguments.
%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    % G = varargin{1};

    %opt = struct('use_s_form', false);
    %opt = merge_options(opt, varargin{:});

    % Create initial reservoir solution
    resSol = initResSol(G, p0, s0, varargin{:});
    % Height of the plume is assumed to be zero at t=0    
    resSol.extSat= repmat(resSol.s, 1, 2);
    resSol.h_max = 0*resSol.s;
    
    if any(strcmp(G.type, 'topSurfaceGrid'))
        resSol.twophaseJacobian = @twophaseJacobianWithVE_s;
    else
        % TODO: figure out region3D option
        resSol.twophaseJacobian = @twophaseJacobianWithVE_coupled;
    end
    
end
