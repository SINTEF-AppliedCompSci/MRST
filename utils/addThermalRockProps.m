function rock = addThermalRockProps( rock, varargin )
%Add thermal properties to an existing rock structure 
%
% SYNOPSIS:
%   rock = addThermalRockProps(rock)
%   rock = addThermalRockProps(rock, 'pn1', pv1)
%
% PARAMETERS:
%   rock   - rock structure created with makeRock.
%  
%   lambdaR - Rock thermal conductivity in W m-1 K-1. typically ~2 
%             Can be either a single, scalar value of a column vector 
%             with one entry per cell. 
%
%   rhoR - rock density. Can be either a single, scalar
%          value of a column vector with one entry per cell. 
% 
%   CpR  - rock heat capacity in J kg-1 K-1. typically ~ 2000. Can be
%          either a single, scalar value of a column vector with one 
%          entry per cell.
%
%   tau  - rock tortuosity (no dimension) - used for the GeothermalWaterNaClModel
%          typically a scalar but can also be a vector column with one
%          entry per cell.
%
% RETURNS:
%   rock - valid rock structure with thermal properties for each active cell 
%          in the grid.
%
% EXAMPLE:
%
% SEE ALSO:
%   'makeThermalRock', 'makeRock'

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

    Watt = joule/second;
    opt = struct('lambdaR', 2*Watt/(meter*Kelvin)       , ... % Thermal conductivity
                 'rhoR'   , 2700*kilogram/meter^3       , ... % Density
                 'CpR'    , 1000*joule/(kilogram*Kelvin), ... % Specific heat capacity
                 'tau'    , 1                           );    % Tourtosity
    opt = merge_options(opt, varargin{:});
    % Thermal conductivity
    nc = length(rock.perm);
    if size(opt.lambdaR, 1) < nc
        opt.lambdaR = repmat(opt.lambdaR, nc, 1);
    end
    rock.lambdaR = opt.lambdaR;
    % Density
    if size(opt.rhoR, 1) < nc
        opt.rhoR = repmat(opt.rhoR, nc, 1);
    end
    rock.rhoR = opt.rhoR;
    % Specific heat capacity and internal energy
    if size(opt.CpR, 1) < nc
        opt.CpR = repmat(opt.CpR, nc, 1);
    end
    rock.CpR = opt.CpR;
    rock.uR  = @(p,T) rock.CpR.*T;
    % Tortuosity
    if size(opt.tau, 1) < nc
        opt.tau = repmat(opt.tau, nc, 1);
    end
    rock.tau = opt.tau;
    
end