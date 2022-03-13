function masses = phaseMassesVEADI(Gt, state, rock, fluid)
% Compute column masses of undissolved gas, fluid, and dissolved gas.
%
% SYNOPSIS:
%   function masses = phaseMassesVEADI(Gt, sol, rock, fluidADI)
%
% DESCRIPTION:
%
% PARAMETERS:
%   Gt       - Top surface grid
%   sol      - solution structure containing a valid 'state' structure field
%   rock     - rock structure corresponding to 'Gt'
%   fluidADI - ADI fluid object.  Used to get densities and compressibilities
%
% RETURNS:
%   masses - Matrix with one row per cell and three columns.  The first
%            column gives the mass of undissolved gas per cell.  The second
%            column gives the mass of fluid (water/oil) per cell.  The third
%            column gives the mass of gas dissolved in fluid per cell.

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

    pvMult = 1; 
    if isfield(fluid, 'pvMultR')
        pvMult =  fluid.pvMultR(state.pressure);
    end
    
    % Defining convenience variables
    p      = state.pressure;
    sF     = state.s(:,1);  % fluid saturation
    sG     = state.s(:,2);  % gas saturation
    SF     = sF .* Gt.cells.H;   % vertically integrated saturation, fluid
    SG     = sG .* Gt.cells.H;   % vertically integrated saturation, gas
    pv     = rock.poro .* Gt.cells.volumes .* pvMult; % pore volumes
    rhoCO2 = fluid.rhoGS .* fluid.bG(p);
    rhoW   = fluid.rhoWS .* fluid.bW(p);

    % compute mass of undissolved gas
    gasPhase =  sum(pv .* rhoCO2 .* SG);
    
    % compute mass of liquid (water or oil)
    watPhase =  sum(pv .* rhoW .* SF);  
    
    % compute mass of dissolved gas
    if(isfield(state,'rs'))
        resDis = fluid.rhoGS .* sum(pv .* state.rs .* fluid.bW(p) .* SF);
    else
        resDis=0;
    end
    masses = [gasPhase,watPhase,resDis];  
end
