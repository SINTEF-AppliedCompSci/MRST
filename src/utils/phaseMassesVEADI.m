function masses = phaseMassesVEADI(Gt, sol, rock, fluidADI)
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
%
    pvMult = 1; 
    if isfield(fluidADI, 'pvMultR')
        pvMult =  fluidADI.pvMultR(sol.state.pressure);
    end
    
    % Defining convenience variables
    p      = sol.state.pressure;
    sF     = sol.state.s(:,1);  % fluid saturation
    sG     = sol.state.s(:,2);  % gas saturation
    SF     = sF .* G.cells.H;   % vertically integrated saturation, fluid
    SG     = sG .* G.cells.H;   % vertically integrated saturation, gas
    pv     = rock.poro .* G.cells.volumes .* pvMult; % pore volumes
    rhoCO2 = fluidADI.rhoG .* fluidADI.bG(p);
    rhoW   = fluidADI.rhoO .* fluidADI.bO(p);

    % compute mass of undissolved gas
    gasPhase =  sum(pv .* rhoCO2 .* SG);
    
    % compute mass of liquid (water or oil)
    watPhase =  sum(pv .* rhoW .* SF);  
    
    % compute mass of dissolved gas
    if(isfield(sol.state,'rs'))
        resDis = fluidADI.rhoG .* sum(pv .* sol.state.rs .* fluidADI.bO(p) .* SF);
    else
        resDis=0;
    end
    masses = [gasPhase,watPhase,resDis];  
end
