function state = validateStateForDiagnostics(state)
%Validate and fix state for flow diagnostics
%
% SYNOPSIS:
%   state = validateStateForDiagnostics(state)
%
% DESCRIPTION:
%   Routine for ensuring that states are suitable for diagnostic routines.
%   Not all solvers produce the "flux" fields directly, but may instead
%   give component fluxes that must be combined to form the flux. This
%   routine will, if possible, ensure that flux exists in the state.
%
% PARAMETERS:
%   state - State produced by one of MRST's solvers.
%
% RETURNS:
%   state - State with flux at wellSol and interfaces, if possible.
%
% NOTE:
%   This routine will throw an error if there is insufficient data to
%   produce the correct fields.
%
% SEE ALSO:
%   `computeTOFandTracer`.

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

% Check whether the state object has total flux and if necessary construct
% this flux from black-oil or similar type of model
if isfield(state, 'flux')
   state.flux = sum(state.flux,2);
else
   error('Reservoir state must provide total Darcy flux');
end
if  isfield(state, 'wellSol') 
    if isfield(state.wellSol, 'flux')
       for i=1:numel(state.wellSol)
          state.wellSol(i).flux = sum(state.wellSol(i).flux, 2);
       end
    elseif isfield(state.wellSol, 'cqs')
       for i=1:numel(state.wellSol)
          state.wellSol(i).flux = sum(state.wellSol(i).cqs, 2);
       end
   end
end
