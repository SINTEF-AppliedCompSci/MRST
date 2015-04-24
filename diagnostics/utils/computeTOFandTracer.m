function D = computeTOFandTracer(state, G, rock,  varargin)
%Compute time-of-flight and tracer distribution using finite-volume scheme.
%
% SYNOPSIS:
%    D  = computeTOFandTracer(state, G, rock)
%    D  = computeTOFandTracer(state, G, rock, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Construct the basis for flow diagnostic by computing
%     1) time-of-flight        :   \nabla·(v T) = \phi,
%     2) reverse time-of-flight:  -\nabla·(v T) = \phi,
%     3) stationary tracer     :  -\nabla·(v C) = 0
%   using a first-order finite-volume method with upwind flux.
%
% REQUIRED PARAMETERS:
%   G     - Grid structure.
%
%   rock  - Rock data structure.
%           Must contain a valid porosity field, 'rock.poro'.
%
%   state - Reservoir and well solution structure either properly
%           initialized from functions 'initResSol' and 'initWellSol'
%           respectively, or the results from a call to function
%           'solveIncompFlow'.  Must contain valid cell interface fluxes,
%           'state.flux'.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   wells - Well structure as defined by function 'addWell'.  May be empty
%           (i.e., wells = []) which is interpreted as a model without any
%           wells.
%
%   src   - Explicit source contributions as defined by function
%           'addSource'.  May be empty (i.e., src = []) which is
%           interpreted as a reservoir model without explicit sources.
%
%   bc    - Boundary condition structure as defined by function 'addBC'.
%           This structure accounts for all external boundary conditions
%           to the reservoir flow.  May be empty (i.e., bc = []) which is
%           interpreted as all external no-flow (homogeneous Neumann)
%           conditions.
% RETURNS:
%   D - struct that contains the basis for computing flow diagnostics:
%       'inj'     - list of injection wells
%       'prod'    - list of production wells
%       'tof'     - time-of-flight and reverse time-of-flight returned
%                   as an (G.cells.num x 2) vector
%       'itracer' - steady-state tracer distribution for injectors
%       'ipart'   - tracer partition for injectors
%       'ptracer' - steady-state tracer distribution for producers
%       'ppart'   - tracer partition for producers

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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


% Process optional parameters
opt = struct('bc', [], 'src', [], 'wells', []);
opt = merge_options(opt, varargin{:});

check_input(G, rock, opt);

state = validateStateForDiagnostics(state);

% Find injectors and producers
wflux = getTotalWellSolFlux(state.wellSol);
iwells = wflux > 0;
D.inj  = find( iwells);
D.prod = find(~iwells);

%Check that we actually have both injectors and producers
if (numel(D.inj) * numel(D.prod) == 0)
    D.tof = Inf(G.cells.num, 2);
    D.itracer = NaN(1, numel(D.inj));
    D.ipart = NaN;
    D.ptracer = NaN(1, numel(D.prod));
    D.ppart = NaN;
    return;
end

% Compute time-of-flight and tracer partition from injectors
t = computeTimeOfFlight(state, G, rock, 'wells', opt.wells, ...
   'tracer', {opt.wells(D.inj).cells});
D.tof     = t(:,1);
D.itracer = t(:,2:end);
[val,D.ipart] = max(D.itracer,[],2); %#ok<*ASGLU>

% Compute time-of-flight and tracer partition from producers
t = computeTimeOfFlight(state, G, rock, 'wells', opt.wells, ...
   'tracer', {opt.wells(D.prod).cells}, 'reverse', true);
D.tof(:,2) = t(:,1);
D.ptracer  = t(:,2:end);
[val,D.ppart] = max(D.ptracer,[],2);
end

%--------------------------------------------------------------------------

function check_input(G, rock, opt)
   assert(~isempty(opt.wells));
   assert(isempty(opt.src), 'Source terms not supported yet');
   assert(isempty(opt.bc),  'Boundary conditions not supported yet');

   assert (isfield(rock, 'poro')         && ...
           numel(rock.poro)==G.cells.num,   ...
          ['The rock input must have a field poro with porosity ',...
           'for each cell in the grid.']);

   assert(min(rock.poro) > 0, 'Rock porosities must be positive numbers.');
end

%--------------------------------------------------------------------------

function flux = getTotalWellSolFlux(wellSol)
    nw = numel(wellSol);
    flux = zeros(1, nw);
    for i = 1:nw
        ws = wellSol(i);
        f = sum(ws.flux);
        if isempty(f) || f == 0 && ~isempty(ws.sign)
            % If we have a sign for well and there is no flux defined, set
            % it to a small value of right sign instead.
            f = sqrt(eps)*ws.sign;
        end
        flux(i) = f;
    end
end
