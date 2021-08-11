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
%     3) stationary tracer     :  ±\nabla·(v C) = 0
%   using a first-order finite-volume method with upwind flux. A majority
%   vote is also used to partition the volume and assign each cell to a
%   unique tracer. Optionally, the routine can also compute time-of-flight
%   values for each influence region by solving localized time-of-flight
%   equations
%         \nabla·(v C_i T) = \phi C_i,
%   where C_i is the tracer concentration of each influence region.
%
%   Optionally, first arrival time is computed by a graph algorithm.
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
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs):
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
%
%   tracerWells - Logical index vector indicating subset of wells for which
%           tracer-fields will be computed. Empty matrix (default) is
%           interpeted as all true, i.e., compute tracer fields for all
%           wells.
%
%   computeWellTOFs - Boolean variable. If true, time-of-flight values are
%           computed individually for each influence region by solving
%           localized time-of-flight equations.
%
%   firstArrival - Boolean variable. If true, compute first-arrival time by
%           a graph algorithm.
%
%   solver - Function handle to solver for use in TOF/tracer equations.
%           Default (empty) is matlab mldivide (i.e., \)
%
%   maxTOF - Maximal TOF thresholding to avoid singular/ill-conditioned
%           systems. Default (empty) is 50*PVI (pore-volumes-injected).
%
%   processCycles - Extend TOF thresholding to strongly connected
%           components in flux graph by considering Dulmage-Mendelsohn
%           decomposition (dmperm). Recommended for highly cyclic flux
%           fields. Only takes effect if 'allowInf' is set to false.
%
% RETURNS:
%   D - struct that contains the basis for computing flow diagnostics:
%       'inj'     - list of injection wells
%       'prod'    - list of production wells
%       'tof'     - time-of-flight and reverse time-of-flight returned
%                   as an array where the first (G.cells.num x 2) elements
%                   contain forward/backward TOF for the whole field, and
%                   any other columns optinally contain TOF values for
%                   individual influence regions
%       'itracer' - steady-state tracer distribution for injectors
%       'ipart'   - tracer partition for injectors
%       'ifa'     - first-arrival time injectors
%       'ptracer' - steady-state tracer distribution for producers
%       'ppart'   - tracer partition for producers
%       'pfa'     - first-arrival time for producers

%{
Copyright 2009-2018 SINTEF Digital, Applied Mathematics.

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
opt = struct('bc', [],                ...
             'src', [],               ...
             'wells', [],             ...
             'tracerWells', [],       ...
             'solver', [],            ...
             'maxTOF', [],            ...
             'processCycles', false,  ...
             'computeWellTOFs', false, ...
             'firstArrival', false);
opt = merge_options(opt, varargin{:});
if opt.firstArrival
  opt.computeWellTOFs = true;
end

check_input(G, rock, opt);

state = validateStateForDiagnostics(state);
% If opt.tracerWells is empty include all
if isempty(opt.tracerWells)
    opt.tracerWells = true(numel(opt.wells), 1);
end
opt.tracerWells = opt.tracerWells(:)';

% Find injectors and producers
wflux = getTotalWellSolFlux(state.wellSol);
iwells = wflux > 0;
D.inj  = find( iwells & opt.tracerWells);
D.prod = find(~iwells & opt.tracerWells);

% Check that we actually have a meaningful TOF scenario. If all wells are
% shut off, return Inf.
sum_flux = cell2mat(arrayfun(@(x) sum(abs(x.flux)), state.wellSol, 'UniformOutput', false));
if (sum(sum_flux) == 0.0)
    D.tof = Inf(G.cells.num, 2);
    D.itracer = NaN(G.cells.num, numel(D.inj));
    D.ipart = NaN(G.cells.num, 1);
    D.ptracer = NaN(G.cells.num, numel(D.prod));
    D.ppart = NaN(G.cells.num, 1);
    return;
end

% Compute time-of-flight and tracer partition from injectors
t = computeTimeOfFlight(state, G, rock, 'wells', opt.wells,  ...
   'tracer', {opt.wells(D.inj).cells}, 'solver', opt.solver, ...
   'maxTOF', opt.maxTOF, 'processCycles', opt.processCycles, ...
   'computeWellTOFs', opt.computeWellTOFs,                   ...
   'firstArrival', opt.firstArrival);
D.tof     = t(:,1);
D.itracer = t(:,2:numel(D.inj)+1);
if opt.computeWellTOFs
    D.itof = t(:, numel(D.inj)+2:2*numel(D.inj)+1);
    if opt.firstArrival
        D.ifa = t(:, 2*numel(D.inj)+2:end);
    end
end
[val,D.ipart] = max(D.itracer,[],2); %#ok<*ASGLU>
% set 'non-traced' cells to zero
D.ipart(val==0) = 0;

% Compute time-of-flight and tracer partition from producers
t = computeTimeOfFlight(state, G, rock, 'wells', opt.wells, ...
   'tracer', {opt.wells(D.prod).cells}, 'reverse', true, ...
   'solver', opt.solver, 'maxTOF', opt.maxTOF, ...
   'processCycles', opt.processCycles, ...
   'computeWellTOFs', opt.computeWellTOFs, ...
   'firstArrival', opt.firstArrival);
D.tof(:,2) = t(:,1);
D.ptracer  = t(:,2:numel(D.prod)+1);
if opt.computeWellTOFs
    D.ptof = t(:, numel(D.prod)+2:2*numel(D.prod)+1);
    if opt.firstArrival
       D.pfa = t(:, 2*numel(D.prod)+2 : end);
    end
end
[val,D.ppart] = max(D.ptracer,[],2);
D.ppart(val==0) = 0;
end

%--------------------------------------------------------------------------

function check_input(G, rock, opt)
   assert(~isempty(opt.wells), 'Wells are required for computeTOFandTracer');
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
        % if status is false, set flux to zero
        if isfield(ws, 'status') && ~isempty(ws.status)
            f = f*ws.status;
        end 
        if isempty(f) || f == 0 && ( isfield(ws, 'sign') && ~isempty(ws.sign) )
            % If we have a sign for well and there is no flux defined, set
            % it to a small value of right sign instead.
            f = sqrt(eps)*ws.sign;
        end
        flux(i) = f;
    end
end
