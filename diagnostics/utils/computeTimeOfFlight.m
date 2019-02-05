function [T, A, q] = computeTimeOfFlight(state, G, rock, varargin)
%Compute time of flight using finite-volume scheme.
%
% SYNOPSIS:
%    T        = computeTimeOfFlight(state, G, rock)
%    T        = computeTimeOfFlight(state, G, rock, 'pn1', pv1, ...)
%   [T, A]    = computeTimeOfFlight(...)
%   [T, A, q] = computeTimeOfFlight(...)
%
% DESCRIPTION:
%   Compute time-of-flight by solving
%
%       v·\nabla T = \phi
%
%   using a first-order finite-volume method with upwind flux. Here, 'v' is
%   the Darcy velocity, '\phi' is the porosity, and T is time-of-flight.
%   The time-of-flight T(x) is the time it takes an inert particles that is
%   passively advected with the fluid to travel from the nearest inflow
%   point to the point 'x' inside the reservoir. For the computation to
%   make sense, inflow points must be specified in terms of inflow
%   boundaries, source terms, wells, or combinations of these.
%
% REQUIRED PARAMETERS:
%   G     - Grid structure.
%
%   rock  - Rock data structure.
%           Must contain a valid porosity field, 'rock.poro'.
%
%   state - Reservoir state structure, must contain valid cell interface
%           fluxes, 'state.flux'. Typically, 'state' will contain the
%           solution produced by a flow solver like 'incompTPFA'.
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
%   reverse - Boolean variable. If true, the reverse time-of-flight will be
%           computed, i.e., the travel time from a cell inside the
%           reservoir to the nearest outflow point (well, source, or
%           outflow boundary). Default: FALSE, i.e, compute forward tof.
%
%   tracer - Cell-array of cell-index vectors for which to solve a
%           stationary tracer equation:
%
%             \nabla·(v C) = 0,      C(inflow,t)=1
%
%           One equation is solved for each vector with tracer injected in
%           the cells with the given indices. Each vector adds one
%           additional right-hand side to the original tof-system. Output
%           given as additional columns in T.
%
%   computeWellTOFs - Boolean variable. If true, time-of-flight values are
%           computed individually for each influence region by solving
%           equations of the form:
%
%              \nabla·(v C_i T) = \phi C_i,    T(inflow)=0
%
%           where C_i denotes the tracer concentration associated with each
%           individual influence region
%
%   firstArrival - Boolean variable. If true, also compute the first-arrival
%           time by a graph algorithm (requires MATLAB R2015b or newer)

%   allowInf - Switch to turn off (true) or on (false) maximal TOF
%           thresholding. Default is false.
%
%   maxTOF - Maximal TOF thresholding to avoid singular/ill-conditioned
%           systems. Default (empty) is 50*PVI (pore-volumes-injected).
%           Only takes effect if 'allowInf' is set to false.
%
%   solver - Function handle to solver for use in TOF/tracer equations.
%           Default (empty) is matlab mldivide (i.e., \)
%
%   processCycles - Extend TOF thresholding to strongly connected
%           components in flux graph by considering Dulmage-Mendelsohn
%           decomposition (dmperm). Recommended for highly cyclic flux
%           fields. Only takes effect if 'allowInf' is set to false.
%
%
%
% RETURNS:
%   T - Cell values of a piecewise constant approximation to time-of-flight
%       computed as the solution of the boundary-value problem
%
%           (*)    v · \nabla T = \phi
%
%       using a finite-volume scheme with single-point upwind approximation
%       to the flux.
%
%   A - Discrete left-hand side of (*), a G.cells.num-by-G.cells.num matrix
%       whose entries are
%
%           A_ij = min(F_ij, 0), and
%           A_ii = sum_j max(F_ij, 0) + max(2*q_i, 0),
%
%       where F_ij = -F_ji is the flux from cell i to cell j
%
%           F_ij = A_ij n_ij · v_ij.
%
%       and n_ij is the outward-pointing normal of cell i for grid face ij.
%       The discretization uses a simple model for cells containing inflow.
%       If q_i denotes the rate and V_i the volume of cell i, then T_i is
%       set to half the time it takes to fill the cell, T_i = V_i/(2q_i).
%
%       OPTIONAL.  Only returned if specifically requested.
%
%   q - Aggregate source term contributions (per grid cell) from wells,
%       explicit sources and boundary conditions.  These are the
%       contributions referred to as 'q_i' in the definition of the matrix
%       elements, 'A_ii'.  Measured in units of m^3/s.
%
%       OPTIONAL.  Only returned if specifically requested.
%
% SEE ALSO:
%   `simpleTimeOfFlight`, `solveIncompFlow`.

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


opt = struct('bc',              [], ...
             'src',             [], ...
             'wells',           [], ...
             'reverse',         false,...
             'allowInf',        false, ...
             'maxTOF',          [], ...
             'tracer',          {{}}, ...
             'solver',          [], ...
             'processCycles',   false, ...
             'computeWellTOFs', false, ...
             'firstArrival',    false);

opt = merge_options(opt, varargin{:});

checkInput(G, rock, opt)
if size(state.flux, 2) > 1
    state = validateStateForDiagnostics(state);
end

tr = opt.tracer;
if ~iscell(tr), tr = {tr}; end

% Find external sources of inflow: q contains the contribution from from
% sources (src) and wells (W), while qb contains the cell-wise sum of
% fluxes from boundary conditions.
[q,qb] = computeSourceTerm(state, G, opt.wells, opt.src, opt.bc);

pv = poreVolume(G, rock);

if opt.reverse
   q  = -q;
   qb = -qb;
   state.flux = -state.flux;
end

% Build upwind flux matrix in which we define v_ji = max(flux_ij, 0) and
% v_ij = -min(flux_ij, 0). Then the diagonal of the discretization matrix
% is obtained by summing rows in the upwind flux matrix. This will give the
% correct diagonal in all cells except for those with a positive fluid
% source. In these cells, the average time-of-flight is set to half the
% time it takes to fill the cell, which means that the diagonal entry
% should be equal twice the fluid rate inside these cells.
i  = ~any(G.faces.neighbors==0, 2);
n  = double(G.faces.neighbors(i,:));
nc = G.cells.num;
qp = max(q+qb, 0);

out = min(state.flux(i), 0);
in  = max(state.flux(i), 0);

% Cell wise total inflow
inflow  = accumarray([n(:, 2); n(:, 1)], [in; -out]);

% The diagonal entries are equal to the sum of outfluxes minus divergence
% which equals the influx plus the positive source terms.
d = inflow + qp;

%--------------------------------------------------------------------------
% Handling of t -> inf (if opt.allowInf == false):
% We locate cells that have sufficiently small influx to reach maxTOF
% locally, i.e., t-t_up >= maxTOF <=> pv/influx >= maxTOF. System is
% modified such that these cells are set to maxTOF and upstream connections
% are removed.
% If no maxTOF is given, it is set to the time it takes to fill 50
% pore volumes (subject to change)
if ~opt.allowInf
    if isempty(opt.maxTOF)
        if ~opt.reverse, str = 'Forward '; else str = 'Backward'; end
        opt.maxTOF = 50*sum(pv)/sum(qp);
        dispif(mrstVerbose, ...
               '%s maximal TOF set to %5.2f years.\n', str, opt.maxTOF/year)
    end

    % Find cells that reach max TOF locally
    maxIn = max(d);
    aboveMax = full(pv./ (max(d, eps*maxIn)) ) > opt.maxTOF;

    % Set aboveMax-cells to maxTOF
    d(aboveMax)  = maxIn;
    pv(aboveMax) = maxIn*opt.maxTOF;
    % remove upstream connections
    in(aboveMax(n(:,2)))  = 0;
    out(aboveMax(n(:,1))) = 0;
end

% Inflow flux matrix
A  = sparse(n(:,2), n(:,1),  in, nc, nc)...
   + sparse(n(:,1), n(:,2), -out, nc, nc);
if ~isempty(d)
    A = -A + spdiags(d, 0, nc, nc);
end

if ~opt.allowInf && opt.processCycles
    [A, pv] = thresholdConnectedComponents(A, pv, maxIn, opt);
end

% Build right-hand sides for tracer equations. Since we have doubled the
% rate in any cells with a positive source, we need to also double the rate
% on the right-hand side here.
numTrRHS = numel(tr);
TrRHS = zeros(nc,numTrRHS);
for i=1:numTrRHS
   TrRHS(tr{i},i) = qp(tr{i});
end

% Time of flight for velocity field.
if isempty(opt.solver)
   T  = A \ [pv TrRHS];
else % if other solver, iterate over RHSs
   T = zeros(size(TrRHS)+[0, 1]);
   T(:, 1) = opt.solver(A, pv);
   for k = 2:size(T, 2)
       T(:, k) = opt.solver(A, TrRHS(:, k-1));
   end
end

% reset all tof > maxTOF to maxTOF
if ~opt.allowInf
    T(T>opt.maxTOF) = opt.maxTOF;
end

% compute individual well-tofs A*(c_i.*tau_i) = c_i.*pv
if opt.computeWellTOFs
    C   = T(:, 2:end);
    pvi = bsxfun(@times, C, pv);
    pvi(pvi<0) = 0;
    if isempty(opt.solver)
        X  = A \ pvi;
    else % if other solver, iterate over RHSs
        X = zeros(size(pvi));
        for k = 1:size(X, 2)
            X(:, k) = opt.solver(A, pvi(:, k));
        end
    end
    X(X<0) = 0;
    % disregard tracer values < sqrt(eps)
    ix     = and(pvi*opt.maxTOF > X, C > sqrt(eps));
    X(~ix) = opt.maxTOF;
    X(ix)  = X(ix)./C(ix);
    %X      = min(X, opt.maxTOF);
    T      = [T, X];

    if opt.firstArrival
      if verLessThan('matlab','8.6')
	    warning('Computing first arrival requires MATLAB R2015b or later');
        T  = [T, X];  % add dummy data equal TOF  
      else
	    n(out<0, [1 2]) = n(out<0, [2 1]);
        % add edges for well-to-connection cell
        nconn = cellfun(@numel, tr);
        wnum  = rldecode((1:numel(tr))', nconn(:));
        wc    = vertcat(tr{:});
        n1    = [n(:,1); wnum+nc];
        n2    = [n(:,2); wc];
        fl    = [inflow(n(:,2)); qp(wc)];
        fl    = max(fl, max(fl)*eps);
        gr    = digraph(n1, n2, pv(n2)./fl);
        X = T(:, 2:(numTrRHS+1));
        for k = 1:size(X,2)
            [~, t] = shortestpathtree(gr, nc+k, 'OutputForm', 'vector');
            X(:,k) = t(1:nc);
        end
        X = min(X, opt.maxTOF);
        T = [T, X];
      end
    end
end

end


function [q, qb] = computeSourceTerm(state, G, W, src, bc)
   qi = [];  % Cells to which sources are connected
   qs = [];  % Actual strength of source term (in m^3/s).

   % Contribution from wells
   if ~isempty(W)
      qi = [qi; vertcat(W.cells)];
      qs = [qs; vertcat(state.wellSol.flux)];
   end

   % Contribution from sources
   if ~isempty(src)
      qi = [qi; src.cell];
      qs = [qs; src.rate];
   end

   % Assemble all source and sink contributions to each affected cell.
   q = sparse(qi, 1, qs, G.cells.num, 1);

   % Contribution from boundary conditions
   if ~isempty(bc)
      ff    = zeros(G.faces.num, 1);

      isDir = strcmp('pressure', bc.type);
      i     = bc.face(isDir);
      if ~isempty(i)
         ff(i) = state.flux(i) .* (2*(G.faces.neighbors(i,1)==0) - 1);
      end

      isNeu = strcmp('flux', bc.type);
      ff(bc.face(isNeu)) = bc.value(isNeu);

      is_outer = ~all(double(G.faces.neighbors) > 0, 2);
      qb = sparse(sum(G.faces.neighbors(is_outer,:), 2), 1, ...
         ff(is_outer), G.cells.num, 1);
   else
      qb = sparse(G.cells.num,1);
   end
end


function [A, pv] = thresholdConnectedComponents(A, pv, maxIn, opt)
    % Find strongly connected components in flux-matrix:
    [p,r,r]= dmperm(A); %#ok
    % Pick components containing more than a single cell
    ix = find(diff(r)>1);
    if ~isempty(ix)
        nc = numel(p);
        % Retrieve cell-indices to components, and construct sparse index-mapping
        c  = arrayfun(@(b,e)p(b:e)', r(ix), r(ix+1)-1, 'UniformOutput', false);
        rc = rldecode( (1:numel(c))', cellfun(@numel, c)');
        C  = sparse(vertcat(c{:}), rc, 1, nc, numel(c));
        % Compute influx to each component
        q_in = full(diag(C'*A*C));
        % Threshold
        compAboveMax = full((C'*pv)./ (max(q_in, eps*maxIn)) ) > opt.maxTOF;

        if any(compAboveMax)
            badCells = (vertcat(c{compAboveMax}));
            dispif(mrstVerbose, 'Found %d strongly connected components, ', nnz(compAboveMax));
            dispif(mrstVerbose, 'total of %d cells, with influx below threshold.\n', numel(badCells));
            % Modify system setting tof to maxTOF and remove upstream connections
            A(badCells,:) = sparse((1:numel(badCells))', badCells, maxIn, numel(badCells), nc);
            pv(badCells)  = maxIn*opt.maxTOF;
        end
    end
end


function checkInput(G, rock, opt)
    assert (~all([isempty(opt.src), isempty(opt.bc), isempty(opt.wells)]), ...
        'Must have inflow described as boundary conditions, sources, or wells');
    assert (isfield(rock, 'poro')         && ...
        numel(rock.poro)==G.cells.num,   ...
        ['The rock input must have a field poro with porosity ',...
        'for each cell in the grid.']);
    assert(min(rock.poro) > 0, 'Rock porosities must be positive numbers.');
    if ~isempty(opt.maxTOF) && opt.allowInf
        warning('Input value for ''maxTOF'' ignored since option ''allowInf'' has value ''true''.')
    end
    if opt.processCycles && opt.allowInf
        warning('Input request to process cycles ignored since option ''allowInf'' has value ''true''.')
    end
end

