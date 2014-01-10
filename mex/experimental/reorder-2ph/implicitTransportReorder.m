function [state, report] = implicitTransportReorder(state, G, tf, rock, fluid, varargin)
% implicitTransportReorder -- Single point upwind solver for Buckley-Leverett flow.
%
% SYNOPSIS:
%   state = implicitTransportReorder(state, G, tf, rock, fluid)
%   state = implicitTransportReorder(state, G, tf, rock, fluid, 'pn1', pv1,...)
%
% DESCRIPTION:
%   Function implicitTransportReorder solves the Buckley--Leverett transport
%   equation
%
%        s_t + · \nabla [f(s)v·n] = f(s)q
%
%   using a first-order mobility weighted upwind discretisation in space
%   and backward Euler discretisation in time. The transport equation is
%   solved on the time interval [0, tf] calling the C-solver implicitupwind
%   via a mex interface. This solver uses a topological sort of the flux
%   graph to permute the nonlinear system into a block-triangular form
%   leading to a Gauss-Seidel type nonlinear solver in which the block
%   systems are computed in sequence, moving gradually downstream.
%
% REQUIRED PARAMETERS:
%   state - Reservoir and well solution structure either properly
%           initialized from function 'initState', or the results from a
%           previous call to function 'solveIncompFlow' and, possibly, a
%           transport solver such as function 'explicitTransport'.
%
%   G     - Grid data structure discretising the reservoir model.
%
%   tf    - End point of time integration interval (i.e., final time).
%           Measured in units of seconds.
%
%   rock  - Rock data structure.  Must contain the field 'rock.poro' and,
%           in the presence of gravity or capillary forces, valid
%           permeabilities measured in units of m^2 in field 'rock.perm'.
%
%   fluid - Fluid data structure as defined in 'fluid_structure'. However,
%           the only fluid that is yet supported over the mex interface is
%           the one defined in initCoreyFluidROC
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%
%   wells - Well structure as defined by function 'addWell'.  This
%           structure accounts for all injection and production well
%           contribution to the reservoir flow. Default value: wells = [],
%           meaning a model without any wells.
%
%   bc    - Boundary condtion structure as defined by function
%          'addBC'.  This structure accounts for all external boundary
%           contributions to the reservoir flow. Default value: bc = []
%           meaning all external no-flow (homogeneous Neumann) conditions.
%
%   src   - Explicit source contributions as defined by function
%           'addSource'. Default value: src = [] meaning no explicit
%           sources exist in the model.
%
% RETURNS:
%   state - Reservoir solution with updated saturation, state.s.
%
% EXAMPLE:
%   See <....>
%
% SEE ALSO:
%   implicitTransport, explicitTransport, initCoreyFluidROC

opt=struct('bc', [], 'src', [], 'wells', [],'substeps',1,...
    'max_iterations',100,'verbose',mrstVerbose());
opt = merge_options(opt, varargin{:});

assert (isempty(opt.src) || ~isempty(opt.src.sat), ...
   'Source terms must have a valid ''sat'' field.');
assert (isempty(opt.bc) || ~isempty(opt.bc.sat), ...
   'Boundary conditions must have a valid ''sat'' field.');
assert (isfield(rock, 'poro')         && ...
    numel(rock.poro)==G.cells.num,   ...
    ['The rock input must have a field poro with porosity ',...
      'for each cell in the grid.']);
assert(min(rock.poro) > 0, 'Rock porosities must be positive numbers.');

q  = computeTransportSourceTerm(state, G, opt.wells, opt.src, opt.bc);

% Make flux matrix
i  = ~any(G.faces.neighbors==0, 2);
n  = double(G.faces.neighbors(i,:));
nc = G.cells.num;
A  = sparse(n(:,1), n(:,2),  max(state.flux(i), 0), nc, nc)...
   + sparse(n(:,2), n(:,1), -min(state.flux(i), 0), nc, nc); clear n i;
A  = -A' + spdiags(sum(A)' + max(q, 0), 0, nc, nc);

fopt = fluid.param();
fparam = [reshape(fopt.n,2,[])', reshape(fopt.sr,2,[])', reshape(fopt.kwm,2,[])']';
fparam = [fopt.mu fparam(:)'];

n_verb=0; if(opt.verbose), n_verb=2; end

opt_re = struct(...
   'max_iterations', int32(opt.max_iterations), ...  % nonlinear iterations
   'substeps',       int32(opt.substeps),       ...
   'tolerance',      1e-6,                      ...  % nonlinear tolerance
   'verbosity',      int32(n_verb),             ...  % level of reporting
   'scalar_solver',  'ridder');                      % nonlinear solver

opt_re.verbosity      = int32(opt_re.verbosity);
opt_re.max_iterations = int32(opt_re.max_iterations);
opt_re.substeps       = int32(opt_re.substeps);

if nargout == 1
    state.s = implicitupwind(-A', poreVolume(G,rock), state.s(:,1), ...
        tf, full(q), fparam, int32(fopt.reg)-1, opt_re);
elseif nargout == 2
    [state.s, report] = implicitupwind(-A', poreVolume(G,rock), ...
        state.s(:,1), tf, full(q), fparam, int32(fopt.reg)-1, opt_re);
else
    error('Too many output arguments');
end
end