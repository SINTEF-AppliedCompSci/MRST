function [gravmat, flux, q, pv, dt] = initTransport(G, state, ...
                                                    rock, fluid, varargin)
%Compute input to transport solver from flow calculation results.
%If specified, also determine timestep for use in explicit solver.
%
% SYNOPSIS:
%   [gravmat, ...
%    flux, q, pv, dt] = initTransport(G, state, rock, fluid)
%   [gravmat, ...
%    flux, q, pv, dt] = initTransport(G, state, rock, fluid, ...
%                                     'pn1', pv1, ...)
%
% PARAMETERS:
%   G       - Grid data structure discretising the reservoir model.
%
%   state   - Reservoir and well solution structure either properly
%             initialized from functions 'initResSol' and 'initWellSol'
%             respectively, or the results from a previous call to function
%             'solveIncompFlow' and, possibly, a transport solver such as
%             function 'implicitTransport'.
%
%   rock    - Rock data structure.  Must contain valid fields 'rock.perm',
%             with permeability measured in units of m^2, and 'rock.poro'.
%
%   fluid   - Data structure describing the fluids in the problem.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%
%               - wells  -- Well structure as defined by functions
%                           'addWell' and 'assembleWellSystem'.  May be
%                           empty (i.e., W = struct([])) which is
%                           interpreted as a model without any wells.
%
%               - bc  -- Boundary condtion structure as defined by function
%                        'addBC'.  This structure accounts for all external
%                        boundary contributions to the reservoir flow.
%                        Default value: bc = [] meaning all external
%                        no-flow (homogeneous Neumann) conditions.
%
%               - src -- Explicit source contributions as defined by
%                        function 'addSource'.
%                        Default value: src = [] meaning no explicit
%                        sources exist in the model.
%
%               - OnlyGrav  --  Whether or not to only compute gravity
%                               contribution (for use in gravity
%                               splitting).  Default value: false.
%
%               - ComputeDt --  Whether or not to compute timestep.
%                               Default value: true.
%
%               - max_dt    --  Maximal timestep (in units of seconds).
%                               Default value = 100 * 86400 (100 days).
%
%               - dt_factor --  Factor to multiply the computed dt by.
%                               Default value = 0.5.
%                               Ignored if ComputeDt = false.
%
%               - Verbose --
%                         Whether or not to emit progress reports during
%                         the process.
%                         Logical.  Default value = false.
%
% RETURNS:
%   flux    - In-flow connection matrix suitable for passing to transport
%             solver.  Content depends on presence or absence of gravity
%             effects:
%
%               - NORM(gravity()) == 0:
%                 flux contains only inflow,
%                 size(flux) = (G.cells.num, G.cells.num)
%               - NORM(gravity()) >  0:
%                 flux contains flow over all non-boundary faces,
%                 size(flux) = (G.cells.num, (G.faces.num-#bndFaces))
%
%             The fluxes are measured in units of m^3/s.
%
%   gravmat - Matrix with gravity contribution for each face.
%             Size G.cells.num-by-(G.faces.num-#bndFaces)).  The gravity
%             flux contributions, if present, are measured in units of
%             m^3/s.
%
%   q       - Aggregate source term suitable for passing to transport
%             solver.  The source term is measured in units of m^3/s.
%
%   pv      - Vector of size G.cells.num-by-1 of cell pore volumes.
%
%   dt      - Timestep (in units of seconds) for use in explicit solver.
%
% SEE ALSO:
%   `explicitTransport`, `implicitTransport`, `solveIncompFlow`, `gravity`.

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


   opt = struct('OnlyGrav' , false,       ...
                'ComputeDt', true,        ...
                'max_dt'   , 100 * day(), ...
                'dt_factor', 0.5,         ...
                'wells'    , [],          ...
                'src'      , [],          ...
                'bc'       , [],          ...
                'verbose'  , false);

   opt = merge_options(opt, varargin{:});
   assert (opt.dt_factor > 0 && opt.max_dt > 0);

   %-----------------------------------------------------------------------
   %% Assemble flux and gravity matrices ----------------------------------
   %
   is_int = all(double(G.faces.neighbors) > 0, 2);
   g_vec  = gravity();

   if norm(g_vec(1 : size(G.nodes.coords,2))) > 0,
      [gravmat, flux] = matrices(G, rock, fluid, state, ...
                                 opt.OnlyGrav, is_int);
   else
      [gravmat, flux] = matrices_nograv(G, state, is_int);
   end

   %-----------------------------------------------------------------------
   %% External sources/sinks (e.g. wells and BC's) ------------------------
   %
   q  = computeTransportSourceTerm(state, G, opt.wells, opt.src, opt.bc);
   pv = poreVolume(G, rock);

   %-----------------------------------------------------------------------
   %% Compute CFL-number/dt -----------------------------------------------
   %
   % CFL-condition:
   %    dt <= min(pv ./ (flux_in * max{f'(s)}_{0<=s<=1})
   %
   % We compute dt for flux_in = darcy flux and flux_in = flux due to
   % gravity, and use min(dt1, dt2).

   if opt.ComputeDt,
      dt = opt.max_dt;
      s  = linspace(0, 1, 1001) .';

      mu = fluid.properties(state);
      kr = fluid.relperm(s, state);

      deriv = @(f) diff(f(mu, kr)) ./ diff(s);

      if nnz(gravmat),
         f  = @(mu,kr) (prod(kr,2) ./ prod(mu,2)) ./ ...
                           sum(bsxfun(@rdivide, kr, mu), 2);

         dt_grav = estimate_dt(gravmat, deriv(f), pv, 0);
         if opt.verbose,
            fprintf('Estimated segregation step:\t%.2e [day]\n', ...
                    convertTo(dt_grav, day()));
         end
         dt = min(dt, dt_grav);
      end
      if ~opt.OnlyGrav
         f  = @(mu,kr) (kr(:,1) ./ mu(1)) ./ ...
                        sum(bsxfun(@rdivide, kr, mu), 2);

         dt_advection = estimate_dt(flux, deriv(f), pv, q);
         if opt.verbose,
            fprintf('Estimated advection step:\t%.2e [day]\n', ...
                    convertTo(dt_advection, day()));
         end
         dt = min(dt, dt_advection);
      end

      dt = dt * opt.dt_factor;
   else
      % Return dt = 0 if no estimate of timestep is requested.
      dt = 0;
   end
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function [gravity_matrix, flux] = matrices_nograv(G, state, is_int)
% matrices_nograv -- Compute flux and (trivial) gravity matrix in absence
%                    of gravity effects.
   nc = G.cells.num;

   % - Gravity matrix -
   gravity_matrix = sparse(G.cells.num, sum(double(is_int)));

   % - Flux matrix -
   % cflux is amount of outflow for all half contacts on all
   % cells.  Negative cell flux is flux into a cell.
   %
   cflux = faceFlux2cellFlux(G, state.flux);
   is_inflow    = cflux < 0;

   inflow_faces = G.cells.faces(is_inflow, 1);
   inflow_flux  = cflux(is_inflow);

   N            = double(G.faces.neighbors(inflow_faces, :));
   cellNo       = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
   sgn      = 2*double(G.faces.neighbors(G.cells.faces(:,1), 1)==cellNo)-1;

   ind          = sgn(is_inflow) > 0;

   N(ind,[1,2]) = N(ind,[2,1]);

   % In distributing the flow pattern, only consider internal faces when
   % constructing the connection matrix.  Flux boundary conditions, if any,
   % are considered sources and handled separately (function 'contrib_bc').
   %
   int_face = prod(N, 2) > 0;
   int_conn = N(int_face, :);
   int_flux = inflow_flux(int_face);

   flux = sparse(int_conn(:,2), int_conn(:,1), int_flux, nc, nc);
end

%--------------------------------------------------------------------------

function [gravity_matrix, flux] = matrices(G, rock, fluid, state, ...
                                           grav_only, is_int)
% matrices -- Compute flux and gravity matrix in presence of gravity.
% is_int - Logical vector of length G.faces.num, true for internal face.

   [nc, nf] = deal(G.cells.num, sum(double(is_int)));
   cellNo   = rldecode(1 : nc, diff(G.cells.facePos), 2) .';

   % Indices to (internal) half-faces.
   cIntFInx = is_int(G.cells.faces(:,1));

   % Global-to-local face map for internal faces.
   G2L      = cumsum(double(is_int));

   % Indices of internal face corresponding to cIntFInx.
   globfix  = G.cells.faces(cIntFInx, 1);

   % Renumbering globfix to 1:numel(globfix)
   locfix   = G2L(globfix);

   % - Gravity matrix -
   % Spread gravity to cellfaces: cellFInx * g_const .* cellFace_normal

   % rho(1) - rho(2) is always correct (even when rho(2) > rho(1)), because
   % during transport (using, e.g., the 'twophaseUpwBEGrav' transport
   % solver) we only modify the *first* saturation component.
   [rho, rho] = fluid.properties(state);
   g   = gravity() * (rho(1) - rho(2));
   dim = size(G.nodes.coords, 2);

   harm          = zeros([G.faces.num, 1]);
   renum         = zeros([G.faces.num, 1]);
   renum(is_int) = 1 : nf;

   % nKg == n' * K * g for all cellfaces.
   [K, r, c] = permTensor(rock, dim);

   assert (size(K,1) == G.cells.num, ...
          ['Permeability must be defined in active cells only.\n', ...
           'Got %d tensors, expected %d (== number of cells).'],   ...
           size(K,1), G.cells.num);

   nKg = sum(G.faces.normals(G.cells.faces(:,1), r) .* ...
             bsxfun(@times, K(cellNo,:), g(c)), 2);

   % Compute harmonic average of nKg on all *internal* faces.
   harm(is_int) = 2 ./ accumarray(renum(globfix), 1 ./ nKg(cIntFInx));

   sgn = 2*double(G.faces.neighbors(G.cells.faces(:,1), 1) == cellNo) - 1;

   % Compute final gravity matrix.
   gravity_matrix = sparse(cellNo(cIntFInx), locfix, ...
                           sgn(cIntFInx) .* harm(globfix), nc, nf);

   % - Flux matrix -
   if grav_only,
      flux = sparse(G.cells.num, nf);
   else
      cflux = faceFlux2cellFlux(G, state.flux);
      flux  = sparse(cellNo(cIntFInx), locfix, cflux(cIntFInx), nc, nf);
   end
end

%--------------------------------------------------------------------------

function dt = estimate_dt(flux_in, df, pv, q)
   % Consider max of inflow or outflow from a cell
   % when estimating time step size.

   if size(flux_in,2) == numel(q)
      flux_out = abs(sum(flux_in,1))' - min(q,0);
      flux_in  = abs(sum(flux_in,2))  + max(q,0);
   else
      flux_out = flux_in;
      flux_out(flux_out < 0) = 0;
      flux_out = abs(sum(flux_out, 2)) - min(q,0);

      flux_in(flux_in > 0) = 0;
      flux_in = abs(sum(flux_in, 2)) + max(q,0);
   end
   flux_max = max(flux_in, flux_out);
   is_flux  = abs(flux_max) > 0;

   if any(is_flux),
      dt = min((pv(is_flux) ./ flux_max(is_flux)) ./ max(abs(df)));
   else
      dt = inf;
   end
end
