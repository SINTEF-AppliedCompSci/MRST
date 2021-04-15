function state = explicitTransport(state, G, tf, rock, fluid, varargin)
%Explicit single-point upstream mobility-weighted transport solver for two-phase flow.
%
% SYNOPSIS:
%   state = explicitTransport(state, G, tf, rock, fluid)
%   state = explicitTransport(state, G, tf, rock, fluid, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Implicit discretization of the transport equation
%
%      s_t + div[f(s)(v + mo K((rho_w-rho_o)g + grad(P_c)))] = f(s)q
%
%   where v is the sum of the phase Darcy fluxes, f is the fractional
%   flow function,
%
%                  mw(s)
%        f(s) = -------------
%               mw(s) + mo(s)
%
%   mi = kr_i/mu_i is the phase mobility of phase i, mu_i and rho_i are the
%   phase viscosity and density, respectively, g the (vector) acceleration
%   of gravity, K the permeability, and P_c(s) the capillary pressure.  The
%   source term f(s)q is a volumetric rate of water.
%
%   We use a first-order upstream mobility-weighted discretization in space
%   and a backward Euler discretization in time. The transport equation is
%   solved on the time interval [0,tf] by calling twophaseJacobian to build
%   a function computing the residual of the discrete system in addition to
%   a function taking care of the update of the solution during the
%   time loop.
%
% REQUIRED PARAMETERS:
%   state - Reservoir and well solution structure either properly
%           initialized from function 'initState', or the results from a
%           previous call to function 'solveIncompFlow' and, possibly, a
%           transport solver such as function 'explicitTransport'.
%
%   G     - Grid data structure discretizing the reservoir model.
%
%   tf    - End point of time integration interval (i.e., final time).
%           Measured in units of seconds.
%
%   rock  - Rock data structure.  Must contain the field 'rock.poro' and,
%           in the presence of gravity or capillary forces, valid
%           permeabilities measured in units of m^2 in field 'rock.perm'.
%
%   fluid - Fluid data structure as defined in 'fluid_structure'.
%
% OPTIONAL PARAMETERS:
%   W         - Well structure as defined by function 'addWell'.  This
%               structure accounts for all injection and production well
%               contribution to the reservoir flow.
%               Default value: W = [], meaning a model without any wells.
%
%   bc        - Boundary condition structure as defined by function
%               'addBC'.  This structure accounts for all external boundary
%               contributions to the reservoir flow.
%               Default value: bc = [] meaning all external no-flow
%               (homogeneous Neumann) conditions.
%
%   src       - Explicit source contributions as defined by function
%               'addSource'. Default value: src = [] meaning no explicit
%               sources exist in the model.
%
%   onlygrav  - Ignore content of state.flux.         Default false.
%
%   computedt - Estimate time step.                   Default true.
%
%   max_dt    - If 'computedt', limit time step.      Default inf.
%
%   dt_factor - Safety factor in time step estimate.  Default 0.5.
%
%   dt        - Set time step manually.  Overrides all other options.
%
%   gravity   - The current gravity in vector form. Defaults to gravity().
%
%   satwarn   - Issue a warning if saturation is more than 'satwarn'
%               outside the default interval of [0,1].
%
% RETURNS:
%   state     - Reservoir solution with updated saturation, state.s.
%
% EXAMPLE:
%   See simple2phWellExample.m
%
% SEE ALSO:
%   `twophaseJacobian`, `implicitTransport`.

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


   opt = struct(...
      'verbose'  , false      , ...
      'onlygrav' , false      , ...
      'computedt', true       , ...
      'max_dt'   , inf        , ...
      'dt_factor', 0.5        , ...
      'wells'    , []         , ...
      'W'        , []         , ...
      'src'      , []         , ...
      'bc'       , []         , ...
      'dt'       , tf         , ...
      'Trans'    , []         , ...
      'gravity'  , gravity()  , ...
      'satwarn'  , sqrt(eps));

   opt = merge_options(opt, varargin{:});
   opt = treatLegacyForceOptions(opt);


   if opt.onlygrav
      flux = state.flux;
      state.flux = zeros(G.faces.num, 1);
   end

   [F,~,gf,q] = twophaseJacobian(G, state, rock, fluid, ...
                        'wells', opt.wells,    ...
                        'src'  , opt.src  ,    ...
                        'bc'   , opt.bc,       ...
                        'Trans', opt.Trans);


   if opt.computedt

      % ---------- Time step estimate from state ---------------
      compi = { 'use_compi', true };
      vsrc = computeTransportSourceTerm(state, G, opt.wells, ...
                                        opt.src, opt.bc, compi{:});
      vsrc = assembleTransportSource(state, fluid, vsrc, G.cells.num, compi{:});

      gflux = getFlux(G, rock,opt);
      getdt = @(state) min([...
                   opt.max_dt, ...
                   opt.dt_factor * estimate_dt(G, state, rock, fluid, ...
                                                state.flux, gflux, vsrc)]);
   else

      % ----------- Constant time step -------------------------
      getdt =@(state) opt.dt;

   end


   s  = state.s(:,1);
   t  = 0;
   dispif(opt.verbose, 'explicitTransport: Computing transport step in %d substeps\n', ...
         ceil(tf/getdt(state)));
   while t < tf
      dt      = min(tf-t, getdt(state));

      s(:)    = s - F(state, state, dt);
      t       = t + dt;

      s       = correct_saturations(s, opt.satwarn);

      state.s = [s, 1-s];

      if isfield(state, 'extSat')
         % Save minimum saturation for use in modeling of relative
         % permeability hysteresis.
         state.extSat(:,1) = min(state.s(:,1), state.extSat(:,1));
         state.extSat(:,2) = max(state.s(:,1), state.extSat(:,2));
      end
   end

   if opt.onlygrav
      state.flux = flux;
   end

   if any(any(isnan(state.s)))
      error('explicitTransport: Transport step failed')
   end
end

%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------

function gflux = getFlux(G, rock, opt)
%harmonic average of one-sided n \cdot K \cdot g on each face

   gvec   = -opt.gravity;
   gflux  = zeros([G.faces.num, 1]);

   dim    = size(G.nodes.coords, 2);
   if(isempty(opt.Trans))
      [K, r, c] = permTensor(rock, dim);

      assert (size(K,1) == G.cells.num, ...
             ['Permeability must be defined in active cells only.\n', ...
              'Got %d tensors, expected %d (== number of cells).'],   ...
              size(K,1), G.cells.num);
   end
   nc     = G.cells.num;
   cellNo = rldecode(1 : nc, diff(G.cells.facePos), 2) .';

   if norm(gvec) > 0

      % nKg == n' * K * g for all cellfaces.
      nKg    = sum(G.faces.normals(G.cells.faces(:,1), r) .* ...
                   bsxfun(@times, K(cellNo,:), gvec(c)), 2);

      % Compute harmonic average of one-sided nKg on all faces.
      gflux = 2 ./ accumarray(G.cells.faces(:,1), 1 ./ nKg, [G.faces.num, 1]);
   end


end


%--------------------------------------------------------------------------

function dt = estimate_dt(G, state, rock, fluid, flux, gflux, sources)
   [rho, kr, mu, dkr] = getIncompProps(state, fluid);

   % Compute cell mobility and its derivative
   sat       = state.s;
   mob       = bsxfun(@rdivide, kr, mu);

   % dkr is Jacobian of kr.  We need derivatives with respect to s(:,1),
   % hence sign of 'dkr(:,end)'.
   dmob = bsxfun(@rdivide, [dkr(:,1), -dkr(:,end)], mu);

   % Compute face density as average of cell densities.
   i = all(G.faces.neighbors > 0, 2);
   N = G.faces.neighbors(i, :);

   % Find simple approximation to the maximal wave speed from advective
   % term in reservoir based on derivative of flux on face.
   df = @(mob, dmob) ...
       (mob(:,2).*dmob(:,1) - mob(:,1).*dmob(:,2))./(sum(mob, 2).^2);
   d  = df(mob, dmob);
   m  = max(abs(d(N)), [], 2);

   wavespeed  = max(abs(m.*flux(i)));

   % Find max wave speed from segregation term
   g      = @(mob) mob(:,1).*mob(:,2)./sum(mob, 2);
   s      = linspace(0,1, 101)';
   ss     = struct('s',[s,1-s]);
   [rho_loc, kr_loc, mu_loc] = getIncompProps(ss, fluid); %#ok<ASGLU>
   m      = bsxfun(@rdivide, kr_loc, min(mu_loc, [], 2));
   dg     = max(abs(diff(g(m)))./diff(s));

   wavespeed  = max([wavespeed; abs(dg.*gflux(i).*max(diff(rho_loc, 1, 2), [], 1))]  );

   % Find max wave speed from advective term for positive sources in
   % interval [s, 1],
   i      = sources > 0;
   if any(i)
      s     = (min(sat(i)) : 0.05 : 1.0)';
      ss    = struct('s',[s,1-s]);
      if numel(s) > 1
         [tmp, kr, mu_w] = getIncompProps(ss, fluid); %#ok<ASGLU>
         mob   = bsxfun(@rdivide,  kr, mu_w);
         f     = bsxfun(@rdivide, mob(:,1), sum(mob, 2));
         maxdf = max(diff(f)./diff(s));
         wavespeed  = max(wavespeed, max(abs(sources(i) .* maxdf)));
      end
   end

   % Find max wave speed from advective term for negative sources in
   % interval [0, s],
   i      = sources < 0;
   if any(i)
      s     = (max(sat(i)) : -0.05 : 0.0)';
      ss=struct('s',[s,1-s]);
      if numel(s) > 1
         [kr, mu] = getIncompProps(ss, fluid);
         mob   = bsxfun(@rdivide,  kr, mu);
         f     = bsxfun(@rdivide, mob(:,1), sum(mob, 2));
         maxdf = max(diff(f)./diff(s));
         wavespeed  = max(wavespeed, max(abs(sources(i) .* maxdf)));
      end
   end

   dt = min(abs(poreVolume(G, rock)/wavespeed));
end

function s = correct_saturations(s, satwarn)
   % Display error if s > 1+satwarn
   % or s < 0 - satwarn
   i = find(s(:,1) > 1 + satwarn);
   if ~isempty(i)
      disp('Saturation exceeds 1 in cells:')
      fprintf('%5d %g\n', [i, s(i,1)] .');
      error('explicitTransport: Error larger than satwarn')
   end

   i = find(s(:,1) < -satwarn);
   if ~isempty(i)
      disp('Saturation less than 0 in cells:')
      fprintf('%5d %g\n', [i, s(i,1)] .');
      error('explicitTransport: Error larger than satwarn')
   end
   % Correct numerical errors
   s(s(:,1) > 1, 1) = 1 - eps;
   s(s(:,1) < 0, 1) = 0;
end
