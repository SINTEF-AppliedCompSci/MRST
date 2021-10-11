function [Res, Jac, gflux, q] = twophaseJacobian(G, state, rock, fluid, varargin)
%Residual and Jacobian of single point upwind solver for two-phase flow.
%
% SYNOPSIS:
%   F       = twophaseJacobian(G, state, rock, fluid)
%   [F,Jac] = twophaseJacobian(G, state, rock, fluid, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Function twophaseJacobian returns function handles for the residual
%   and its Jacobian matrix for the implicit upwind-mobility weighted
%   discretization of
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
%   Using a first-order upstream mobility-weighted discretization in space
%   and a backward Euler discretization in time, the residual of the
%   nonlinear system of equations that must be solved to move the solution
%   state.s from time=0 to time=tf, is obtained by calling F(s,s0,dt),
%   defined as
%      F(s, s0, dt) = s - s0 + dt/pv·[H(s) - max(q,0) - min(q,0)·f],
%
%      H(s) = sum_i [f_i(dflux_i + mo_i*(gflux_i)+pcflux)]  (faces i of cell)
%
%   where f_i = mw_i/(mw_i+mo_i), mo_i and mw_i are phase upwind
%   mobilities at face i and dflux, gflux, and pcflux are Darcy flux,
%   gravity flux equal face_normal*Kg(rho1-rho2), and capillary flux,
%   respectively. The other flux function f is the fractional flow function
%   evaluated with cell mobilities. Likewise, the Jacobian matrix is
%   obtained using the function Jac.
%
% REQUIRED PARAMETERS:
%   resSol  - Reservoir solution structure containing valid
%             saturation resSol.s with one value for each cell in
%             the grid.
%
%   G       - Grid data structure discretizing the reservoir model.
%
%   rock    - Struct with fields perm and poro.  The permeability field
%             ('perm') is only referenced when solving problems involving
%             effects of gravity and/or capillary pressure and need not be
%             specified otherwise.
%
%   trans    - two point flux transmissibilities. If given this will be
%              used not rock
%
%   fluid   - Data structure describing the fluids in the problem.
%
% OPTIONAL PARAMETERS:
%
%   gravity   - The current gravity in vector form. Defaults to gravity().
%
%
% RETURNS:
%   F       - Residual
%
%   Jac     - Jacobian matrix (with respect to s) of residual.
%
% EXAMPLE:
%
%
%
% SEE ALSO:
%   `implicitTransport`, `explicitTransport`.

% TODO:
%   - implement gravity effects for pressure boundary and wells
%   - multipliers for gravity flux.  Handle this in 'getFlux'

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

% NOTE we solve
%              __
%        s_t + \/· [D·f(s) + (G + P)·f(s)·mo(s)] = f(s) Q
%
% for constant vector and scalar fields D (Darcy flux), G (gravity flux),
% P (capillary flux), and Q (source term).

   opt = struct('verbose', mrstVerbose, 'gravity', gravity(), ...
                'wells', [], 'src', [], 'bc', [], 'Trans', [],'dhfz',[]);
   opt = merge_options(opt, varargin{:});

   assert ((size(state.s,2) < 3) || all(state.s(:,3) == 0), ...
           'Function ''%s'' is for two-phase flow only', mfilename);

   assert(all(isfinite(state.flux)), 'Passed state contained non-finite fluxes');
   % All source terms, i.e., boundary conditions, source terms and wells.
   compi = { 'use_compi', true };
   q = computeTransportSourceTerm(state, G, opt.wells, ...
                                  opt.src, opt.bc, compi{:});
   q = assembleTransportSource(state, fluid, q, G.cells.num, compi{:});
   assert(all(isfinite(q)))
   % Extract (constant) fluid densities.
   rho = getIncompProps(state, fluid);

   % Compute the gravitational potential (rho_w - rho_o)n·Kg across each
   % face in the grid.
   [cellNo, cellFaces] = getCellNoFaces(G);
   [gflux, pc_flux, pcJac] = getFlux(G, cellNo, cellFaces, rock, rho, fluid, opt);

   % Bind flux terms state.flux and gflux.  The call
   % [iw,io]=findUpwindCells(..) computes cell index of upwind cell for the
   % water phase (iw) and oil phase (io).
   findUpwindCells = @(mob, state) upwindIndices(G, state.flux, gflux, ...
                                                 pc_flux(state), mob);

   pv = poreVolume(G, rock);

   % Return function handles
   Res = @(xr, xr_0, dt) Residual(xr, xr_0, dt, fluid, state.flux, gflux, ...
                                pc_flux, findUpwindCells, G, pv, q);

   Jac = @(xr, xr_0, dt)     Jacobian(xr, xr_0,    dt, fluid, state.flux, gflux, ...
                                pc_flux, pcJac, findUpwindCells, G, pv, q);
end

%--------------------------------------------------------------------------

function J = Jacobian (resSol, resSol_0, dt, fluid, ...
                       dflux, gflux, pcflux, pcJac, ...
                       findPhaseUpwindCells, G, pv, q)
%  At an interface, we compute the Jacobain of F as follows,
%
%    dF      dt      mo       1
%    --- =   -- · -------  -----  (dflux + mo*gflux)
%    dmw     pv   mw + mo  mw + mo
%
%            dt  fo
%        =   --  --  (dflux + mo*gflux)
%            pv  mt
%
%
%    dF      dt  fw
%    --- = - --  --  (dflux - mw*gflux)
%    dmo     pv  mt
%
%    dF          dt  fw
%    ------- =   --  --  mo
%    dpcflux     pv
%
%
%    dF    dF  dmw    dF  dmo   dF      dpcflux
%    --  = --- ---  + --- --- + ------- -------
%    dS    dmw dS     dmo dS    dpcflux   dS
%
%  If the mobilities mw and mo are evaluated in grid cells iw and io,
%  respectively, then the respective (matrix) contributions of each
%  term is added to elements (i, iw) and (i, io) of the Jacobian matrix.

   if isfield(resSol, 'extSat')
      % Save extremal values of saturation for hysteresis purposes
      % First column contains minima, second maxima
      resSol.extSat(:,1) = min(resSol.s(:,1), resSol_0.extSat(:,1));
      resSol.extSat(:,2) = max(resSol.s(:,1), resSol_0.extSat(:,2));
   end
   [tmp, kr, mu, dkr]  = getIncompProps(resSol, fluid); %#ok<ASGLU>
   m = bsxfun(@rdivide, kr,  mu);    % mobility in each cell

   % We need derivatives with respect to s(:,1) only.
   dkr(:,end) = -dkr(:,end);

   % Mobility derivative.  Two-phase only.  Extract diagonal of Jacobian.
   dm = bsxfun(@rdivide, dkr(:, [1, end]), mu);

   clear mu sat kr dkr

   neighbors = getNeighbourship(G, 'Topological', true);
   internal   = all(neighbors~=0, 2);
   ic1  = double(neighbors(internal,1));
   ic2  = double(neighbors(internal,2));
   clear internal_faces


   % For each face k, iw(k) and io(k) is the cell index in which we
   % must evaluate the water and oil mobilities, respectively.
   [iw, io] = findPhaseUpwindCells(m, resSol);
   if ~any(gflux) && ~isfield(fluid, 'pc')
      assert(all(iw == io));
   end
   m_face   = [m(iw, 1) m(io,2)];
   f_face   = bsxfun(@rdivide, m_face, sum(m_face,2));

   % Compute Jacobi contributions from interfaces
   dmw_fo    = f_face(:,2).*dm(iw,1);
   dmo_fw    = f_face(:,1).*dm(io,2);

   % contribution from capillary forces (pcflux) must be evaluated from
   % saturation at each step
   dvw_ds1 =  (dflux(internal) + (pcflux(resSol) + gflux(internal)) ...
              .*m_face(:,2)) .* dmw_fo ./ sum(m_face,2);
   dvw_ds2 = -(dflux(internal) - (pcflux(resSol) + gflux(internal)) ...
              .*m_face(:,1)) .* dmo_fw ./ sum(m_face,2);

   dFdSiw1  = dt./pv(ic1).*dvw_ds1;
   dFdSiw2  = dt./pv(ic2).*dvw_ds1;
   dFdSio1  = dt./pv(ic1).*dvw_ds2;
   dFdSio2  = dt./pv(ic2).*dvw_ds2;


   clear dmw_fo dmo_fw dvw_ds1 dvw_ds2


   % In cell: differentiate term s-s0 + dt/pv(max(q,0) + min(q,0).*f) with
   % respect to cell saturation.
   f   = bsxfun(@rdivide, m, sum(m,2));
   df  = (f(:,2).*dm(:,1) - f(:,1).*dm(:,2))./sum(m,2);
   d   = 1-dt./pv.*(min(q, 0).*df);


   id  = (1:G.cells.num)';

   % Assemble contribututions.
   %           diagonal  water mobility      oil mobility

   J = sparse(...
      [id;        ic1;      ic2;        ic1;      ic2], ... %row
      [id;         iw;       iw;         io;       io], ... %column
      [d;     dFdSiw1; -dFdSiw2;    dFdSio1; -dFdSio2], ... %value
      G.cells.num, G.cells.num);
  % Get capillary pressure
  [pc, dpc] = getIncompCapillaryPressure(resSol, fluid);
  if ~isempty(pc)
      d_p  =      dt.* f_face(:,1).*m_face(:,2).*pcJac(internal);

      % d pcflux
      %  --------
      %     dS
      % for a face e_ij gives four contributions to the Jacobi matrix:

      dFdSipc1_1 =  d_p.*dpc(ic1)./pv(ic1);  % J_ii
      dFdSipc1_2 =  d_p.*dpc(ic2)./pv(ic1);  % J_ij
      dFdSipc2_1 =  d_p.*dpc(ic1)./pv(ic2);  % J_ji
      dFdSipc2_2 =  d_p.*dpc(ic2)./pv(ic2);  % J_jj


      J = J + sparse(...
         [ic1;                ic1;          ic2;        ic2], ... %row
         [ic1;                ic2;          ic1;        ic2], ... %column
         [dFdSipc1_1; -dFdSipc1_2;  -dFdSipc2_1; dFdSipc2_2], ... %value
         G.cells.num, G.cells.num);
   end

end

%--------------------------------------------------------------------------

function F = Residual (resSol, resSol_0, dt, fluid, ...
                       dflux, gflux, pcflux,        ...
                       findPhaseUpwindCells, G, pv, q)
%  DISCRETISATION
%
%  Compute  F given by a residual function Fk for each grid cell
%
%   Fk(s, s0, dt) = s - s0 + dt/pv·[H(s) - max(q,0) - min(q,0)·f],
%
%   H(s) = sum_i [f_i(dflux_i + mo_i*(gflux_i)+pcflux)]  (faces i of cell)
%
%  where f_i = mw_i/(mw_i+mo_i), mo_i and mw_i are phase upwind
%  mobilities at face i and dflux, gflux, and pcflux are Darcy flux,
%  gravity flux (=face_normal*Kg(rho1-rho2)), and capillary flux,
%  respectively. The other flux function f is the fractional flow function
%  evaluated with cell mobilities.
%
%  Advancing an implicit upwind mobility weighted scheme one time step
%  amounts to solving F(s, s0, dt) = 0, given s0 and dt.
%
%  The upwind mobility weighted flux f_i is simply computed by evaluating
%  the water mobility mw in the upwind cell wrt water flux
%
%    wflux = fw(dflux + mo*gflux).
%
%  Similarily, we evaluate mo in the upwind cell wrt oil flux
%
%    oflux = fo(dflux - mw*gflux).
%
%  The upwind cell for each phase at each internal grid face is found
%  using findPhaseUpwindCells, and computeConstData.  The gravity flux and
%  capillary flux are computed in getFlux.

   % Compute cell mobilities
   s0 = resSol_0.s;
   if isfield(resSol, 'extSat'),
      % Save extremal values of saturation for hysteresis purposes
      % First column contains minima, second maxima
      resSol.extSat(:,1) = min(resSol.s(:,1), resSol_0.extSat(:,1));
      resSol.extSat(:,2) = max(resSol.s(:,1), resSol_0.extSat(:,2));
   end
   [tmp, kr, mu]  = getIncompProps(resSol, fluid); %#ok<ASGLU>
   
%    mu  = fluid.properties(resSol);
%    sat = fluid.saturation(resSol);
%    kr  = fluid.relperm(sat, resSol);
   m   = bsxfun(@rdivide, kr, mu);
   clear mu sat kr tmp

   neighbors = getNeighbourship(G, 'Topological', true);
   internal = all(neighbors~=0, 2);
   ic1  = neighbors(internal,1);
   ic2  = neighbors(internal,2);
   clear internal_faces

   % Compute water source term
   f_cell = m(:,1) ./ sum(m,2);
   Q      = max(q,0) + min(q,0).*f_cell;
   clear f_cell

   % Compute the upwind mobility weighted flow of water:
   % For each face k, iw(k) and io(k) is the cell index in which we
   % must evaluate the water and oil mobilities, respectively.
   [iw, io] = findPhaseUpwindCells(m, resSol);
   if ~any(gflux) && ~isfield(fluid, 'pc')
      assert(all(iw == io));
   end
   m_face  = [m(iw, 1) m(io,2)];
   f_face  = m_face(:,1)./sum(m_face,2);

   v_water = f_face.*(dflux(internal)+(gflux(internal) + pcflux(resSol)).*(m_face(:,2)));
   clear io iw m m_face f_face

   F   = resSol.s(:,1) - s0(:,1);
   F   = F + dt.*(1./pv).*(accumarray([ic1; ic2], [v_water; -v_water], ...
                                    size(F)) - Q);
end

%--------------------------------------------------------------------------

function [gflux, pc_flux, pcJac]  = getFlux(G, cellNo, cellFace, rock, rho, fluid, opt)
% -------------------------------------------------------------------------
%
%                        TECHNICAL DESCRIPTION
%
% -------------------------------------------------------------------------
%
%  gflux = harmonic average of (n·K·gravity·(rho1-rho2)) in each cell-

   g     = opt.gravity * (rho(1) - rho(2));
   [N, neighborship] = deal(getNeighbourship(G, 'Topological', true));

   % Number of interfaces (faces + NNC)
   nif = size(N, 1);

   gflux = zeros([size(N, 1), 1]);
   dim   = G.griddim;

   if (norm(g) > 0) || isfield(fluid, 'pc')
      if isempty(opt.Trans)
         [K, r, c] = permTensor(rock, dim);

         assert (size(K,1) == G.cells.num, ...
                ['Permeability must be defined in active cells only.\n', ...
                 'Got %d tensors, expected %d (== number of cells).'],   ...
                 size(K,1), G.cells.num);
      end
   end

   if norm(g) > 0
      if isempty(opt.Trans)
        % nKg == n' * K * g *(rho1 - rho2) for all cellfaces.
        nKg    = sum(G.faces.normals(G.cells.faces(:,1), r) .* ...
                   bsxfun(@times, K(cellNo,:), g(c)), 2);
        % Compute harmonic average of nKg(rho1-rho2) on all faces.
        gflux = 2 ./ accumarray(G.cells.faces(:,1), 1 ./ nKg, [G.faces.num, 1]);

      else
        i = any(N==0, 2);
        % N(i, :) = [c_i 0] or [0 c_i], change to [c_i c_i] to make
        % C(outerFaces) = 0
        N(i, :) = repmat(sum(N(i, :), 2), [1, 2]);
        % face transmissibility
        harm_c  = 1 ./ accumarray(cellFace, 1./opt.Trans, [nif, 1]);
        if isempty(opt.dhfz)
            C       = G.cells.centroids(N(:,2),:) - G.cells.centroids(N(:,1),:);
            gflux   = harm_c.*(C*g');
        else
            %if(false)% work for normal and proper periodic grids
            %    dzf=G.cells.z(N(:,2)) - G.cells.z(N(:,1));
            %else % work for general case
               sgn=2*(G.faces.neighbors(G.cells.faces(:,1),1)==cellNo)-1;
               dzf=accumarray(G.cells.faces(:,1),opt.dhfz.*sgn);
            %end
            gflux=harm_c.*dzf;
            gflux=norm(g)*gflux;
            gflux(i)=0;
        end
      end
   end

   if isfield(fluid, 'pc')
      % Capillary pressure
      i = any(N == 0, 2);

      % N(i, :) = [c_i 0] or [0 c_i], change to [c_i c_i] to make
      % C(innerFaces) = 0
      N(i, :) = repmat(sum(N(i, :), 2), [1, 2]);

      if isempty(opt.Trans)
         C   = G.cells.centroids(N(:,1),:) - G.cells.centroids(N(:,2),:);
         nC  = sqrt(sum(C .* C, 2));

         % Weigthing of permeability for each cellface
         nKC = sum(G.faces.normals(G.cells.faces(:,1), r).* ...
                   K(cellNo,:).* C(G.cells.faces(:,1), c), 2);

         % Permeability contribution for each face by taking the
         % harmonic average of nKC on all faces and divide by |C| = nC
         % Area weighted since n is area weighted.
         harm_c = 2 ./ accumarray(G.cells.faces(:,1), 1 ./ nKC, ...
                                  [G.faces.num, 1]) ./ nC;

         % Two point approximation to grad pc
         d_pc = @(cell_pc) ...
            (cell_pc(N(~i,1),:) - cell_pc(N(~i,2),:)) ./ nC(~i);

         % Flux contribution of capillary pressure for internal faces
         % pc_flux_ij = A_ij*K_avg*(pc(s_i)-pc(s_j))/nC = harm_c*d_pc,
         % A_ij = area of face.
         if isfield(fluid, 'pc') || isfield(fluid, 'pcOW')
            pc_flux = @(rSol) harm_c(~i) .* d_pc(getIncompCapillaryPressure(rSol, fluid));
         else
            pc_flux = @(rSol) zeros(sum(i), 1);
         end

         % The constant contribution to the Jacobian from the term
         % d/ds( K grad pc) which in flux formulation is
         % d/ds_i (harm_c(pc_flux)
         pcJac = harm_c ./ nC;
      else
         % Use transmissibilites
         if numel(opt.Trans) == numel(cellFace)
            harm_c = 1 ./ accumarray(cellFace, 1./opt.Trans, ...
                                     [G.faces.num, 1]);
         elseif numel(opt.Trans) == size(N, 1)
            harm_c = opt.Trans;
         else
            error('Unsupported Size of Transmissibility Field');
         end

         % Two point approximation to grad pc
         d_pc_mod = @(cell_pc) (cell_pc(N(~i,2),:) - cell_pc(N(~i,1),:));

         % Flux contribution of capillary pressure for internal faces
         % pc_flux_ij = A_ij*K_avg*(pc(s_i)-pc(s_j))/nC = harm_c*d_pc,
         % A_ij = area of face.
         pc_flux = @(rSol) harm_c(~i) .* d_pc_mod(fluid.pc(rSol));

         % The constant contribution to the Jacobian from the term
         % d/ds( K grad pc) which in flux formulation is
         % d/ds_i (harm_c(pc_flux)
         pcJac = -harm_c;
      end
   else
      pc_flux = @(rSol) zeros(sum(all(neighborship~=0, 2)), 1);
      pcJac = zeros(size(gflux));
   end
end

%--------------------------------------------------------------------------

function [iw, io] = upwindIndices(G, dflux, gflux, pcflux, mob)
   neighbors = getNeighbourship(G, 'Topological', true);
   intern = all(neighbors ~= 0, 2);
   gflux  = gflux(intern) + pcflux;
   dflux  = dflux(intern);
   iw     = nan(sum(intern), 1);
   io     = nan(sum(intern), 1);


   % Upwind direction for 'water': When v and g have the same sign, the
   %  sign of water phase flux is independent of mobilities.
   a = ~(dflux < 0) & ~(gflux < 0);
   b = ~(dflux > 0) & ~(gflux > 0);
   c = a | b;

   N       = neighbors(intern,:);
   N(b, :) = N(b, [2,1]);
   iw(c)   = N(c, 1);
   clear a b c N

   % Upwind direction for 'oil': When v and g have the opposite sign, the
   %  sign of the oil phase flux is independent of mobilities.
   a = ~(dflux < 0) & ~(gflux > 0);
   b = ~(dflux > 0) & ~(gflux < 0);
   c =  a | b;

   N       = neighbors(intern,:);
   N(b, :) = N(b, [2,1]);
   io(c)   = N(c, 1);
   clear a b c N

   if any(any(isnan([iw, io])))
      mw = nan(size(iw));
      mo = mw;


      mw(~isnan(iw)) = mob(iw(~isnan(iw)), 1);
      mo(~isnan(io)) = mob(io(~isnan(io)), 2);

      vw = dflux + gflux.*mo;
      vo = dflux - gflux.*mw;

      N             = neighbors(intern,:);
      N(vw<0, :)    = N(vw<0, [2,1]);
      iw(isnan(iw)) = N(isnan(iw), 1);

      N             = neighbors(intern,:);
      N(vo<0, :)    = N(vo<0, [2,1]);
      io(isnan(io)) = N(isnan(io), 1);
   end
end
