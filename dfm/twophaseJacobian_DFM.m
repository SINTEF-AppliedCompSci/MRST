function [Res, Jac] = twophaseJacobian_DFM(G, state, rock, fluid, varargin)
%Residual and Jacobian of single point upwind solver for two-phase flow.
%
% SYNOPSIS:
%   resSol = twophaseJacobian_DFM(G, state, rock, fluid)
%   resSol = twophaseJacobian_DFM(G, state, rock, fluid, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Function twophaseJacobian returns function handles for the residual
%   and its Jacobian matrix for the implicit upwind-mobility weighted
%   dicretization of
%              __
%        s_t + \/· [f(s)(v·n + mo(rho_w - rho_o)n·Kg)] = f(s)q
%
%
%   where v·n is the sum of the phase Dary fluxes, f is the fractional
%   flow function,
%
%                  mw(s)
%        f(s) = -------------
%               mw(s) + mo(s)
%
%   mi = kr_i/mu_i is the phase mobiliy of phase i, mu_i and rho_i are the
%   phase viscosity and density, respectivelym, g the (vector) acceleration
%   of gravity and K the permeability.  The source term f(s)q is a
%   volumetric rate of water.
%
%   Using a first-order upwind discretisation in space and a backward Euler
%   discretisation in time,  the residual of the nonlinear system of
%   equations that must be solved to move the solution state.s from time=0
%   to time=tf, are obtained by calling F(state, s0, dt) which yields
%   Likewise, the Jacobian matrix is obtained using the function Jac.
%
% NOTE:
%   This file has been modified from the original twophaseJacobian to
%   account for non-neighbor connections, as described by the field
%   G.cells.neighbors. The modification is intended for grids with hybrid
%   cells, but might be of use for other applications as well. Also note
%   that all modifications needed to apply the implicit transport solver to
%   a grid with cell-to-cell connections are done within this file.
%
% REQUIRED PARAMETERS:
%   resSol  - Reservoir solution structure containing valid
%             saturation resSol.s with one value for each cell in
%             the grid.
%
%   G       - Grid data structure discretising the reservoir model.
%
%   rock    - Struct with fields perm and poro.  The permeability field
%             ('perm') is only referenced when solving problems involving
%             effects of gravity and/or capillary pressure and need not be
%             specified otherwise.
%
%   fluid   - Data structure describing the fluids in the problem.
%
% OPTIONAL PARAMETERS:
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
%   - implement gravity effects for hybrid cells
%   - implement capillary pressure effects for hybrid cells
%{
Copyright 2009, 2010, 2011 SINTEF ICT, Applied Mathematics.

Portions Copyright 2011-2012 University of Bergen.

This file is part of DFM module of The MATLAB Reservoir Simulation Toolbox
(MRST).

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
%        s_t + \/· [D·f(s) + D·f(s)·mo(s)] = f(s) Q
%
% for constant vector and scalar fields fields D, G and Q.


   opt = struct('verbose', mrstVerbose, 'wells', [], 'src', [], 'bc', []);
   opt = merge_options(opt, varargin{:});

   assert (size(state.s,2)<3 || all(state.s(:,3)==0));


   % All source terms, i.e., boundary conditions, source terms and wells.
   q = computeTransportSourceTerm(state, G, opt.wells, opt.src, opt.bc);

   % Extract (constant) fluid densities.
   [rho, rho]     = fluid.properties(state);

   % Compute the gravitational potential (rho_w - rho_o)n·Kg across each
   % face in the grid.
   [gflux, pc_flux, pcJac, gfluxc2c] = getFlux(G, rock, rho, fluid);

   % Bind flux terms state.flux and gflux.  The call
   % [iw,io]=findUpwindCells(..) computes cell index of upwind cell for the
   % water phase (iw) and oil phase (io).
   findUpwindCells = @(mob, state) upwindIndices(G, state.flux, state.fluxc2c, gflux,gfluxc2c, ...
                                                 pc_flux(state), mob);

   pv = poreVolume(G, rock);

   % Return function handles
   Res = @(xr, s0, dt) Residual(xr, s0, dt, fluid, state.flux, state.fluxc2c,gflux, gfluxc2c, pc_flux,...
                                findUpwindCells, G, pv, q);

   Jac = @(xr, dt)     Jacobian(xr, dt,     fluid, state.flux, state.fluxc2c,gflux, gfluxc2c,pc_flux, pcJac, ...
                                findUpwindCells, G, pv, q);
end

%--------------------------------------------------------------------------
function J = Jacobian (resSol, dt, fluid, dflux, dflux2,gflux,gflux2, pcflux, pcJac, ...
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

   mu        = fluid.properties(resSol);     % viscosity
   sat       = fluid.saturation(resSol);     % cell saturation
   [kr, dkr] = fluid.relperm(sat, resSol);
   m         = bsxfun(@rdivide, kr,  mu);    % mobility in each cell

   % We need derivatives with respect to s(:,1) only.
   dkr(:,end) = -dkr(:,end);

   % Mobility derivative.  Two-phase only.  Extract diagonal of Jacobian.
   dm = bsxfun(@rdivide, dkr(:, [1, end]), mu);

   clear mu sat kr dkr

   internal   = all(G.faces.neighbors~=0, 2);
   ic1  = double(G.faces.neighbors(internal,1));
   ic2  = double(G.faces.neighbors(internal,2));
   clear internal_faces

   % add for cell-cell connections
   if isfield(G.cells,'neighbors')
    hic = G.cells.neighbors;
    ic1 = [ic1 ; hic(:,1)];
    ic2 = [ic2 ; hic(:,2)];
   end

   dflux=[dflux(internal);dflux2];
   gflux=[gflux(internal);gflux2]+pcflux(resSol);

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
   dvw_ds1 =  (dflux + (gflux) ...
              .*m_face(:,2)) .* dmw_fo ./ sum(m_face,2);
   dvw_ds2 = -(dflux - (gflux) ...
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

  if isfield(fluid, 'pc')
      % added:
      [dpc, dpc]       = fluid.pc(resSol);

      d_p  =      dt.* f_face(:,1).*m_face(:,2).*pcJac;

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
function F = Residual (resSol, s0, dt, fluid, dflux, dflux2,gflux, gflux2, pcflux,...
                       findPhaseUpwindCells, G, pv, q)
%  DISCRETISATION
%
%  Compute  F given by a residual function Fk for each grid cell
%
%   Fk(s, s0, dt) = s - s0 + dt/pv·[H(s) - max(q,0) - min(q,0)·f],
%
%   H(s) = sum_i [f_i(dflux_i + mo_i*gflux_i)]  (faces i of cell)
%
%  where f_i = mw_i/(mw_i+mo_i), mo_i and mw_i are phase upwind
%  mobilities at face i and dflux and gflux are Darcy flux and gravity flux
%  g(rho1-rho2)*face_normal, respectively. The other flux function f is
%  the fractional flow function evaluated with cell mobilities.
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
%  The upwind cell for each phase at each internal grid face are found
%  using findPhaseUpwindCells, and computeConstData.  The gravity flux is
%  computed in getFlux.

   % Compute cell mobilities
   mu  = fluid.properties(resSol);
   sat = fluid.saturation(resSol);
   kr  = fluid.relperm(sat, resSol);
   m   = bsxfun(@rdivide, kr, mu);
   clear mu sat kr


   internal = all(G.faces.neighbors~=0, 2);
   ic1  = G.faces.neighbors(internal,1);
   ic2  = G.faces.neighbors(internal,2);
   clear internal_faces

   % Modified to allow cell2cell connections
   if isfield(G.cells,'neighbors')
       hic = G.cells.neighbors;
       ic1 = [ic1 ; hic(:,1)];
       ic2 = [ic2 ; hic(:,2)];
   end

   dflux=[dflux(internal);dflux2];
   gflux=[gflux(internal);gflux2];

   % Compute water source term
   f_cell = m(:,1) ./ sum(m,2);
   Q      = max(q,0) + min(q,0).*f_cell;
   clear f_cell

   % Compute the upwind mobility weighted flow of water:
   % For each face k, iw(k) and io(k) is the cell index in which we
   % must evaluate the water and oil mobilities, respectively.
   [iw, io] = findPhaseUpwindCells(m, resSol);
   if ~any(gflux) && ~isfield(fluid, 'pc'),
      assert(all(iw == io));
   end
   m_face  = [m(iw, 1) m(io,2)];
   f_face  = m_face(:,1)./sum(m_face,2);

   v_water = f_face.*(dflux+(gflux + pcflux(resSol)).*(m_face(:,2)));
   clear io iw m m_face f_face

   F   = resSol.s(:,1) - s0(:,1);
   F   = F + dt.*(1./pv).*(accumarray([ic1; ic2], [v_water; -v_water], ...
                                    size(F)) - Q);

end


function [gflux, pc_flux, pcJac, gfluxc2c]  = getFlux(G, rock, rho, fluid)
% -------------------------------------------------------------------------
%
%                        TECHNICAL DESCRIPTION
%
% -------------------------------------------------------------------------
%
%  gflux = harmonic average of (n·K·gravity·(rho1-rho2)) in each cell-

   g  = gravity();
   g  = g(1:G.griddim) * (rho(1) - rho(2));
   gflux  = zeros([G.faces.num, 1]);

   % for the moment no gravity flux in the hybrid cells.
   if isfield(G.cells,'neighbors')
       gfluxc2c = zeros(size(G.cells.neighbors,1),1);
   else
       gfluxc2c = [];
   end

   dim = size(G.nodes.coords, 2);

   if (norm(g) > 0) || isfield(fluid, 'pc'),
      [K, r, c] = permTensor(rock, dim);
      nc        = G.cells.num;
      cellNo    = rldecode(1 : nc, diff(G.cells.facePos), 2) .';
      if isfield(G.cells,'neighbors')
          nh = length(G.hybridNeighbors.n);
          hc = rldecode((1:nh)',diff(G.hybridNeighbors.facePos));
          cellNo2 = G.faces.neighbors(G.hybridNeighbors.faces,1);
      end
   end

   if norm(g) > 0,
      % nKg == n' * K * g *(rho1 - rho2) for all cellfaces.
      nKg    = sum(G.faces.normals(G.cells.faces(:,1), r) .* ...
                   bsxfun(@times, K(cellNo,:), g(c)), 2);
      % Compute harmonic average of nKg(rho1-rho2) on all faces.
      gflux = 2 ./ accumarray(G.cells.faces(:,1), 1 ./ nKg, [G.faces.num, 1]);

      % This part is not thoroughly tested
      if isfield(G.cells, 'neighbors')
          nKg_hyb = sum(G.faces.normals(G.hybridNeighbors.faces, r) .* ...
              bsxfun(@times, K(cellNo2,:), g(c)), 2);
          sumnKg_hyb = accumarray(hc, abs(nKg_hyb), [nh, 1]);
          gfluxc2c = 2.*prod(nKg_hyb(G.hybridNeighbors.neighbors).*(ones(length(gfluxc2c),1)*[1,-1]),2)./rldecode(sumnKg_hyb,G.hybridNeighbors.n);
      end

   end

   if isfield(fluid, 'pc')
      % Capillary pressure
      N = G.faces.neighbors;
      i = any(N==0, 2);

      % N(i, :) = [c_i 0] or [0 c_i], change to [c_i c_i] to make
      % C(innerFaces) = 0
      N(i, :) = repmat(sum(N(i, :), 2), [1, 2]);
      C   = G.cells.centroids(N(:,1),:) - G.cells.centroids(N(:,2),:);
      nC  = sqrt(sum(C.*C, 2));

      % Weigthing of permeability for each cellface
      nKC = sum(G.faces.normals(G.cells.faces(:,1), r).* ...
                K(cellNo,:).* C(G.cells.faces(:,1), c), 2);

      % Permeability contribution for each face by taking the
      % harmonic average of nKC on all faces and divide by |C| = nC
      % Area weighted since n is area weighted.
      harm_c = 2 ./ accumarray(G.cells.faces(:,1), 1./ nKC, [G.faces.num, 1])./nC;

      N = N(~i, :);
      nC = nC(~i);
      harm_c = harm_c(~i);

      if isfield(G.cells,'neighbors')

          C2 = G.cells.centroids(cellNo2,:) - G.faces.centroids(G.hybridNeighbors.faces,:);
          nC2  = sqrt(sum(C2.*C2, 2));
          nC2 = sum(nC2(G.hybridNeighbors.neighbors),2);
          nKC2 = sum(G.faces.normals(G.hybridNeighbors.faces, r).* ...
                K(cellNo2,:).* C2(:,c), 2);
          sumnKC2= accumarray(hc, nKC2, [nh, 1]);
          harm_c2 = prod(reshape(nKC2(G.hybridNeighbors.neighbors),sum(G.hybridNeighbors.n),2),2)./rldecode(sumnKC2,G.hybridNeighbors.n);
          harm_c2 = 2.* harm_c2 ./ nC2;

          N2 = G.cells.neighbors;
          N = [N ; N2];
          nC = [nC ; nC2];
          harm_c = [harm_c; harm_c2];
      end


      % Two point approximation to grad pc
      d_pc   = @(cell_pc) ( cell_pc(N(:,1),:)-cell_pc(N(:,2),:))./nC;

      % Flux contribution of capillary pressure for internal faces
      % pc_flux_ij = A_ij*K_avg*(pc(s_i)-pc(s_j))/nC = harm_c*d_pc,
      % A_ij = area of face.
      pc_flux = @(rSol)  harm_c.*d_pc(fluid.pc(rSol));

      % The constant contribution to the Jacobian from the term
      % d/ds( K grad pc) which in flux formulation is
      % d/ds_i (harm_c(pc_flux)
      pcJac  = harm_c./nC;

   else
       num = sum(all(G.faces.neighbors~=0,2));
       % Modified to allow cell2cell connections
       if isfield(G.cells,'neighbors')
           num = num + length(G.cells.neighbors);
       end
      pc_flux = @(rSol) zeros(num, 1);
      pcJac = zeros(size(gflux));
   end
end

%--------------------------------------------------------------------------
function [iw, io] = upwindIndices(G, dflux, dfluxc2c, gflux, gfluxc2c,pcflux, mob)
   intern = all(G.faces.neighbors ~= 0, 2);

   gflux  = [gflux(intern) ; gfluxc2c] + pcflux;
   dflux  = [dflux(intern) ; dfluxc2c];
   iw     = nan(length(gflux), 1);
   io     = nan(length(gflux), 1);


   %% Upwind direction for 'water': When v and g have the same sign, the
   %  sign of water phase flux is independent of mobilities.
   a = ~(dflux < 0) & ~(gflux < 0);
   b = ~(dflux > 0) & ~(gflux > 0);
   c = a | b;

   intern = all(G.faces.neighbors ~= 0, 2);
   N       = G.faces.neighbors(intern,:);

   % Modified to allow cell2cell connections
   if isfield(G.cells,'neighbors')
       N = [N ; G.cells.neighbors];
   end

   N(b, :) = N(b, [2,1]);
   iw(c)   = N(c, 1);
   clear a b c N

   %% Upwind direction for 'oil': When v and g have the opposite sign, the
   %  sign of water phase flux is independent of mobilities.
   a = ~(dflux < 0) & ~(gflux > 0);
   b = ~(dflux > 0) & ~(gflux < 0);
   c =  a | b;

   N       = G.faces.neighbors(intern,:);

   % Modified to allow cell2cell connections
   if isfield(G.cells,'neighbors')
       N = [N ; G.cells.neighbors];
   end
   N(b, :) = N(b, [2,1]);
   io(c)   = N(c, 1);
   clear a b c N

   if any(any(isnan([iw, io]))),
      mw = nan(size(iw));
      mo = mw;


      mw(~isnan(iw)) = mob(iw(~isnan(iw)), 1);
      mo(~isnan(io)) = mob(io(~isnan(io)), 2);

      vw = dflux + gflux.*mo;
      vo = dflux - gflux.*mw;

      N             = G.faces.neighbors(intern,:);

      % Modified to allow cell2cell connections
      if isfield(G.cells,'neighbors')
          N = [N ; G.cells.neighbors];
      end
      N(vw<0, :)    = N(vw<0, [2,1]);
      iw(isnan(iw)) = N(isnan(iw), 1);

      N             = G.faces.neighbors(intern,:);

      % Modified to allow cell2cell connections
      if isfield(G.cells,'neighbors')
          N = [N ; G.cells.neighbors];
      end
      N(vo<0, :)    = N(vo<0, [2,1]);
      io(isnan(io)) = N(isnan(io), 1);
   end
end
