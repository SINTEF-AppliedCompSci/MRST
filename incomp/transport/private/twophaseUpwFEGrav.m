function [resSol, report] = twophaseUpwFEGrav(resSol, G, tf, q, flux, grav, ...
                                    pv, fluid,  varargin)
%Explicit single point upwind solver for two-phase flow, including gravity.
%
% SYNOPSIS:
%   resSol = twophaseUpwFEGrav(resSol, G, tf, q, flux, grav, porvol,
%                              fluid)
%   resSol = twophaseUpwFEGrav(resSol, G, tf, q, flux, grav, porvol,
%                              fluid, G, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Function twophaseUpwFEGrav solves the Buckley-Leverett transport
%   equation
%
%        s_t + f(s)_x = q
%
%   using a first-order upwind discretisation in space and a forward Euler
%   discretisation in time.  The transport equation is solved on the time
%   interval [0,tf].
%
%   The upwind forward Euler discretisation of the Buckley-Leverett model
%   can be written as:
%
%     s^(n+1) = s^n - (dt./pv)*((H(s^n) - max(q,0) - min(q,0)*f(s^n))
%
%   where
%        H(s) = (flux + grav*diag(A_o*lam_o(s))
%               *A_w*lam_w(s)./(A_w*lam_w(s)+A_o*lam_o(s)),
%
%   pv is the porevolume, lam_l is the mobility for face l, f is
%   Buckely-Leverett fractional flow function, while A_o and A_w are index
%   matrices that determine the upstream mobility.
%
% PARAMETERS:
%   resSol  - Reservoir solution structure containing valid water
%             saturation resSol.s(:,1) with one value for each cell
%             in the grid.
%
%   G       - Grid data structure discretising the reservoir model.
%
%   tf      - End point of time integration interval (i.e., final time),
%             measured in units of seconds.
%
%   q       - Accumulated sources (e.g., contributions from wells and
%             boundary conditions).  One scalar value for each cell in the
%             model.  The source rates are assumed to be measured in units
%             of m^3/s.
%
%   flux    - Inflow matrix of fluxes into each cell, created
%             by function initTransport.  Size
%             G.cells.num-by-(G.faces.num-boundary faces)
%
%             Specifically, the entry flux(i,j) is the flux over
%             (active) face j from cell i.  The fluxes are assumed to be
%             measured in units of m^3/s.
%
%   grav    - Matrix with gravity contribution for each face, created by
%             function initTransport.  Same size as flux.
%
%   porvol  - Reservoir pore volumes, measured in units of m^3.  One scalar
%             value for each cell in the model.
%
%   fluid   - Data structure describing the fluids in the problem.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%
%               - verbose - Whether or not time integration progress should
%                           be reported to the screen.
%                           Default value: verbose = false.
%
%               - dt      - Internal timesteps, measured in units of
%                           seconds, for instance determined by function
%                           initTransport.  Default value = tf.
%                           NB: The explicit scheme is only stable provided
%                           that dt satisfies a CFL time step restriction.
%
%               - satWarn - Tolerance level for saturation warning.
%                            Default value: satWarn = sqrt(eps).
% RETURNS:
%   resSol - Reservoir solution with updated saturations, resSol.s.
%
% SEE ALSO:
%   `twophaseUpwFE`, `initTransport`, `explicitTransport`, `implicitTransport`.

% TODO:
%   - implement gravity effects for pressure boundary and wells

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


assert (size(resSol.s,2)<3 || all(resSol.s(:,3)==0));

% Verify that the flux matrix is correctly sized.  The output from
% function initTransport depends on presence or absence of gravity effects.
%
if size(flux) ~= size(grav),
   error('twophaseUpwFEGrav:InputSize:flux',                    ...
        ['Use twophaseUpwFEGrav only in presence of gravity, ', ...
         'else use twophaseUpwFE.'])
end

opt = struct('dt', tf, 'verbose', false, 'satWarn', sqrt(eps));
opt = merge_options(opt, varargin{:});
dt  = opt.dt;
dt  = min(dt, tf);

assert (dt > 0 && tf > 0);

[v_darcy, neighbors, normals, upw_inx, g_vec] = initFaceMob(G, resSol, ...
                                                           flux, grav);

H = @(f,dz,dt) ((dt ./ pv).*(dz - max(q,0) - min(q,0).*f));
t = 0;

if opt.verbose,
%   h = waitbar(0, ['Solving transport using ' ...
%                   num2str(ceil(tf/dt)), ' time steps...']);
end

report.timesteps = ceil(tf/dt);

[mu, rho] = fluid.properties(resSol);

while t < tf,
   % Compute cell mobility and fractional flow.
   sat = fluid.saturation(resSol);
   kr  = fluid.relperm(sat, resSol);

   mob = bsxfun(@rdivide, kr, mu);
   f_w = mob(:,1) ./ sum(mob,2);

   % Initialize face mobility.
   if nnz(grav)
      [iw, io] = findFaceMobIx(upw_inx, mob, g_vec, ...
                                  v_darcy, neighbors, normals, rho);
      faceMob = [mob(iw,1) mob(io,2)];
   else
      faceMob = mob(upw_inx, :);
   end

   fw_face = faceMob(:,1) ./ sum(faceMob,2);
   %Remove possible NaNs
   fw_face(sum(faceMob,2)==0) = 0;

   dz = flux*fw_face + grav*(faceMob(:,2) .* fw_face);

   % Explicit computation of saturation
   resSol.s(:,1) = resSol.s(:,1) - H(f_w, dz, dt);

   % Correct possible bad saturations
   if true %opt.verbose,
      %waitbar(t/tf,h);
      sMax = max(resSol.s(:,1));
      if sMax> 1 + opt.satWarn,
         inx = find(resSol.s(:,1) > 1+opt.satWarn);
         disp('Saturation exceeds 1 in cells from:')
         fprintf(1,'%d %g\n', [inx, resSol.s(inx,1)]');
      end
      sMin = min(resSol.s(:,1));
      if sMin < -opt.satWarn,
         inx = find(resSol.s(:) < -opt.satWarn);
         disp('Saturation less than 0 in cells from:')
         fprintf(1,'%d %g\n', [inx, resSol.s(inx,1)]');
      end
   end
   resSol.s(resSol.s(:,1)>1,1) = 1 - eps;
   resSol.s(resSol.s(:,1)<0,1) = 0;

   if isfield(resSol, 'minSat')
      % Save minimum saturation for use in modeling of relative
      % permeability hysteresis
      resSol.minSat = min(resSol.s(:,1), resSol.minSat);
   end

   t  = t + dt;
   dt = min(dt, tf - t);
end

%if opt.verbose, close(h), end

if size(resSol.s,2) > 1,
   % Update oil saturation:
   resSol.s(:,2) = 1 - resSol.s(:,1);
end

end




function [w_inx, o_inx] = findFaceMobIx( ...
                     u_inx, mob, g_vec, v_darcy, neighbors, normals, rho)
%Initialize upwind saturation index matrices for face mobility.
%Uses known cell mobilities and output from function initFaceMob.
%
% SYNOPSIS:
%   [iw, io] = findFaceMobMat(u_inx, mob, g_vec, v_darcy, ...
%                               neighbors, normals)
%
%
% PARAMETERS:
%   u_inx     - Index to cells for initial upwind face mobility not
%               considering gravity effects. Size = (#internal faces).
%
%   mob       - Cell mobility. (mob = fluid.mob(resSol)).
%
%   v_darcy   - Darcy flux for internal faces.
%
%   neighbors - Neighbors for internal faces.
%
%   normals   - Vector of face normals for internal faces.
%
%   g_vec     - Gravity contributions = abs[K(rho_1-rho_2)*g*n_z], for
%               each face, where n_z is the z-component of the normal.
%
% RETURNS:
%   iw       - Index for converting from cell mobility to face
%               mobility for water-phase (or primary phase).
%               facemob_w = mob(iw,1).
%
%   io       - Index for converting from cell mobility to face
%               mobility for oil-phase (or secondary phase).
%
% SEE ALSO:
%   `initFaceMob`, `twophaseUpwFEGrav`, `twophaseUpwBEGrav`.

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

% -------------------------------------------------------------------------
%
%                        TECHNICAL DESCRIPTION
%
% -------------------------------------------------------------------------
%
%   Given to cells i, j sharing a face e_ij we need to decide whether to
%   use the saturation from cell i (s_i) or cell j (s_j) when calculating
%   the mobility for e_ij. Let v_ij be the darcy flux over e_ij from cell i
%   to cell j. For cases without gravity, upstream weighting means that we
%   define the saturation of e_ij as
%
%           s(e_ij) = { s_i   if v_ij >= 0,
%                     { s_j   if v_ij <  0.
%
%   However, for a two-phase problem with gravity, the phase fluxes v_1 and
%   v_2 need not have the same direction as the Darcy flux: either v_1 or
%   v_2 can have oposite direction. Thus, in order to determine the
%   upstream mobility for each edge we need to know sign(v_1) and
%   sign(v_2).
%
%   Let f_k be the fractional flow of phase k and let
%      f_k(e_ij) := f_k(s_k(e_ij)), and  mob_k(s_k(e_ij)) := mob_k(e_ij).
%   Furhtermore let z be the depth and direction of gravity. The phase
%   fluxes for phases 1 and 2 over e_ij is given as
%      v_1(e_ij) = f_1(e_ij)*[v + K*(rho_1-rho_2)* mob_2(e_ij) *g*grad(z)]
%      v_2(e_ij) = f_2(e_ij)*[v - K*(rho_1-rho_2)* mob_1(e_ij) *g*grad(z)]
%
%   Observe that
%      sign(v_1) = sign[v + K*(rho_1-rho_2)* mob_2(e_ij) *g*grad(z)].
%
%   Since this equation involves mob_2(e_ij), it means that we need to know
%   the upstream face mobility for one of the phases in order to calculate
%   the sign of phase velocity for the other. We will explain how this is
%   possible below. First, let v_h be the phase velocity of the heaviest
%   phase and v_l the lightest phase (thus rho_h > rho_l), and let n_z be
%   the z-component of the normal pointing from cell i to j. We initialize
%   both phases with upstream weighting given by the Darcy velocity
%   (u_inx), and therefore need to change weighting for the faces where
%   sign(phase velocity)~=sign(Darcy velocity).
%
%   The possible scenarios for an edge e_ij are:
%
%   1. v >= 0 & n_z > 0  =>  v_h >= 0
%      sign(v_l) = sign(v - g_vec(e_ij)*mob_h(i)), if v_l < 0: change.
%   2. v <  0 & n_z > 0  =>  v_l <  0
%      sign(v_h) = sign(v + g_vec(e_ij)*mob_l(j)), if v_h > 0: change.
%   3. v >= 0 & n_z < 0  =>  v_l >= 0
%      sign(v_h) = sign(v + g_vec(e_ij)*mob_l(i)), if v_h < 0: change.
%   4. v <  0 & n_z < 0  =>  v_h <  0
%      sign(v_l) = sign(v - g_vec(e_ij)*mob_h(j)), if v_l > 0: change.
%
%   Notice that faces where n_z == 0 are not included since g_vec == 0 for
%   such faces. The result of the alogrithm is two arrays of cell-numbers,
%   h_inx and l_inx, giving the index to the upwind weighting for each
%   face.

   faceMob = mob(u_inx, :);

   h_inx = u_inx;
   l_inx = u_inx;

   if rho(1)-rho(2) > 0 % 1. phase is heaviest
         h = 1; l = 2;
   else % 2. phase is heaviest
         h = 2; l = 1;
   end



   if nnz(g_vec),

      pos_v = ~(v_darcy < 0);  % v_darcy >= 0
      pos_n =   normals > 0;
      neg_n = ~pos_n;

      if norm(gravity()) > 0,
         inx1 =  pos_v & pos_n;
         inx2 = ~pos_v & pos_n;
         inx3 =  pos_v & neg_n;
         inx4 = ~pos_v & neg_n;
      else
         inx1 =  pos_v & neg_n;
         inx2 = ~pos_v & neg_n;
         inx3 =  pos_v & pos_n;
         inx4 = ~pos_v & pos_n;
      end

      % Recognize the heaviest phase:
      if rho(1) > rho(2), % 1. phase is heaviest
         h = 1; l = 2;
      else % 2. phase is heaviest
         h = 2; l = 1;
      end

      %Four possible scenarios for a face:
      %1: v >= 0, n > 0
      if any(inx1)
         faces = find(inx1);
         cells2 = neighbors(inx1,2);
         a = (v_darcy(inx1) - faceMob(inx1,h).*g_vec(inx1)) < 0;
         l_inx(faces(a)) = cells2(a);
      end

      %2: v < 0, n > 0
      if any(inx2)
         faces = find(inx2);
         cells1 = neighbors(inx2,1);
         a = (v_darcy(inx2) + faceMob(inx2,l).*g_vec(inx2)) > 0;
         h_inx(faces(a))= cells1(a);
      end

      %3: v >= 0, n < 0
      if any(inx3)
         faces = find(inx3);
         cells2 = neighbors(inx3,2);
         a = (v_darcy(inx3) + faceMob(inx3,l).*g_vec(inx3)) < 0;
         h_inx(faces(a)) = cells2(a);
      end

      %4: v < 0 , n < 0
      if any(inx4)
         faces = find(inx4);
         cells1 = neighbors(inx4,1);
         a = (v_darcy(inx4) - faceMob(inx4,h).*g_vec(inx4)) > 0;
         l_inx(faces(a)) = cells1(a);
      end

   end %end if g_vec

   if h == 1 % Water (or the primary phase) is the heaviest phase.
      w_inx = h_inx;
      o_inx = l_inx;
   else
      w_inx = l_inx;
      o_inx = h_inx;
   end
end


