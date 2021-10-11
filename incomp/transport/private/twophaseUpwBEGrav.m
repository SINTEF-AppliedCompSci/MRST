function [resSol,report] = twophaseUpwBEGrav(resSol, G, tf, q, flux, grav, ...
                                    pv, fluid, varargin)
%Implicit single point upwind solver for two-phase flow, including gravity.
%
% SYNOPSIS:
%   resSol = twophaseUpwBEGrav(resSol, G, tf, q, flux, grav, pv, fluid)
%   resSol = twophaseUpwBEGrav(resSol, G, tf, q, flux, grav, pv, fluid, ...
%                              'pn1', pv1, ...)
%
% DESCRIPTION:
%   Function twophaseUpwBEGrav solves the Buckley-Leverett transport
%   equation
%
%        s_t + f(s)_x = q
%
%   using a first-order upwind discretisation in space and a backward Euler
%   discretisation in time.  The nonlinear system of equations that must be
%   solved to move the solution from time=0 to time=tf, are solved using a
%   Newton-Raphson algorithm with line search to increase robustness.
%
%   In the case of failing to compute a solution using only a single step
%   over the interval [0,tf], an alternative strategy involving sub-steps
%   and step size control is employed.
%
% REQUIRED PARAMETERS:
%   resSol  - Reservoir solution structure containing valid
%             saturation resSol.s with one value for each cell in
%             the grid.
%
%   G       - Grid data structure discretising the reservoir model.
%
%   tf      - End point of time integration interval (i.e., final time),
%             measured in units of seconds.
%
%   q       - Accumulated sources (typically contributions from wells).
%             One scalar value for each cell in the grid.  The source rates
%             are assumed to be measured in units of m^3/s.
%
%   flux    - Inflow matrix of fluxes into each cell, created
%             by function initTransport.  Size
%             G.cells.num-by-(G.faces.num-boundary faces).  The fluxes are
%             assumed to be measured in units of m^3/s.
%
%   grav    - Matrix with gravity contribution for each face, created by
%             function initTransport.  Same size as flux.
%
%   pv      - Reservoir pore volumes, measured in units of m^3.  One scalar
%             value for each cell in the model.
%
%   fluid   - Data structure describing the fluids in the problem.
%
% OPTIONAL PARAMETERS:
%
%   verbose  - Whether or not time integration progress should be reported
%              to the screen.
%              Default value: verbose = false.
%
%   nltol    - Absolute tolerance of iteration.  The numerical solution
%              must satisfy the condition
%
%                NORM(S-S0 + dt/porvol(out - in) - Q, INF) <= nltol
%
%              at all times in the interval [0,t].
%              Default value: nltol = 1.0e-6.
%
%   lstrials - Maximum number of trials in linesearch method.  Each new
%              trial corresponds to halving the step size along the search
%              direction.
%              Default value: lstrials = 20.
%
%   maxnewt  - Maximum number of inner iterations in Newton-Raphson method.
%              Default value: maxnewt = 25.
%
%   tsref    - Maximum time step refinement power.  The minimum time step
%              allowed is tf / 2^tsref.
%              Default value: tsref = 12.
%
% RETURNS:
%   resSol   - Reservoir solution with updated resSol.s.
%
% SEE ALSO:
%   `twophaseUpwBE`, `initTransport`, `implicitTransport`, `explicitTransport`.

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

[F, Jac, linsrch, verbose, ...
 nltol, maxnwt, mints] = setup(tf, q, flux, grav, pv, fluid, varargin{:});

[v_darcy, neighbors, ...
 normals, up_inx, grav_vec] = initFaceMob(G, resSol, flux, grav);

[rho, rho] = fluid.properties(resSol);
findMobM = @(mob) findFaceMobMat(up_inx, mob, grav_vec, v_darcy, ...
                                 neighbors, normals, rho);

report = struct('success',      true,...
                'iterations',   0, ...
                'sub_steps',    0, ...
                'failed_steps', 0);

sz = size(resSol.s(:,1));
sn = resSol;

dispif(verbose, '\n\n');
dispif(verbose, [repmat('-', [1, 70]), '\n']);
dispif(verbose, '  Time interval           iter   relax    residual \n');
dispif(verbose, [repmat('-', [1, 70]), '\n']);

[t, dt, dtprev, count] = deal(0.0, tf, -tf, 0);

%--------------------------------------------------------------------------
% Main time loop (ideally, just a single step) ----------------------------
%
while t < tf && dt > mints,
   dt = min(dt, tf - t);

   % Update faceMob index matrices:
   % NB: A_w, and A_o could be updated inside the inner NR algorithm
   % (before computing the Jacobian matrix) to increase accuracy, but doing
   % so would decrease the efficiency greatly.
   if t == 0 || nnz(grav),
      mob = mobilities(resSol, fluid);
      [A_w, A_o] = findMobM(mob);
   end

   %-----------------------------------------------------------------------
   % Outer controlling loop (sub step size &c) ----------------------------
   %
   redo_newton = true;
   while redo_newton,
      dispif(verbose, '[%f, %f]:', t / tf, (t + dt) / tf); nc = 0;

      s0 = resSol.s(:,1);
      sn = resSol; sn.s(:,1) = 0.5;

      %--------------------------------------------------------------------
      % Initialise NR algorithm -------------------------------------------
      %
      res     = F(sn, s0, dt, A_w, A_o);

      err     = norm(res, inf);
      nwtfail = err > nltol;

      linfail = false;
      it      = 0;
      %--------------------------------------------------------------------
      % Execute inner NR algorithm ----------------------------------------
      %
      while nwtfail && ~linfail && it < maxnwt,

         ds = - Jac(sn, dt, A_w, A_o) \ res;

         [sn, res, alph, linfail] = linsrch(sn, s0, ds, dt, err, A_w, A_o);

         it      = it + 1;
         err     = norm(res, inf);
         nwtfail = err > nltol;

         dispif(verbose, repmat('\b', [1, nc]));
         nc = dispif(verbose, '\t%4d\t %2.2f\t %5.5e', it, alph, err);
         report.iterations = report.iterations +1;
      end

      %--------------------------------------------------------------------
      % Determine if NR succeeded or not ----------------------------------
      %
      count = count - 1;
      if nwtfail,
         dispif(verbose, '\tReducing step.');

         if count > 0 && dtprev > 0,
            dt = dtprev;
         else
            dt = dt / (1.5 + 0.5);
         end
         count = 5;
         report.failed_steps = report.failed_steps + 1;
      else
         redo_newton = false;
         t = t + dt; dtprev = -dt;
         if it == 0,
            dispif(verbose, '\t%4d\t %s\t %5.5e', it, '-', err);
            dispif(verbose, '\t NB: err <= ntol.');
         end

         if it <= 5  && count < 0 && t < tf, % Arbitrary threshold.
            dispif(verbose, '\tIncreasing step.');

            dtprev = dt;
            dt     = min((1 + 0.5) * dt, tf - t);
            count  = 5;
         end
      end
      dispif(verbose, '\n');
   end

   % Time step [t, t+dt] was successful
   report.sub_steps = report.sub_steps + 1;
   resSol = sn;
end

if ~(t < tf),
   dispif(verbose, 'We''re done (%f)\n', t);
else
   dispif(verbose, 'Unable to integrate to final time, t=%f\n', tf);
   resSol.s(:,1) = NaN(sz);
   report.success = false;
end

if size(resSol.s,2) > 1,
   % Update oil saturation:
   resSol.s(:,2) = 1 - resSol.s(:,1);
end
end
%-----------------------------------------------------------------------
% Private helpers follow.
%-----------------------------------------------------------------------

function [F, Jac, linsrch, verbose, nltol, maxnwt, mints] = ...
      setup(tf, q, flux, grav, pv, fluid, varargin)

nc = numel(pv);

prm  = struct('verbose',  false,  ...  % emit progress reports
              'nltol',    1.0e-6, ...  % non-linear residual tolerance
              'lstrials', 20,     ...  % max no of line search trials
              'maxnewt',  25,     ...  % max no. of NR iterations
              'tsref',    12,     ...  % time step refinement
              'resred',   0.99);       % residual reduction factor

prm = merge_options(prm, varargin{:});

verbose = prm.verbose;
nltol   = prm.nltol;
maxnwt  = prm.maxnewt;
mints   = pow2(tf, -prm.tsref);        % minimum time step size

q  = q (:);
pv = pv(:);

% System F(s) = 0 of non-linear equations (and the system's Jacobian
% matrix) defining saturation equilibrium at a given time step.
%
% Here, we use a Buckley-Leverett model with
%
%            s²/µw              s²                µw
%    f(s) = ---------------- = ------------ ,  mr=---
%            s²/µw+(1-s)²/µo    s²+mr*(1-s)²      µo
%
%
%            2 mr*s(1-s)
%   df(s) = ---------------
%           (s²+mr*(1-s)²)²
%
% With the matrices flux, grav, A_w and A_o, we can write an upwind
% backward Euler discretisation of the Buckley-Leverett model as
%
%     s^(n+1) = s^n - (dt./pv)*(H(s^n+1) - max(q,0) - min(q,0)*f(s^n+1))
%
% where H(s) = (flux + grav*diag(A_o*lam_o(s))
%                 *A_w*lam_w(s)./(A_w*lam_w(s)+A_o*lam_o(s)),
%
% lam_l is the mobility for phase l, f is Buckely-Leverett fractional
% flow function, while A_o and A_w are index matrices that determine
% the upstream mobility.
%
% The target function is
%
%    F(s) = s-s^{n-1} + (dt./pv)*(H(s) - max(q, 0)- min(q,0)*f(s))
%
% and the Jacobian is
%
%   dF(s) =  I + (dt./pv)(H'(s)  - min(q,0)*df(s))

   PV   = sparse(1 : nc, 1 : nc, pv);
   flux = PV \ flux;
   grav = PV \ grav;

   in   = PV \ max(q, 0);
   out  = PV \ min(q, 0);

   % Return function handles
   F   = @(varargin) Residual(varargin{:}, flux, grav, in, out, fluid);
   Jac = @(varargin) Jacobian(varargin{:}, flux, grav,     out, fluid);

   % Bind prm.resred * err, @(sat) F(sat, s0, dt) and prm.lstrials
   % to function call to reduce clutter in main code.
   linsrch = @(state, s0, ds, dt, err, A_w, A_o)        ...
               linesearch(state, ds, prm.resred * err,  ...
                          @(state) F(state, s0, dt, A_w, A_o), ...
                          prm.lstrials);
end

%--------------------------------------------------------------------------

function F = Residual(state, s0, dt, Aw, Ao, flux, grav, in, out, fluid)
   mob  = mobilities(state, fluid);
   mobf = [Aw * mob(:,1), Ao * mob(:,2)];

   % Evaluate fractional flow function
   f  = mob (:,1) ./ sum(mob , 2);   % In cells.
   ff = mobf(:,1) ./ sum(mobf, 2);   % Across faces (upwind).

   nf = size(mobf, 1);
   M  = flux + grav*sparse(1 : nf, 1 : nf, mobf(:,2));

   s = fluid.saturation(state);
   F = s(:,1) - s0 + dt.*(M*ff - out.*f - in);
end

%--------------------------------------------------------------------------

function J = Jacobian(state, dt, Aw, Ao, flux, grav, out, fluid)
   [nf, nc] = size(Aw);

   [mob, dmob] = mobilities(state, fluid);

   mobf = [Aw * mob(:,1), Ao * mob(:,2)];
   dw   = Aw * sparse(1 : nc, 1 : nc, dmob(:,1));
   do   = Ao * sparse(1 : nc, 1 : nc, dmob(:,2));

   Lt  = sum(mob, 2);
   f   = bsxfun(@rdivide, mob, Lt);
   df  = (f(:,2).*dmob(:,1) - f(:,1).*dmob(:,2)) ./ Lt;

   Ltf = sum(mobf, 2);
   ff  = bsxfun(@rdivide, mobf, Ltf);
   dff = sparse(1 : nf, 1 : nf, ff(:,2) ./ Ltf) * dw - ...
         sparse(1 : nf, 1 : nf, ff(:,1) ./ Ltf) * do;

   M = flux + grav*sparse(1 : nf, 1 : nf, mobf(:,2));
   J = speye(nc) + ...
       dt .* (M    * dff                              + ...
              grav * sparse(1:nf, 1:nf, ff(:,1)) * do - ...
              sparse(1 : nc, 1 : nc, out .* df));
end

%--------------------------------------------------------------------------

function varargout = mobilities(state, fluid)
   mu  = fluid.properties(state);
   s   = fluid.saturation(state);
   [kr{1:nargout}] = fluid.relperm(s, state);

   varargout = cellfun(@(n) bsxfun(@rdivide, n, mu), kr, ...
                       'UniformOutput', false);
end

%--------------------------------------------------------------------------

function [resSol, res, alph, fail] = linesearch(resSol, ds, target, F, ni)
%
% Basic idea: search for a step size 'alpha' in direction 'ds', subject to
% the restriction that alpha be in [0,1], which ensures that the objective
% function 'F' decreases.  That is: F(s + alpha*ds) < F(s).
%
% In the current implementation, alpha is reduced in a geometric sequence.
% A more sophisticated approach would ensure a certain minimum reduction as
% well.
%
   minSat = 0;     % Minimum (water) saturation
   maxSat = 1;     % Maximum (water) saturation
   capSat = @(sat) min(max(minSat, sat), maxSat);

   alph = 0;
   i    = 0;
   fail = true;

   % Geometric line search: seems pretty robust
   while fail && (i < ni),
      sn.s = capSat(resSol.s(:,1) + pow2(ds, alph));
      res  = F(sn);

      alph = alph - 1;
      i    = i + 1;
      fail = ~(norm(res, inf) < target);
   end

   alph = pow2(alph + 1);      % Undo last (unneeded) scaling.
   resSol.s(:,1) = sn.s;
end


function [A_w, A_o] = ...
      findFaceMobMat(u_inx, mob, g_vec, v_darcy, neighbors, ...
                     normals, rho)
%Initialize upwind saturation index matrices for face mobility.
%Uses known cell mobilities and output from function initFaceMob.
%
% SYNOPSIS:
%   [A_w, A_o] = findFaceMobMat(u_inx, mob, g_vec, v_darcy, ...
%                               neighbors, normals, fluid)
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
%   fluid     - Data structure describing the fluids in the problem.
%
% RETURNS:
%   A_w       - Index matrix for converting from cell mobility to face
%               mobility for water-phase (or primary phase).
%               facemob_w = A_w * mob.
%
%   A_o       - Index matrix for converting from cell mobility to face
%               mobility for oil-phase (or secondary phase).
%
% SEE ALSO:
%   `initFaceMob`, `twophaseUpwFEGrav`, `twophaseUpwBEGrav`.

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
%   face. These arrays are used to build the upwind saturation index
%   matrices A_w and A_o.

   faceMob = mob(u_inx, :);

   h_inx = u_inx;
   l_inx = u_inx;

   % Recognize the heaviest (h) and lightest (l -- ell) phase.
   if rho(1) > rho(2),
      h = 1; l = 2;
   else
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

   % Initialize upwind saturation index matrices
   n = numel(w_inx);
   A_w = sparse(1:n, double(w_inx), 1, n, size(mob,1));

   if o_inx == w_inx
      A_o = A_w;
   else
      A_o = sparse(1:n, double(o_inx), 1, n, size(mob,1));
   end
end
