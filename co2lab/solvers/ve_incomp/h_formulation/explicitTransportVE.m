function [state, report] = explicitTransportVE(state, G_top, tf, rock,...
                                                     fluid,  varargin)
%Explicit single point upwind solver for two-phase flow using VE equations.
%
% SYNOPSIS:
% [state, dt_v] = explicitTransportVE(state, G_top, tf, rock, fluid)
% [state, dt_v] = explicitTransportVE(state, G_top, tf, rock,...
%                                                     fluid,  'pn1', pv1)
%
% DESCRIPTION:
%   Function explicitTransportVE solves the Buckley-Leverett transport
%   equation
%
%        h_t + f(h)_x = q
%
%   using a first-order upwind discretisation in space and a forward Euler
%   discretisation in time.  The transport equation is solved on the time
%   interval [0,tf].
%
%   The upwind forward Euler discretisation of the Buckley-Leverett model
%   for the Vertical Equilibrium model can be written as:
%
%     h^(n+1) = h^n - (dt./pv)*((H(h^n) - max(q,0) - min(q,0)*f(h^n))
%
%   where
%        H(h) = f_up(h)(flux + grav*lam_nw_up*(z_diff+rho_diff*h_diff(h)))
%
%   z_diff, h_diff are two point approximations to grad_x z, and grad_x h,
%   f_up and lam_nw_up are the Buckely-Leverett fractional flow function
%   and the mobility for the non-wetting phase, respectively, evaluated for
%   upstream mobility:
%
%        f_up = *A_w*lam_w(h)./(A_w*lam_w(h)+A_nw*lam_nw(h))
%        lam_nw_up = diag(A_nw*lam_nw(h)
%
%   pv is the porevolume, lam_x is the mobility for phase x, while A_nw
%   and A_w are index matrices that determine the upstream mobility.
%
%   If h_diff is evaluated at h^(n+1) instead of h^n we get a semi implicit
%   method.
%
%
% PARAMETERS:
%   state   - Reservoir solution structure containing valid water
%             saturation state.h(:,1) with one value for each cell
%             in the grid.
%
%   G_top   - Grid data structure discretising the top surface of the
%             reservoir model, as defined by function 'topSurfaceGrid'.
%
%   tf      - End point of time integration interval (i.e., final time),
%             measured in units of seconds.
%
%   fluid   - Data structure as defined by function 'initVEFluid'.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%
%   wells, src, bc - Source terms
%
%   verbose       - Whether or not time integration progress should be
%                   reported to the screen.
%                   Default value: verbose = mrstVerbose.
%
%   dt            - Internal timesteps, measured in units of seconds.
%                   Default value = tf. NB: The explicit scheme is only
%                   stable provided that dt satisfies a CFL time step
%                   restriction.
%
%   time_stepping - Either use a standard CFL condition ('simple'), Coats
%                   formulae ('coats'), or a heuristic bound that allows
%                   for quite optimistic time steps ('dynamic').
%                   Default value: 'simple'
%
%   heightWarn    - Tolerance level for saturation warning.
%                   Default value: satWarn = sqrt(eps).
%
%   computeDt     - Whether or not to compute timestep from CLF condition.
%                   Default value: computeDt = true.
%
%   intVert       - Whether or not integrate permeability from 0 to h. If
%                   false, use average permeability value of column
%                   (rock2D). Default value: intVert = true.
%
%   intVert_poro  - Whether or not to compute pore volume using 3D model
%                   instead of average in z-dir.
%                   Default value: intVert_poro = false.
%
%   preComp       - Use precomputed values calculated in initTransportVE to
%                   speed up computation. Default value: preComp = [].
%
% RETURNS:
%   state  - Reservoir solution with updated saturations, state.h.
%
%   report - Structure reporting timesteps: report.dt
%
% SEE ALSO:
%   `twophaseUpwFE`, `initTransport`, `explicitTransport`, `implicitTransport`.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

opt = struct('dt', tf, 'computedt', true, ...
             'verbose', mrstVerbose,...
             'intVert', false, 'intVert_poro', false, ...
             'wells', [], 'src', [], 'bc', [], ...
             'preComp', [],'heightWarn', sqrt(eps), ...
             'time_stepping','simple','cfl_fac',0.5);

opt = merge_options(opt, varargin{:});
dt  = opt.dt; dt  = min(dt, tf);
assert (dt > 0  && tf > 0);

if ~isempty(opt.wells)
   opt.time_stepping = 'simple';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INIT VARIABLES BEFORE LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
first_time = true; % indication for dynamic time step
is_int   = all(double(G_top.faces.neighbors) > 0, 2);
cells    = G_top.faces.neighbors(is_int,:);
rho_diff = (fluid.rho(1) - fluid.rho(2));

if ~ isfield(state, 'h_max'), state.h_max = state.h; end

%%% Assign precomputed values to variables
if isempty(opt.preComp) %
   opt.preComp = initTransportVE(G_top, averageRock(rock, G_top));
end
grav    = opt.preComp.grav;   g_vec  = opt.preComp.g_vec;
z_diff  = opt.preComp.z_diff; flux   = opt.preComp.flux(state);
pv      = opt.preComp.pv;     K_face = opt.preComp.K_face;
n_z     = opt.preComp.n_z;
opt     = rmfield(opt, 'preComp');
if(strcmp(opt.time_stepping,'coats'))
   K_face_tmp = K_face;
end
%%% Compute sources
q = sources(G_top, state, is_int, 'wells', opt.wells, 'src', opt.src, ...
            'bc', opt.bc);

%%% Compute time step satisfying CFL condition
if (strcmp(opt.time_stepping,'simple') && opt.computedt)
   dt = compute_dt(dt, flux, grav, G_top, fluid, pv, z_diff, ...
                   K_face,q, rock, opt);
end

%%% Precompute values for integration of permeability and inversion of poro
if opt.intVert
   % Find relperm for column (when h = H) by vertical integration
   kr_H   = integrateVertically(rock.perm(:,1), inf, G_top);
   K_face = ones(numel(K_face),1);
end
if opt.intVert_poro
   % precompute values used for inverting porosity
   pv_3D(G_top.columns.cells)=rock.poro(G_top.columns.cells)...
      .*rldecode(G_top.cells.volumes,diff(G_top.cells.columnPos));

   cum_vol3D = accumulateVertically(pv_3D(G_top.columns.cells)'...
               .*G_top.columns.dz,G_top.cells.columnPos);
end

%%% Make function handle to compute pv h^n+1 =pv h^n - H(f_w, dz, dt))
H = @(f,dz) (dz - max(q,0) - min(q,0).*f);

if strcmp(opt.time_stepping,'simple') && opt.verbose
  disp(['Solving transport using ', num2str(ceil(tf/dt)), ...
        ' time steps of size ', num2str(dt/day), ' [day]']);
end

report.timesteps = ceil(tf/dt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN TRANSPORT LOOP ------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0;

while t < tf,
   %%% Compute diffusion/ "capillary pressure"-gradient term
   % compute explicit discretisation of grad h = h_diff:
   h_diff = n_z(cells(:,1)).*fluid.pc(state.h(cells(:,1)))- ...
            n_z(cells(:,2)).*fluid.pc(state.h(cells(:,2)));

   %%% Compute mobilites/flux function
   % Compute cell mobility and fractional flow.
   if opt.intVert
      [kr,dkr]  = integrateVertically(rock.perm(:,1), state.h, G_top);

      if any(state.h_max-state.h > 0) %account for residual CO2
         kr_h_max = integrateVertically(rock.perm(:,1), state.h_max, G_top);
         mob = [fluid.kwm(1)*kr./fluid.mu(1), ...
            fluid.kwm(2)*(kr_H-kr-(kr_h_max-kr)*fluid.res_gas)./fluid.mu(2)];
         dmob = [fluid.kwm(1)*dkr./fluid.mu(1), ...
            fluid.kwm(2)*(dkr*fluid.res_gas)./fluid.mu(2)];
      else
         mob = [fluid.kwm(1)*kr./fluid.mu(1), ...
                fluid.kwm(2)*(kr_H-kr)./fluid.mu(2)];
         dmob = [fluid.kwm(1)*dkr./fluid.mu(1), ...
                fluid.kwm(2)*dkr./fluid.mu(2)];
      end
      f_w = mob(:,1)./sum(mob, 2);
   else
      mob = fluid.mob(state);
      dmob=1./repmat(fluid.mu,G_top.cells.num,1);
      f_w = mob(:,1)./sum(mob, 2);
   end
   % Compute face mobility (Mobility weighting)
   [faceMob,dfaceMob] = computeMobility(G_top, state.flux(is_int), ...
                           -g_vec.*(z_diff+h_diff).*K_face*rho_diff , mob, dmob);

   fw_face = faceMob(:,1) ./ sum(faceMob,2);
   fw_face(sum(faceMob,2)==0) = 0; %Remove possible NaNs

   %%% Build/solve system
   df =  rho_diff*K_face.*faceMob(:,2) .* fw_face;

   %------------------------------------------------------------------------
   %%% Explicit
   % discretise the full equation explicitly.
   %------------------------------------------------------------------------
   dz = flux*fw_face - grav*(df.*(z_diff+h_diff));

   % Use Coats time stepping
   if(strcmp(opt.time_stepping,'coats'))
      % since on faces dont use flux and grav but flux and g_vec
      cellNo   = rldecode(1 : G_top.cells.num, diff(G_top.cells.facePos), 2) .';
      int_flux=state.flux(is_int);
      v_co2   =  (int_flux.*fw_face-g_vec.*(rho_diff*K_face.*fw_face.*faceMob(:,2).*(z_diff+h_diff)));
      v_water =  (int_flux.*(1-fw_face)+g_vec.*(rho_diff*K_face.*(1-fw_face).*faceMob(:,1).*(z_diff+h_diff)));
      tmob=sum(faceMob,2);
      %dfaceMob=1./repmat(fluid.mu,sum(is_int),1);
      %the pluss sign compared to Coats is due to definition of derivative
      ff= faceMob(:,2).*dfaceMob(:,1).*abs(v_co2)./(tmob.*faceMob(:,1))+...
         faceMob(:,1).*dfaceMob(:,2).*abs(v_water)./(tmob.*faceMob(:,2))+...
         -(g_vec.*(rho_diff*K_face.*faceMob(:,2).*fw_face));
      hf_int=is_int(G_top.cells.faces(:,1));
      ff_all=zeros(G_top.faces.num,1);
      ff_all(is_int)=ff;
      dt_tmp=pv(cellNo(hf_int))./(ff_all(G_top.cells.faces(hf_int,1)));
      % coats claims no safety is needed, but experienced oscilation
      dt_tmp=opt.cfl_fac*min(dt_tmp)*(1-(fluid.res_gas+fluid.res_water));

      if dt_tmp > 0
         dt=min(tf-t,dt_tmp);
         fprintf(1,'Dynamic coats time step :\t%.2e [day]\n', ...
                  convertTo(dt, day()));
      else
          dt = compute_dt(dt, flux, grav, G_top, fluid, pv, z_diff, ...
                          K_face_tmp,q, rock, opt);
          opt.cfl_fac = 1;
          dt=min(tf-t,dt);
          fprintf(1,'Use simple timestepping:\t%.2e [day]\n', ...
                  convertTo(dt, day()));
      end
   elseif(any(strcmp(opt.time_stepping,{'dynamic','simple'})))
      if first_time
         pv_res=pv.*(1-fluid.res_water-fluid.res_gas);
         first_time = false;
         dt_para = estimate_parabolic_timestep(dt, G_top, df.*g_vec, ...
                                               n_z, pv_res, opt);
      end
      if(strcmp(opt.time_stepping,'dynamic'))
         pv_res=pv.*(1-fluid.res_water-fluid.res_gas);% use safe value
         dt_dyn = dynamic_dt_from_dh(G_top, H(f_w, dz)./(pv_res), opt, state);
         dt = opt.cfl_fac*dt_dyn;
         dt = min(dt,2*dt_para);
         dt = min(dt,tf-t);
         if opt.verbose,
            fprintf(1,'Dynamic internal time step :\t%.2e [day]\n', ...
                  convertTo(dt, day()));
         end
      elseif(strcmp(opt.time_stepping,'simple'))
         dt=min(dt,dt_para);
         dt=min(dt,tf-t);
      else
         error(['Time stepping ',opt.time_stepping,' not implemented']);
      end
   end

   %%% Hysteresis: account for residual saturation in pore volume
   %              after updating the new total volume of CO2 in the cell
   if(opt.intVert_poro)
      vol= integrateVertically(pv_3D', state.h(:,1), G_top);
      vol_t_max= integrateVertically(pv_3D', state.h_max(:,1), G_top);
      % total volume of CO2 = free CO2 + residual CO2
      vol_t = (vol_t_max-vol)*fluid.res_gas+vol*(1-fluid.res_water);

      % calculate new total volume of CO2
      vol_t = vol_t-dt.*H(f_w, dz);

      % find cells where the total volume of CO2 now is larger than
      % vol_t_max from the previous time step, i.e. where h_max must be
      % updated
      hyst_ind=vol_t>vol_t_max*(1-fluid.res_water);
      % update max total volume in these cells
      vol_t_max(hyst_ind)=vol_t(hyst_ind)./(1-fluid.res_water);

      % calculate volume of free CO2
      vol_new=(vol_t-vol_t_max*fluid.res_gas)./(1-fluid.res_water-fluid.res_gas);
      % update h = height/thickness of free C02
      state.h = invertVerticalFunction(vol_new,G_top,G_top.columns.z,cum_vol3D);
      state.h_max=min(max(state.h_max,state.h), G_top.cells.H);

   else
      % total volume of CO2
      vol=pv.*((state.h_max-state.h)*fluid.res_gas+state.h*(1-fluid.res_water));
      % new total volume of C02
      vol_new=vol-dt.*H(f_w, dz);

      vol_max=state.h_max.*pv*(1-fluid.res_water);
      % cells where h_max has increased during the timestep
      hyst_ind=vol_new>vol_max;
      state.h_max(hyst_ind,:) = vol_new(hyst_ind,:)./...
                                (pv(hyst_ind).*(1-fluid.res_water));
      state.h_max = min(state.h_max, G_top.cells.H);
      % update height/thickness of free CO2
      state.h = (vol_new-state.h_max.*pv.*fluid.res_gas)./ ...
                 (pv.*(1-fluid.res_water-fluid.res_gas));
   end

   % Correct possible bad heights
   state = correct_heights(state, G_top, opt);
   % Save maximum height for use in modeling of relative
   % permeability hysteresis
   state.h_max = max(state.h(:,1), state.h_max);

   t  = t + dt;
   dt = min(dt, tf - t);
end

if size(state.h,2) > 1,
   % Update water/brine height:
   state.h(:,2) = G_top.cells.H - state.h(:,1);
end % end transport loop
if(any(strcmp(opt.time_stepping,{'dynamic'})))
   report.dt_end=dt;
   report.dt_dyn=dt_dyn;
   report.dt_para=dt_para;
end
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function dt = compute_dt(dt, flux, gravmat, G_top, fluid, pv, z_diff, ...
                         K_face, q, rock, opt)
% Estimate time step that satisfies CFL condition. dt is estimated from the
% segregation part (gravity) and the advection part (darcy flux), and we
% take dt to be the minimum of these two multiplied by a factor (0.4).

if(opt.intVert_poro) % compute new pv, pv = min(pv) in each column
   % make array indexing which 2D cell a 3D cell belongs to.
   cell_2d = rldecode((1:G_top.cells.num)',diff(G_top.cells.columnPos));
   % find index to 3D cells
   ind = mcolon(G_top.cells.columnPos(1:end-1), ...
                G_top.cells.columnPos(2:end)-1);
   cell_3d = G_top.columns.cells(ind);
   % find minimal pv in each column (use cell with min(poro) in each column)
   pv =  accumarray(cell_2d,rock.poro(cell_3d).* ...
                    G_top.cells.volumes(cell_2d),[],@min);
end

   gravmat = gravmat*(fluid.rho(1)-fluid.rho(2));
   gravmat = gravmat*sparse(1:numel(z_diff),1:numel(z_diff),K_face.*z_diff);

   h  = linspace(0, max(G_top.cells.H),  1001) .';

   deriv = @(f) diff(f(fluid.mob_avg(struct('h',h), ...
                         max(G_top.cells.H)))) ./ diff(h);

   if nnz(gravmat),
      f  = @(mob) (mob(:,1).*mob(:,2)./sum(mob,2));

      dt_grav = estimate_dt(gravmat, deriv(f), pv, 0);
      if opt.verbose,
         fprintf('\nEstimated segregation step:\t%.2e [day]\n', ...
                  convertTo(dt_grav, day()));
      end
      dt = min(dt, dt_grav);
   end
   if nnz(flux)
      f  = @(mob) (mob(:,1) ./ sum(mob, 2));

      dt_advection = estimate_dt(flux, deriv(f), pv, q);
      if opt.verbose,
         fprintf('Estimated advection step:\t%.2e [day]\n', ...
                  convertTo(dt_advection, day()));
      end
      dt = min(dt, dt_advection);
   end

   dt = dt * opt.cfl_fac*(1-(fluid.res_gas+fluid.res_water));
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

%--------------------------------------------------------------------------
function dt_para = estimate_parabolic_timestep(dt, G_top, vals_face, n_z, pv_res, opt) %#ok
   if any(vals_face)
      % compute an estimate of the parabolic time step (grad pc)
      is_int = all(double(G_top.faces.neighbors) > 0, 2);
      nc = G_top.cells.num;
      % build semi implicit matrix:
      A = sparse(double(G_top.faces.neighbors(is_int,1)), ...
         double(G_top.faces.neighbors(is_int,2)), vals_face, nc, nc);
      A = A + A';
      %A = (A*spdiags(n_z, 0, nc, nc) - spdiags(sum(A,2).*n_z, 0, nc, nc));
      A = (A - spdiags(sum(A,2), 0, nc, nc));
      % A = A
      try
         d = eigs(A,1,'LA',struct('disp',0,'maxit',100));
         dt_para = 1./(d/min(pv_res));
      catch %#ok
         warning('Did not find parabolic time step restriction'); %#ok
         dt_para = 1e99;
      end
   else
      dt_para = Inf;
   end
   if opt.verbose,
       fprintf(1, ['Estimated parabolic step: ', ...
                  '\t%.2e [day]\n'], convertTo(dt_para, day()));
   end
   if opt.verbose && (dt_para<dt)
      warning(['Time step will be limited by parabolic term.']); %#ok
   end
end
%--------------------------------------------------------------------------
function dt = dynamic_dt_from_dh(G_top,dh,opt,state)
% estimate dynamic timestep from dh = H(f_w, dz)./pv_res

internal = sum(G_top.faces.neighbors==0,2)==0;
cellNo   = rldecode(1 : G_top.cells.num, diff(G_top.cells.facePos), 2) .';

dhf = diff(state.h(G_top.faces.neighbors(internal,:)),1,2);
%nb  = reshape(G_top.faces.neighbors(internal,:), [],1);
%val = reshape(repmat(dhf,1,2), [],1);
%dhmax=accumarray(nb,abs(val),[G_top.cells.num 1],@max);

dhff = zeros(G_top.faces.num,1);
dhff(internal) = dhf;
%dhmax=accumarray(cellNo,dhff(G_top.cells.faces(:,1)),[],@max);%toc
dhmax = accumarray(cellNo,abs(dhff(G_top.cells.faces(:,1))));
dhmax=dhmax/2;
%ind=find(abs(dh)>0 & dhmax>1e-6);
if(~isempty(opt.src))
   dhmax(opt.src.cells)=-(state.h(opt.src.cells)-G_top.cells.H(opt.src.cells));
end
if(~isempty(opt.wells))
   dhmax(opt.wells.cells)=-(state.h(opt.wells.cells)-G_top.cells.H(opt.wells.cells));
end
ind = find(dh>0);
dhmax(ind) = min(abs(dhmax(ind)),state.h(ind));
ind = find(dh<0);
dhmax(ind) = min(abs(dhmax(ind)),G_top.cells.H(ind)-state.h(ind));
ind = find(abs(dh)>0 & dhmax>1e-6);
dt  = min(dhmax(ind)./abs(dh(ind)));
assert(dt>0)
end

%--------------------------------------------------------------------------

function state = correct_heights(state, G_top, opt)
   % Make sure 0 <= state.h <= H for all cells
   inx1 = find(state.h(:,1)-G_top.cells.H > 0);
   inx2 = find(state.h(:) < 0);

   if opt.verbose,
      sMax = max(state.h(:,1)-G_top.cells.H);
      if sMax > opt.heightWarn
         disp('h exceeds H in cells from: [inx h]')
         fprintf(1,'%d %g \n', [inx1, state.h(inx1,1)]');
      end
      sMin = min(state.h(:,1));
      if sMin <  -opt.heightWarn,
         disp('Heigth less than 0 in cells from: [inx h]')
         fprintf(1,'%d %g \n', [inx2, state.h(inx2,1)]');
      end
   end

   state.h(inx1,1) = G_top.cells.H(inx1) - eps;
   state.h(inx2) = 0;
end

%--------------------------------------------------------------------------

function q = sources(G_top, state, is_int,  varargin)
% Contributions from wells, sources and boundary conditions
   opt = struct('wells'    , [],          ...
                'src'      , [],          ...
                'bc'       , []);

   opt = merge_options(opt, varargin{:});
   %-----------------------------------------------------------------------
   % External sources/sinks (e.g. wells and BC's) ------------------------
   %
   qi = [];  % Cells to which sources are connected
   qh = [];  % Actual strength of source term (in m^3/s).

   if ~isempty(opt.wells),
      assert (isfield(state, 'wellSol'));
      [ii, h] = contrib_wells(opt.wells, state.wellSol,G_top);
      qi = [qi; ii];
      qh = [qh; h];
   end

   if ~isempty(opt.src), assert (~isempty(opt.src.h))
      [ii, h] = contrib_src(opt.src,G_top);
      qi = [qi; ii];
      qh = [qh; h];
   end

   if ~isempty(opt.bc), assert (~isempty(opt.bc.h))
      [ii, h] = contrib_bc(G_top, state, opt.bc, is_int);
      qi = [qi; ii];
      qh = [qh; h];
   end

   %-----------------------------------------------------------------------
   % Assemble final phase flux and source contributions in SPARSE format -
   %
   q  = sparse(qi, 1, qh, G_top.cells.num, 1);
end

%--------------------------------------------------------------------------
function [qi, qh] = contrib_wells(W, wellSol,g_top)
   % Contributions from wells as defined by 'addWell'.

   nperf = cellfun(@numel, { W.cells }) .';

   qi = vertcat(W.cells);
   qh = vertcat(wellSol.flux);

   % Injection perforations have positive flux (into reservoir).
   %
   %assert(nperf==1);
   comp      = rldecode(vertcat(W.h), nperf);
   inj_p     = qh > 0;
   qh(inj_p) = qh(inj_p).*comp(inj_p)./g_top.cells.H(qi(inj_p));
end

%--------------------------------------------------------------------------

function [qi, qh] = contrib_src(src,g_top)
   % Contributions from explicit sources defined by (e.g.) 'addSource'.

   qi = src.cell;
   qh = src.rate;

   % Injection sources have positive rate into reservoir.
   %
   in = find(src.rate > 0);
   if ~isempty(in),
      qh(in) = qh(in) .* src.h(in,1)./g_top.cells.H(qi(in));
   end
end

%--------------------------------------------------------------------------

function [qi, qh] = contrib_bc(G, resSol, bc, is_int)
   % Contributions from boundary conditions as defined by 'addBC'.

   qh = zeros([G.faces.num, 1]);
   dF = false([G.faces.num, 1]);

   isDir = strcmp('pressure', bc.type);
   isNeu = strcmp('flux',     bc.type);

   dF(bc.face(isDir))      = true;
   cfIx                    = dF(G.cells.faces(:,1));

   cflux = faceFlux2cellFlux(G, resSol.flux);
   qh(G.cells.faces(cfIx,1)) = -cflux(cfIx);
   qh(bc.face(isNeu))        = bc.value(isNeu);

   % Injection BC's have positive rate (flux) into reservoir.
   %
   is_inj = qh > 0;
   if any(is_inj),
      qh(is_inj) = qh(is_inj) .* bc.h(is_inj(bc.face), 1);
   end

   is_outer = ~is_int;

   qi = sum(G.faces.neighbors(is_outer,:), 2);
   qh = qh(is_outer);
end

function [facemob,dfacemob] = computeMobility(G, v_darcy, g_vec, mob, dmob)
   %
   % Given fixed dflux and gflux, we can compute upwind directions for
   % phase w given that dflux and gflux have the same sign.  If dflux and
   % gflux have opposite signs, then the upwind directions for phase o can be
   % computed.
   %
   % If either dflux or gflux is zero, the phase velocity does not depend on
   % the phase mobilities. Then, the sequence induced by the present function
   % is as good as any other.
   %
   % We store these precomputed values to speed up calculations of upwind
   % direcions in findFaceMobIx.

   intern     = all(G.faces.neighbors ~= 0, 2);

   a = ~(v_darcy<0) & ~(g_vec<0);
   b = ~(v_darcy>0) & ~(g_vec>0);
   vec = a | b;
   nvec = find(~vec);
   vec = find(vec);

   iw = zeros(sum(intern), 1);
   io = iw;

   ineighbors = G.faces.neighbors(intern,:);%ineighbors(b, [2,1]);
   ineighbors(b, :) = ineighbors(b,[2,1]);
   iw( vec) = ineighbors( vec, 1);

   ineighbors = G.faces.neighbors(intern,:);%ineighbors(b, [2,1]);
   b = ~(v_darcy>0) & ~(g_vec<0);
   ineighbors(b, :) = ineighbors(b,[2,1]);
   io(nvec) = ineighbors(nvec, 1);

   mw = zeros(size(iw));
   dmw = zeros(size(iw));
   mo = mw;
   dmo = dmw;

   mw( vec) = mob(iw( vec), 1);
   mo(nvec) = mob(io(nvec), 2);
   dmw( vec) = dmob(iw( vec), 1);
   dmo(nvec) = dmob(io(nvec), 2);

   vw = v_darcy +  g_vec.*mo;
   vo = v_darcy -  g_vec.*mw;

   ineighbors = G.faces.neighbors(intern,:);
   ineighbors(vw<0, :) = ineighbors(vw<0, [2,1]);
   iw(nvec) = ineighbors(nvec, 1);

   ineighbors = G.faces.neighbors(intern,:);
   ineighbors(vo<0, :) = ineighbors(vo<0, [2,1]);
   io(vec) = ineighbors(vec, 1);

   mw(nvec) = mob(iw(nvec),1);
   mo(vec)  = mob(io(vec),2);
   dmw(nvec) = dmob(iw(nvec),1);
   dmo(vec)  = dmob(io(vec),2);
   facemob= [mw,mo];
   dfacemob= [dmw,dmo];
end

