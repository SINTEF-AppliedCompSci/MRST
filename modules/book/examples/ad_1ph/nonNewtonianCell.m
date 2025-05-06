%% Non-Newtonian fluid
% In this example, we will demonstrate how one can easily extend the
% compressible single-phase pressure solver to include the effect of
% non-Newtonian fluids modelled using a simple power law in which the
% effective viscosity depends on the norm of the velocity

%% Define geometric quantitites
% Grid that represents the reservoir geometry
[nx,ny,nz] = deal( 10,  10, 10);
[Lx,Ly,Lz] = deal(200, 200, 50);
G = cartGrid([nx, ny, nz], [Lx, Ly, Lz]);
G = computeGeometry(G);

% Discrete operators
N  = double(G.faces.neighbors);
intInx = all(N ~= 0, 2);
N  = N(intInx, :);
n = size(N,1);
C = sparse( [(1:n)'; (1:n)'], N, ones(n,1)*[-1 1], n, G.cells.num);
aC = bsxfun(@rdivide,0.5*abs(C),G.faces.areas(intInx))';
grad = @(x) C*x;
div  = @(x) -C'*x;
cavg = @(x) aC*x;
favg = @(x) 0.5 * (x(N(:,1)) + x(N(:,2)));
clear aC C N;

%% Rock model and transmissibilities
rock = makeRock(G, 30*milli*darcy, 0.3);

cr   = 1e-6/barsa;
p_r  = 200*barsa;
pv_r = poreVolume(G, rock);
pv   = @(p) pv_r .* exp( cr * (p - p_r) );
clear pv_r;

hT = computeTrans(G, rock);
cf = G.cells.faces(:,1);
nf = G.faces.num;
T  = 1 ./ accumarray(cf, 1 ./ hT, [nf, 1]);
T  = T(intInx);
clear hT;

%% Basic fluid model
c     = 1e-3/barsa;
rho_r = 850*kilogram/meter^3;
rhoS  = 750*kilogram/meter^3;
rho   = @(p) rho_r .* exp( c * (p - p_r) );
if exist('fluidModel', 'var') 
   mu0 = fluidModel.mu0;
   nmu = fluidModel.nmu;
   Kc  = fluidModel.Kc;
   Kbc = (Kc/mu0)^(2/(nmu-1))*36*((3*nmu+1)/(4*nmu))^(2*nmu/(nmu-1));
   if nmu==1, Kbc = 0; end
else
   mu0 = 100*centi*poise;
   nmu = 0.25;
   Kc  = .1;
   Kbc = (Kc/mu0)^(2/(nmu-1))*36*((3*nmu+1)/(4*nmu))^(2*nmu/(nmu-1));
end

%% Initial vertical equilibrium
gravity reset on, g = norm(gravity);
[z_0, z_max] = deal(0, max(G.cells.centroids(:,3)));
equil  = ode23(@(z,p) g .* rho(p), [z_0, z_max], p_r);
p_init = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);
clear equil z_0 z_max;

%% Constant for the simulation
numSteps = 52;
totTime  = 365*day;
dt       = totTime / numSteps;
tol      = 1e-5;
maxits   = 100;

%% Flow equations
phiK  = rock.perm.*rock.poro;
gradz = grad(G.cells.centroids(:,3));
v     = @(p, eta) ...
        -(T./(mu0*favg(eta))).*( grad(p) - g*favg(rho(p)).*gradz );
etaEq = @(p, eta) ...
        eta - ( 1 + Kbc* cavg(v(p,eta)).^2 ./phiK ).^((nmu-1)/2);
presEq= @(p, p0, eta, dt)  ...
        (1/dt)*(pv(p).*rho(p) - pv(p0).*rho(p0)) + div(favg(rho(p)).*v(p, eta));

%% Well model
nperf = 8;
I = repmat(2, [nperf, 1]);
J = (1 : nperf).' + 1;
K = repmat(5, [nperf, 1]);
cellInx = sub2ind(G.cartDims, I, J, K);
W = addWell([ ], G, rock, cellInx, 'Name', 'P1', 'Dir', 'x' );

% Define well equations
wc = W(1).cells; % connection grid cells
WI = W(1).WI;    % well-indices
dz = W(1).dZ;    % connection depth relative to bottom-hole

p_conn = @(bhp) ...
   bhp + g*dz.*rho(bhp);
q_conn = @(p, eta, bhp) ...
   WI .* (rho(p(wc)) ./ (mu0*eta(wc))) .* (p_conn(bhp) - p(wc));
rateEq = @(p, eta, bhp, qS) ...
   qS - sum(q_conn(p, eta, bhp))/rhoS;
ctrlEq = @(bhp) ...
   bhp - 300*barsa;

%% Initialize for solution loop
nc = G.cells.num;
[p_ad, eta_ad, bhp_ad, qS_ad] = ...
   initVariablesADI(p_init, ones(nc,1), p_init(wc(1)), 0);
[pIx, etaIx, bhpIx, qSIx] = ...
   deal(1:nc, nc+1:2*nc, 2*nc+1, 2*nc+2);
sol = repmat(struct('time',[],'pressure',[],'eta',[], ...
                    'bhp',[],'qS',[]), [numSteps+1,1]);
sol(1) = struct('time', 0, 'pressure', value(p_ad), ...
                'eta', value(eta_ad), ...
                'bhp', value(bhp_ad), 'qS', value(qS_ad));
[etamin, etawmin, etamean] = deal(zeros(numSteps,1));

%% Time loop
t = 0; step = 0;
while t < totTime

   % Increment time
   t = t + dt;
   step = step + 1;
   fprintf('Time step %d: Time %.2f -> %.2f days\n', ...
      step, convertTo(t - dt, day), convertTo(t, day));

   % Main Newton loop
   p0  = value(p_ad); % Previous step pressure
   [resNorm,nit] = deal(1e99, 0);
   while (resNorm > tol) && (nit < maxits)

      % Newton loop for eta (effective viscosity)
      [resNorm2,nit2] = deal(1e99, 0);
      eta_ad2 = initVariablesADI(eta_ad.val);
      while (resNorm2 > tol) && (nit2 <= maxits)
         eeq = etaEq(p_ad.val, eta_ad2);
         res = eeq.val;
         eta_ad2.val = eta_ad2.val - (eeq.jac{1} \ res);

         resNorm2 = norm(res);
         nit2     = nit2+1;
      end
      if nit2 > maxits
         error('Local Newton solves did not converge')
      else
         eta_ad.val = eta_ad2.val;
      end

      % Add source terms to homogeneous pressure equation:
      eq1     = presEq(p_ad, p0, eta_ad, dt);
      eq1(wc) = eq1(wc) - q_conn(p_ad, eta_ad, bhp_ad);

      % Collect all equations
      eqs = {eq1, etaEq(p_ad, eta_ad), ...
         rateEq(p_ad, eta_ad, bhp_ad, qS_ad), ctrlEq(bhp_ad)};

      % Concatenate equations and solve for update:
      eq  = cat(eqs{:});
      J   = eq.jac{1};  % Jacobian
      res = eq.val;     % residual
      upd = -(J \ res); % Newton update
      % Update variables
      p_ad.val   = p_ad.val   + upd(pIx);
      eta_ad.val = eta_ad.val + upd(etaIx);
      bhp_ad.val = bhp_ad.val + upd(bhpIx);
      qS_ad.val  = qS_ad.val  + upd(qSIx);

      resNorm = norm(res);
      nit     = nit + 1;
   end

 %  clf,
 %  plotCellData(G,eta_ad.val,'FaceAlpha',.3,'EdgeAlpha', .1);
 %  view(3); colorbar; drawnow

   if nit > maxits
      error('Newton solves did not converge')
   else % store solution
      sol(step+1)  = struct('time', t, 'pressure', value(p_ad), ...
         'eta', value(eta_ad), ...
         'bhp', value(bhp_ad), 'qS', value(qS_ad));
   end
   etamin (step) = min(eta_ad.val);
   etawmin(step) = min(eta_ad.val(wc));
   etamean(step) = mean(eta_ad.val);
      
end

%%
clf
[ha,hr,hp] = ...
   plotyy([sol(2:end).time]/day, [sol(2:end).qS]*day, ...
          [sol(2:end).time]/day, mean([sol(2:end).pressure]/barsa), ...
          'stairs', 'plot');
%set(ha,'FontSize',16);
set(hr,'LineWidth', 2);
set(hp,'LineStyle','none','Marker','o','LineWidth', 1);
xlabel('time [days]');
ylabel(ha(1), 'rate [m^3/day]');
p=get(gca,'Position'); p(1)=p(1)-.01; p(2)=p(2)+.02; set(gca,'Position',p);
ylabel(ha(2), 'avg pressure [bar]');

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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
