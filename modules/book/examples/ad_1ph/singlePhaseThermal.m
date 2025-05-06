%% Single-phase compressible AD solver with thermal effects
% The purpose of the example is to give demonstrate rapid prototying in
% MRST using the automatic differentiation (AD) class. To this end, we
% extend the <singlePhaseAD.m flow simulator for compressible single-phase
% flow> to include temperature effects, which are modeled by incorporating
% a second conservation equation for energy. Except for this, the
% computational setup is the same with a single horizontal well draining
% fluids from a simple box-geometry reservoir.
mrstModule add ad-core

%% Set up model: grid, perm and poro
[nx,ny,nz] = deal( 10,  10, 10);
[Dx,Dy,Dz] = deal(200, 200, 50);
G = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
G = computeGeometry(G);

%% Define rock model
rock = makeRock(G, 30*milli*darcy, 0.3);

cr   = 1e-6/barsa;
p_r  = 200*barsa;
pv_r = poreVolume(G, rock);
pv   = @(p) pv_r .* exp( cr * (p - p_r) );
sv   = @(p) G.cells.volumes - pv(p);

%% Define model for compressible fluid
mu0   = 5*centi*poise;
cmup  = 2e-3/barsa;
cmut  = 1e-3;
T_r   = 300;
mu    = @(p,T)  mu0*(1+cmup*(p-p_r)).*exp(-cmut*(T-T_r));

beta  = 1e-3/barsa;
alpha = 5e-3;
rho_r = 850*kilogram/meter^3;
rho_S = 750*kilogram/meter^3;
rho   = @(p,T) rho_r .* (1+beta*(p-p_r)).*exp(-alpha*(T-T_r) );

%% Quantities for energy equation
Cp = 4e3;
Cr = 2*Cp;
Hf = @(p,T) Cp*T+(1-T_r*alpha).*(p-p_r)./rho(p,T);
Ef = @(p,T) Hf(p,T) - p./rho(p,T);
Er = @(T)   Cp*T;

%% Assume a single horizontal well
nperf = floor(G.cartDims(2)*(ny-2)/ny);
I = repmat(2, [nperf, 1]);
J = (1 : nperf).' + 1;
K = repmat(nz/2, [nperf, 1]);
cellInx = sub2ind(G.cartDims, I, J, K);
W = addWell([ ], G, rock, cellInx, 'Name', 'P1', 'Dir', 'y');

%% Impose vertical equilibrium
gravity reset on; g = norm(gravity);
[z_0, z_max] = deal(0, max(G.cells.centroids(:,3)));
equil  = ode23(@(z,p) g.* rho(p,T_r), [z_0, z_max], p_r);
p_init = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil
T_init = ones(G.cells.num,1)*T_r;

%% Plot well and initial pressure
figure
show = true(G.cells.num,1);
cellInx = sub2ind(G.cartDims, ...
   [I-1; I-1; I; I;   I(1:2)-1], ...
   [J  ; J;   J; J;   nperf+[2;2]], ...
   [K-1; K;   K; K-1; K(1:2)-[0; 1]]);
show(cellInx) = false;
plotCellData(G,p_init/barsa, show,'EdgeColor','k');
plotWell(G,W, 'height',10);
view(-125,20), camproj perspective

%% Compute transmissibilities
N   = double(G.faces.neighbors);
intInx = all(N ~= 0, 2);
cf  = G.cells.faces(:,1);
nf  = G.faces.num;
N   = N(intInx, :);                           % Interior neighbors

hT  = computeTrans(G, rock);                  % Half-transmissibilities
Tp  = 1 ./ accumarray(cf, 1 ./ hT, [nf, 1]);  % Harmonic average
Tp  = Tp(intInx);                             % Restricted to interior

kap = 4;                                      % Heat conduction of granite
tmp = struct('perm',kap*ones(G.cells.num,1)); % Temporary rock object
hT  = computeTrans(G, tmp);                   % Half-transmissibilities
Th  = 1 ./ accumarray(cf, 1 ./ hT, [nf, 1]);  % Harmonic average
Th  = Th(intInx);                             % Restricted to interior

%% Define discrete operators
n    = size(N,1);
C    = sparse( [(1:n)'; (1:n)'], N, ones(n,1)*[-1 1], n, G.cells.num);
grad = @(x) C*x;
div  = @(x) -C'*x;
avg  = @(x) 0.5 * (x(N(:,1)) + x(N(:,2)));
upw  = @(x,flag) x(N(:,1)).*double(flag)+x(N(:,2)).*double(~flag);

%% Define flow equations
% Writing in functional form means that v(p,T) is evaluated two
% times more than strictly needed if the whole definition of
% discrete equations written as a single function
gdz = grad(G.cells.centroids)*gravity()';
v   = @(p,T)  -(Tp./mu(avg(p),avg(T))).*(grad(p) - avg(rho(p,T)).*gdz);
pEq = @(p,T, p0,T0, dt) ...
     (1/dt)*(pv(p).*rho(p,T) - pv(p0).*rho(p0,T0)) ...
      + div( avg(rho(p,T)).*v(p,T) );
hEq = @(p, T, p0, T0, dt) ...
     (1/dt)*(pv(p ).*rho(p, T ).*Ef(p ,T ) + sv(p ).*Er(T ) ...
           - pv(p0).*rho(p0,T0).*Ef(p0,T0) - sv(p0).*Er(T0)) ...
     + div( upw(Hf(p,T),v(p,T)>0).*avg(rho(p,T)).*v(p,T) ) ...
     + div( -Th.*grad(T));

%% Define well equations.
wc = W(1).cells; % connection grid cells
assert(numel(wc)==numel(unique(wc))); % each cell should only appear once
WI  = W(1).WI;    % well-indices
dz  = W(1).dZ;    % connection depth relative to bottom-hole
bhT = ones(size(wc))*200; % temperature of wells not used if production

p_conn  = @(bhp,bhT) bhp + g*dz.*rho(bhp,bhT); %connection pressures
q_conn  = @(p,T, bhp) ...
    WI .* (rho(p(wc),T(wc))./mu(p(wc),T(wc))) .* (p_conn(bhp,bhT)-p(wc));

rateEq = @(p,T, bhp, qS)  qS-sum(q_conn(p, T, bhp))/rho_S;
ctrlEq = @(bhp)           bhp-100*barsa;

%% Initialize for solution loop
[p_ad, T_ad, bhp_ad, qS_ad] = ...
    initVariablesADI(p_init,T_init, p_init(wc(1)), 0);
nc = G.cells.num;
[pIx, TIx, bhpIx, qSIx] = deal(1:nc, nc+1:2*nc, 2*nc+1, 2*nc+2);

numSteps = 78;                  % number of time-steps
totTime  = 365*day*1.5;         % total simulation time
dt       = totTime / numSteps;  % constant time step
tol      = 1e-5;                % Newton tolerance
maxits   = 10;                  % max number of Newton its

sol = repmat(struct('time',[], 'pressure',[], 'bhp',[], 'qS',[], ...
    'T',[], 'qH',[]),[numSteps+1,1]);
sol(1)  = struct('time', 0, 'pressure', value(p_ad), ...
    'bhp', value(bhp_ad), 'qS', value(qS_ad), 'T', value(T_ad),'qH', 0);

%% Main loop
t = 0; step = 0;
hwb = waitbar(0,'Simulation..');
while t < totTime
   t = t + dt;
   step = step + 1;
   fprintf('\nTime step %d: Time %.2f -> %.2f days\n', ...
      step, convertTo(t - dt, day), convertTo(t, day));

   % Newton loop
   resNorm = 1e99;
   p0  = value(p_ad); % Previous step pressure
   T0  = value(T_ad);
   nit = 0;
   while (resNorm > tol) && (nit < maxits)

      % Add source terms to homogeneous pressure equation:
      eq1      = pEq(p_ad,T_ad, p0, T0,dt);
      qw       = q_conn(p_ad,T_ad, bhp_ad);
      eq1(wc)  = eq1(wc) - qw;
      hq       = Hf(bhp_ad,bhT).*qw;    %inflow not in this example
      Hcells   = Hf(p_ad,T_ad);
      hq(qw<0) = Hcells(wc(qw<0)).*qw(qw<0);
      eq2      = hEq(p_ad,T_ad, p0, T0,dt);
      eq2(wc)  = eq2(wc) - hq;

      % Collect all equations. Scale residual of energy equation
      eqs = {eq1, eq2/Cp, rateEq(p_ad,T_ad, bhp_ad, qS_ad), ctrlEq(bhp_ad)};

      % Concatenate equations and solve for update:
      eq  = cat(eqs{:});
      J   = eq.jac{1};  % Jacobian
      res = eq.val;     % residual
      upd = -(J \ res); % Newton update

      % Update variables
      p_ad.val   = p_ad.val   + upd(pIx);
      T_ad.val   = T_ad.val   + upd(TIx);
      bhp_ad.val = bhp_ad.val + upd(bhpIx);
      qS_ad.val  = qS_ad.val  + upd(qSIx);

      resNorm = norm(res);
      nit     = nit + 1;
      fprintf('  Iteration %3d:  Res = %.4e\n', nit, resNorm);
   end

   if nit > maxits
      error('Newton solves did not converge')
   else % store solution
      sol(step+1)  = struct('time', t, 'pressure', value(p_ad), ...
                            'bhp', value(bhp_ad), 'qS', value(qS_ad), ...
			    'T', value(T_ad), 'qH', sum(value(hq)));
      waitbar(t/totTime,hwb);
   end
end
close(hwb)
%% Plot production rate
clf, set(gca,'FontSize',20);
stairs([sol(2:end).time]/day,-[sol(2:end).qS]*day,'LineWidth',2);
xlabel('days');ylabel('m^3/day')

%% Plot pressure evolution
figure; clf;
steps = [2 5 10 25];
%steps = floor(1:(numel(sol)/6-eps):numel(sol));
p = vertcat(sol(:).pressure);
cax=[min(p) max(p)]./barsa;
for i=1:4
   subplot(2,2,i);
   set(gca,'Clipping','off');
   plotCellData(G, sol(steps(i)).pressure/barsa, show, 'EdgeColor',.5*[1 1 1]);
   plotWell(G, W, 'FontSize',12);
   view(-125,20), camproj perspective
   caxis(cax);
   axis tight off; zoom(1.1)
   text(200,170,-8,[num2str(round(steps(i)*dt/day)) ' days'],'FontSize',12);
end
colorbar('horiz','Position',[.1 .05 .8 .025]);
colormap(jet(55));

%%
figure(),clf;
steps = [2 5 10 numel(sol)];
T = vertcat(sol(:).T); %-T_r;
cax=[min(T) max(T)];
for i=1:numel(steps)
   subplot(2,2,i);
   set(gca,'Clipping','off');
   plotCellData(G, sol(steps(i)).T, show, 'EdgeColor',.5*[1 1 1]);
   plotWell(G, W, 'FontSize',12);
   view(-125,20), camproj perspective
   caxis(cax);
   axis tight off; zoom(1.1)
   text(200,170,-8,[num2str(round(steps(i)*dt/day)) ' days'],'FontSize',12);
end
h=colorbar('horiz','Position',[.1 .05 .8 .025]);
colormap(jet(55));
%%
ns = numel(sol);
Tw = nan(ns,numel(wc));
Pw = nan(ns,numel(wc));
[t,Pm,PM,Pa,Tm,TM,Ta] = deal(nan(ns,1));
for i=1:numel(sol)
   t(i) = sol(i).time/day;
   Tw(i,:)=sol(i).T(wc);
   Tm(i) = min(sol(i).T);
   TM(i) = max(sol(i).T);
   Ta(i) = mean(sol(i).T);
   Pw(i,:)=sol(i).pressure(wc)./barsa;
   Pm(i) = min(sol(i).pressure./barsa);
   PM(i) = max(sol(i).pressure./barsa);
   Pa(i) = mean(sol(i).pressure./barsa);
end
figure; plot(t,Pm,t,Pa,t,PM,t,Pw,'.k','LineWidth',2);
legend('min(p)', 'avg(p)', 'max(p)', 'wells');
figure; plot(t,Tm,t,Ta,t,TM,t,Tw,'.k','LineWidth',2);
legend('min(T)', 'avg(T)', 'max(T)', 'wells');

%% Compute the three different expansion temperatures
[p,T] = initVariablesADI(p_r,T_r);
dp   = Pm(end)*barsa-p_r;

% Joule-Thomson
hf   = Hf(p,T);
dHdp = hf.jac{1};
dHdT = hf.jac{2};
Tjt  = T_r - dHdp*dp/dHdT;
hold on, plot(t([1 end]), [Tjt Tjt],'--k'); hold off
text(t(5), Tjt+.25, 'Joule-Tompson');

% linearized adiabatic temperature
hf   = Ef(p,T) + value(p)./rho(p,T);
dHdp = hf.jac{1};
dHdT = hf.jac{2};
Tab  = T_r - dHdp*dp/dHdT;
hold on, plot(t([1 end]), [Tab Tab],'--k'); hold off
text(t(2), Tab-.35, 'Adiabatic expansion');

% free expansion
hf = Ef(p,T);
dHdp = hf.jac{1};
dHdT = hf.jac{2};
Tfr  = T_r - dHdp*dp/dHdT;
hold on, plot(t([1 end]), [Tfr Tfr],'--k'); hold off
text(t(end/2), Tfr+.25, 'Free expansion');
set(gca,'XLim',t([1 end]));

%%
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
