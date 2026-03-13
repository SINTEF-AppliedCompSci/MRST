%% Single-phase compressible AD solver
% The purpose of the example is to give the first introduction to how one
% can use the automatic differentiation (AD) class in MRST to write a flow
% simulator for a compressible single-phase model. For simplicity, the
% reservoir is assumed to be a rectangular box with homogeneous properties
% and no-flow boundaries. Starting from a hydrostatic initial state, the
% reservoir is produced from a horizontal well that will create a zone of
% pressure draw-down. As more fluids are produced, the average pressure in
% the reservoir drops, causing a gradual decay in the production rate.

%% Set up model geometry
[nx,ny,nz] = deal( 10,  10, 10);
[Dx,Dy,Dz] = deal(200, 200, 50);
G = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
G = computeGeometry(G);

plotGrid(G); view(3); axis tight

%% Define rock model
rock = makeRock(G, 30*milli*darcy, 0.3);

cr   = 1e-6/barsa;
p_r  = 200*barsa;
pv_r = poreVolume(G, rock);
pv   = @(p) pv_r .* exp( cr * (p - p_r) );

p = linspace(100*barsa,220*barsa,50);
plot(p/barsa, pv_r(1).*exp(cr*(p-p_r)),'LineWidth',2);


%% Define model for compressible fluid
mu    = 5*centi*poise;
c     = 1e-3/barsa;
rho_r = 850*kilogram/meter^3;
rhoS  = 750*kilogram/meter^3;
rho   = @(p) rho_r .* exp( c * (p - p_r) );

plot(p/barsa,rho(p),'LineWidth',2);

%% Assume a single horizontal well
nperf = 8;
I = repmat(2, [nperf, 1]);
J = (1 : nperf).' + 1;
K = repmat(5, [nperf, 1]);
cellInx = sub2ind(G.cartDims, I, J, K);
W = addWell([ ], G, rock, cellInx, 'Name', 'P1', 'Dir', 'y' );

%% Impose vertical equilibrium
gravity reset on, g = norm(gravity);
[z_0, z_max] = deal(0, max(G.cells.centroids(:,3)));
equil  = ode23(@(z,p) g .* rho(p), [z_0, z_max], p_r);
p_init = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil

%% Plot well and initial pressure
clf
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
N  = double(G.faces.neighbors);
intInx = all(N ~= 0, 2);
N  = N(intInx, :);                          % Interior neighbors
hT = computeTrans(G, rock);                 % Half-transmissibilities
cf = G.cells.faces(:,1);
nf = G.faces.num;
T  = 1 ./ accumarray(cf, 1 ./ hT, [nf, 1]); % Harmonic average
T  = T(intInx);                             % Restricted to interior

%% Define discrete operators
n = size(N,1);
C = sparse( [(1:n)'; (1:n)'], N, ones(n,1)*[-1 1], n, G.cells.num);
grad = @(x)C*x;
div  = @(x)-C'*x;
avg  = @(x) 0.5 * (x(N(:,1)) + x(N(:,2)));
spy(C)

%% Define flow equations
gradz  = grad(G.cells.centroids(:,3));
v      = @(p)  -(T/mu).*( grad(p) - g*avg(rho(p)).*gradz );
presEq = @(p,p0,dt) (1/dt)*(pv(p).*rho(p) - pv(p0).*rho(p0)) ...
                      + div( avg(rho(p)).*v(p) );

%% Define well equations
wc = W(1).cells; % connection grid cells
WI = W(1).WI;    % well-indices
dz = W(1).dZ;    % connection depth relative to bottom-hole

p_conn  = @(bhp)  bhp + g*dz.*rho(bhp); %connection pressures
q_conn  = @(p,bhp) WI .* (rho(p(wc)) / mu) .* (p_conn(bhp) - p(wc));

rateEq = @(p,bhp,qS)  qS-sum(q_conn(p, bhp))/rhoS;
ctrlEq = @(bhp)       bhp-100*barsa;

%% Initialize for solution loop
[p_ad, bhp_ad, qS_ad] = initVariablesADI(p_init, p_init(wc(1)), 0);
nc = G.cells.num;
[pIx, bhpIx, qSIx] = deal(1:nc, nc+1, nc+2);

numSteps = 52;                  % number of time-steps
totTime  = 365*day;             % total simulation time
dt       = totTime / numSteps;  % constant time step
tol      = 1e-5;                % Newton tolerance
maxits   = 10;                  % max number of Newton its

sol = repmat(struct('time',[],'pressure',[],'bhp',[],'qS',[]),[numSteps+1,1]);
sol(1)  = struct('time', 0, 'pressure', value(p_ad), ...
   'bhp', value(bhp_ad), 'qS', value(qS_ad));

%% Main loop
t = 0; step = 0;
hwb = waitbar(t,'Simulation ..');
while t < totTime
   t = t + dt;
   step = step + 1;
   fprintf('\nTime step %d: Time %.2f -> %.2f days\n', ...
      step, convertTo(t - dt, day), convertTo(t, day));
   % Newton loop
   resNorm = 1e99;
   p0  = value(p_ad); % Previous step pressure
   nit = 0;
   while (resNorm > tol) && (nit <= maxits)
      % Add source terms to homogeneous pressure equation:
      eq1     = presEq(p_ad, p0, dt);
      eq1(wc) = eq1(wc) - q_conn(p_ad, bhp_ad);
      % Collect all equations
      eqs = {eq1, rateEq(p_ad, bhp_ad, qS_ad), ctrlEq(bhp_ad)};
      % Concatenate equations and solve for update:
      eq  = cat(eqs{:});
      J   = eq.jac{1};  % Jacobian
      res = eq.val;     % residual
      upd = -(J \ res); % Newton update
      % Update variables
      p_ad.val   = p_ad.val   + upd(pIx);
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
                            'bhp', value(bhp_ad), 'qS', value(qS_ad));
      waitbar(t/totTime,hwb);
   end
end
close(hwb);
%% Plot production rate and pressure decay
clf,
[ha,hr,hp] = plotyy(...
   [sol(2:end).time]/day, -[sol(2:end).qS]*day, ...
   [sol(2:end).time]/day, mean([sol(2:end).pressure]/barsa), 'stairs', 'plot');
set(ha,'FontSize',16);
set(hr,'LineWidth', 2);
set(hp,'LineStyle','none','Marker','o','LineWidth', 1);
set(ha(2),'YLim',[100 210],'YTick',100:50:200);
xlabel('time [days]');
ylabel(ha(1), 'rate [m^3/day]');
ylabel(ha(2), 'avg pressure [bar]');

%% Plot pressure evolution
clf;
steps = [2 5 10 20];
for i=1:4
   subplot(2,2,i);
   set(gca,'Clipping','off');
   plotCellData(G, sol(steps(i)).pressure/barsa, show,'EdgeColor',.5*[1 1 1]);
   plotWell(G,W);
   view(-125,20), camproj perspective
   caxis([115 205]);
   axis tight off; zoom(1.1)
   text(200,170,-8,[num2str(round(steps(i)*dt/day)) ' days'],'FontSize',14);
end
h=colorbar('horiz','Position',[.1 .05 .8 .025]);
colormap(jet(55));

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
