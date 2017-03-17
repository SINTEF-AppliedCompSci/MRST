
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

p = linspace(100*barsa,220*barsa,50)';
s = linspace(0,1,50)';
plot(p/barsa, pv_r(1).*exp(cr*(p-p_r)),'LineWidth',2);


%% Define model for two-phase compressible fluid
% Define a water phase
muW    = 1*centi*poise;
cw     = 1e-5/barsa;
rho_rw = 960*kilogram/meter^3;
rhoWS  = 1000*kilogram/meter^3;
rhoW   = @(p) rho_rw .* exp( cw * (p - p_r) );
krW = @(S) S.^2;

% Define a lighter, more viscous oil phase with different relative
% permeability function
muO   = 5*centi*poise;
co      = 1e-4/barsa;
rho_ro = 850*kilogram/meter^3;
rhoOS  = 750*kilogram/meter^3;
krO = @(S) S.^3;

rhoO   = @(p) rho_ro .* exp( co * (p - p_r) );
figure;
plot(p/barsa, [rhoW(p), rhoO(p)],'LineWidth',2);
legend('Water density', 'Oil density')

figure;
plot(p/barsa, [krW(s), krO(s)],'LineWidth',2);
legend('krW', 'krO')
%% Impose vertical equilibrium
gravity reset on, g = norm(gravity);
[z_0, z_max] = deal(0, max(G.cells.centroids(:,3)));
equil  = ode23(@(z,p) g .* rhoO(p), [z_0, z_max], p_r);
p_init = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil
sW_init = zeros(G.cells.num, 1);
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
grad = @(x)C*x; % Discrete gradient
div  = @(x)-C'*x; % Discrete divergence
avg  = @(x) 0.5 * (x(N(:,1)) + x(N(:,2))); % Averaging
upw = @(flag, x) flag.*x(N(:, 1)) + ~flag.*x(N(:, 2)); % Upwinding
spy(C)

gradz  = grad(G.cells.centroids(:,3));
%% Initialize for solution loop
[p_ad, sW_ad] = initVariablesADI(p_init, sW_init);
nc = G.cells.num;
pIx = 1:nc;
sIx = (nc+1):(2*nc);

numSteps = 100;                  % number of time-steps
totTime  = 365*day;             % total simulation time
dt       = totTime / numSteps;  % constant time step
tol      = 1e-5;                % Newton tolerance
maxits   = 15;                  % max number of Newton its

injIndex = 1;
prodIndex = G.cells.num;

inRate = 1*sum(pv(p_init))/totTime;
outRate = 0.5*inRate;

sol = repmat(struct('time',[],'pressure',[], 's', []),[numSteps+1,1]);
sol(1)  = struct('time', 0, 'pressure', double(p_ad), ...
    's', double(sW_ad));

%% Main loop
t = 0; step = 0;
hwb = waitbar(t,'Simulation ..');
while t < totTime,
   t = t + dt;
   step = step + 1;
   fprintf('\nTime step %d: Time %.2f -> %.2f days\n', ...
      step, convertTo(t - dt, day), convertTo(t, day));
   % Newton loop
   resNorm = 1e99;
   p0  = double(p_ad); % Previous step pressure
   sW0 = double(sW_ad);
   nit = 0;
   while (resNorm > tol) && (nit <= maxits)
      % Evaluate properties
      rW = rhoW(p_ad);
      rW0 = rhoW(p0);
      rO = rhoO(p_ad);
      rO0 = rhoO(p0);
      
      % Define pressure drop over interface for both phases
      dp = grad(p_ad);
      dpW = dp - g*avg(rW).*gradz; % Water
      dpO = dp - g*avg(rO).*gradz; % Oil
      % Pore volume of cells at current pressure and previous pressure 
      vol0 = pv(p0);
      vol = pv(p_ad);
      % Mobility: Relative permeability over constant viscosity
      mobW = krW(sW_ad)./muW;
      mobO = krO(1-sW_ad)./muO;
      % Define phases fluxes. Density and mobility is taken upwinded (value
      % on interface is defined as taken from the cell where the phase flow
      % is currently coming from). This gives more physical solutions than
      % averaging or downwinding.
      vW = -upw(double(dpW) <= 0, rW.*mobW).*T.*dpW;
      vO = -upw(double(dpO) <= 0, rO.*mobO).*T.*dpO;
      % Conservation of water and oil
      water = (1/dt).*(vol.*rW.*sW_ad - vol0.*rW0.*sW0) + div(vW);
      oil   = (1/dt).*(vol.*rO.*(1-sW_ad) - vol0.*rO0.*(1-sW0)) + div(vO);
      % Insert volumetric source term multiplied by density
      water(injIndex) = water(injIndex) - inRate.*rhoWS;
      % Set production cells to fixed pressure of 200 bar and zero water
      water(prodIndex) = p_ad(prodIndex) - 200*barsa;
      oil(prodIndex) = sW_ad(prodIndex);
      % Collect all equations
      eqs = {oil, water};
      % Concatenate equations and solve for update:
      eq  = cat(eqs{:});
      J   = eq.jac{1};  % Jacobian
      res = eq.val;     % residual
      upd = -(J \ res); % Newton update
      % Update variables
      p_ad.val   = p_ad.val   + upd(pIx);
      sW_ad.val = sW_ad.val + upd(sIx);
      sW_ad.val = min(sW_ad.val, 1);
      sW_ad.val = max(sW_ad.val, 0);
      
      resNorm = norm(res);
      nit     = nit + 1;
      fprintf('  Iteration %3d:  Res = %.4e\n', nit, resNorm);
   end

   if nit > maxits,
      error('Newton solves did not converge')
   else % store solution
      sol(step+1)  = struct('time', t, ...
                            'pressure', double(p_ad), ...
                            's', double(sW_ad));
      waitbar(t/totTime,hwb);
   end
end
close(hwb);

%% Plot pressure evolution

for i = 1:numSteps
    figure(1); clf
    subplot(2, 1, 1)
    plotCellData(G, sol(i).pressure);
    title('Pressure')
    view(30, 40);
    subplot(2, 1, 2)
    plotCellData(G, sol(i).s);
    caxis([0, 1])
    view(30, 40);
    title('sW')
    drawnow
end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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