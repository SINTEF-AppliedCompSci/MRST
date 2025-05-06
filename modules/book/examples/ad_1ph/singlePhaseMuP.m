%% Pressure-dependent viscosity
% In this example, we will demonstrate how one can easily extend the
% compressible single-phase pressure solver to include the effect of
% pressure-dependent viscosity using either arithmetic averaging of the
% viscosity or harmonic averaging of the fluid mobility.

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
grad = @(x)C*x;
div  = @(x)-C'*x;
avg  = @(x) 0.5 * (x(N(:,1)) + x(N(:,2)));

%% Rock model and transmissibilities
% rock = makeRock(G, 30*milli*darcy, 0.3);
mrstModule add spe10
rock = getSPE10rock(41:50,101:110,1:10);

cr   = 1e-6/barsa;
p_r  = 200*barsa;
pv_r = poreVolume(G, rock);
pv   = @(p) pv_r .* exp( cr * (p - p_r) );

hT = computeTrans(G, rock);
cf = G.cells.faces(:,1);
nf = G.faces.num;
T  = 1 ./ accumarray(cf, 1 ./ hT, [nf, 1]);
T  = T(intInx);

%% Fluid model
mu0   = 5*centi*poise;
c     = 1e-3/barsa;
rho_r = 850*kilogram/meter^3;
rhoS  = 750*kilogram/meter^3;
rho   = @(p) rho_r .* exp( c * (p - p_r) );

gravity reset on, g = norm(gravity);
[z_0, z_max] = deal(0, max(G.cells.centroids(:,3)));
equil  = ode23(@(z,p) g .* rho(p), [z_0, z_max], p_r);
p_init = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil

%% Assume a single horizontal well
nperf = 8;
I = repmat(2, [nperf, 1]);
J = (1 : nperf).' + 1;
K = repmat(5, [nperf, 1]);
cellInx = sub2ind(G.cartDims, I, J, K);
W = addWell([ ], G, rock, cellInx, 'Name', 'P1', 'Dir', 'x' );

%% Main loop
numSteps = 52;
totTime  = 365*day;
dt       = totTime / numSteps;
tol      = 1e-5;
maxits   = 10;
mu_const = [0 2e-3 5e-3]/barsa;
[mypres,myrate] = deal(nan(numSteps+1,2*numel(mu_const)));

n = 1;
for method=1:2

   for i=1:numel(mu_const)

      mu = @(p) mu0*(1+mu_const(i)*(p-p_r));

      gradz = grad(G.cells.centroids(:,3));
      switch method
         case 1
            v  = @(p)  -(T./mu(avg(p))).*( grad(p) - g*avg(rho(p)).*gradz );
         case 2
            hf2cn = getCellNoFaces(G);
            nhf = numel(hf2cn);
            hf2f  = sparse(double(G.cells.faces(:,1)),(1:nhf)',1);
            hf2if = hf2f(intInx,:);
            hlam = @(mu,p) 1./(hf2if*(mu(p(hf2cn))./hT));
            %
            v  = @(p) -hlam(mu,p).*( grad(p) - g*avg(rho(p)).*gradz );
      end

      presEq = @(p, p0, dt) (1/dt)*(pv(p).*rho(p) - pv(p0).*rho(p0)) ...
         + div( avg(rho(p)).*v(p) );

      % Define well equations
      wc = W(1).cells; % connection grid cells
      WI = W(1).WI;    % well-indices
      dz = W(1).dZ;    % connection depth relative to bottom-hole

      p_conn  = @(bhp)  bhp + g*dz.*rho(bhp); %connection pressures
      q_conn  = @(p,bhp) WI .* (rho(p(wc)) ./ mu(p(wc))) .* (p_conn(bhp) - p(wc));

      rateEq = @(p,bhp,qS)  qS-sum(q_conn(p, bhp))/rhoS;
      ctrlEq = @(bhp)       bhp-100*barsa;

      % Initialize for solution loop
      [p_ad, bhp_ad, qS_ad] = initVariablesADI(p_init, p_init(wc(1)), 0);
      nc = G.cells.num;
      [pIx, bhpIx, qSIx] = deal(1:nc, nc+1, nc+2);
      sol = repmat(struct('time',[],'pressure',[],'bhp',[],'qS',[]),[numSteps+1,1]);
      sol(1)  = struct('time', 0, 'pressure', value(p_ad), ...
         'bhp', value(bhp_ad), 'qS', value(qS_ad));

      % Time loop
      t = 0; step = 0;
      hwb = waitbar(0,'Simulation ..');
      while t < totTime
         t = t + dt;
         step = step + 1;
         fprintf('Time step %d: Time %.2f -> %.2f days\n', ...
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
          end

         if nit > maxits
            error('Newton solves did not converge')
         else % store solution
            sol(step+1)  = struct('time', t, 'pressure', value(p_ad), ...
               'bhp', value(bhp_ad), 'qS', value(qS_ad));
           waitbar(t/totTime,hwb)
         end
      end
      close(hwb)
      myrate(2:end,n) = -[sol(2:end).qS]*day;
      mypres(:,n) = mean([sol(:).pressure]/barsa);
      mytime = [sol(1:end).time]/day;
      n = n+1;
   end
end

%%
figure(1);
myrate(1,:)=NaN;
stairs(mytime, myrate(:,1:3),'LineWidth',2);
set(gca,'FontSize',16);
xlabel('time [days]'); ylabel('rate [m^3/day]');
legend('c_\mu=0', 'c_\mu=0.002', 'c_\mu=0.005');

figure(2);
plot(mytime, mypres(:,1:3),'o');
set(gca,'FontSize',16);
xlabel('time [days]'); ylabel('avg pressure [bar]');
legend('c_\mu=0', 'c_\mu=0.002', 'c_\mu=0.005');

figure(3);
myrate(1,:)=0;
plot(mytime, cumsum(myrate(:,1:3)),'LineWidth',2);
set(gca,'FontSize',16);
xlabel('time [days]'); ylabel('cummulative product [m^3]');
legend('c_\mu=0', 'c_\mu=0.002', 'c_\mu=0.005', 'Location', 'SouthEast');

figure(4);
p = [100*barsa,205*barsa];
mup = zeros(numel(p), numel(mu_const));
for i=1:numel(mu_const)
   mu = @(p) mu0*(1+mu_const(i)*(p-p_r));
   mup(:,i) = mu(p);
end
plot(convertTo(p,barsa),convertTo(mup,centi*poise),'LineWidth',2);
set(gca,'FontSize',16);
xlabel('pressure [bar]'); ylabel('viscosity [cP]');
legend('c_\mu=0', 'c_\mu=0.002', 'c_\mu=0.005', 'Location', 'SouthEast');
axis tight

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
