%% Single phase flow simulation using AD
% This example goes through the steps of setting up a single-phase
% simulation with a single horizontal well using the automatic
% differentiation framework.


% Required modules
mrstModule add ad-fi ad-props

% Setup 10x10x10 grid of 200x200x50 m model.
nx = 10;    ny = 10;    nz = 10;
Dx = 200;   Dy = 200;   Dz = 50;
G = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
G = computeGeometry(G);

%% Setup rock properties
%

% Assume homogeneous/isotropic rock.
permX = 30*milli*darcy;
poro  = 0.3;
rock.perm = repmat(permX, [G.cells.num, 1]);
rock.poro = repmat(poro , [G.cells.num, 1]);

% Set rock compressibility:
cr = 1e-6/barsa;

%%
% In the case of non-zero rock-compressiblity cr, the input rock porosity
% is taken as reference at a given reference pressure p_r. The grid pore
% volume (pv) becomes a function of pressure given by the differential
% equation
%                 cr = (d pv/d p)/pv
% which results in
%                 pv(p) = pv_r e^( cr(p-p_r) )
% in which pv_r is the reference pore volume (rock.poro x volume) and p_r
% is the reference pressure. We assume the reference pressure is 200 Bar:

pv_r = poreVolume(G, rock);
p_r  = 200*barsa;

% Finally, the pressure-dependent function for pore-volumes becomes:

pv   = @(p) pv_r .* exp( cr .* (p - p_r) );

%% Fluid (oil) properties
% We Assume constant viscosity:
mu   = 5*centi*poise;
% Assume that the oil compressibility can be approximated as constant in
% the reservoir:
c    = 1e-3/barsa;
% When the compressibility is constant, the fluid density becomes a
% function of pressure given by the differential equation
%                 c = (d rho/d p)/rho
% which results in
%                 rho(p) = rho_r e^( c(p-p_r) )
% in which rho_r is the reference fluid density at reference pressure p_r.
% We assume that rho_r = 850 kg/m^3 at p_r = 200 Bar:
p_r   = 200*barsa;
rho_r = 850*kilogram/meter^3;

% Finally define the pressure dependent function for rho:
rho   = @(p) rho_r .* exp( c .* (p - p_r) );

% Computing surface volume rates requires fluid density at surface
% conditions. We assume it is 750 kg/m^3:
rhoS = 750*kilogram/meter^3;

%% Single horizontal well in J direction
% We consider a well with 8 connections in the J direction. 

W = [];
nperf = 8;
I = repmat(2, [nperf, 1]);
J = (1 : nperf).' + 1;
K = repmat(5, [nperf, 1]);

% Convert IJK-indices to linear index (as used in G)
cellInx = sub2ind(G.cartDims, I, J, K);

W = addWell(W, G, rock, cellInx, 'Name', 'producer', 'InnerProduct', 'ip_tpf');

% Plotting
f = [figure(1), figure(2)];
set(0, 'CurrentFigure', f(1));
clf
plotGrid(G, 'FaceColor', 'g', 'FaceAlpha', .3, 'EdgeColor', 'w');
plotWell(G, W);
axis off;
set(f(1), 'Color', 'w'); 
camproj perspective;
view(3);

%% Initial conditions
% We assume that the reservoir is initially at equilibrium. This means that
% the following condition must be satisfied:
%          dp/dz = g rho,
% where g is the gavitational accelleration. This relation can be solved
% analytically for p, but alternatively one can solve the above ODE with
% 'initial condtition' p(z_0) = p_r:

gravity reset on
g = norm(gravity);
[z_0, z_max] = deal(0, max(G.cells.centroids(:,3)) + 5);
press  = ode23(@(z,p) g .* rho(p), [z_0, z_max], p_r);
% We then interpolate onto the grid using cell centers:
p_init = reshape(deval(press, G.cells.centroids(:,3)), [], 1);  clear press

%% Setting up components needed for the simulation
% Since we impose no-flow boundary conditions in this example, we restrict
% connections to interior faces only.
N  = double(G.faces.neighbors);
intInx = all(N ~= 0, 2);
N  = N(intInx, :);
% We will be using the two-point flux approximation. First the one-sided
% transmissibilities are computed, then the harmonic average is taken to
% obtain the two-sided transmissibilites.
hT = computeTrans(G, rock);
cf = G.cells.faces(:,1);
nf = G.faces.num;
T  = 1 ./ accumarray(cf, 1 ./ hT, [nf, 1]);
T  = T(intInx);
% In setting up the equations, we need discrete forms of the divergence and
% gradient operators, and we represent these as multiplication by sparse
% matrices. In particular, we construct the 'gradient matrix' C as folows:
n = size(N,1);
C = sparse( [(1:n)'; (1:n)'], N, ones(n,1)*[1 -1], n, G.cells.num);
% The discrete gradient and divergence operators are now given by
grad = @(x)-C*x;
div  = @(x)C'*x;

% Additionally, we will need to take the average of cell-based quantities
% in neighboring cells and define the following function
avg = @(x) 0.5 * (x(N(:,1)) + x(N(:,2)));

%% Pressure and well equations:
% The pressure equation (without well contributions) is given by:
%
% $\frac{d}{dt}(\phi\rho)+\nabla\cdot(\rho v)=0, \quad v = -\frac{K}{\mu}\nabla(p-g\rho z)$
%
% In discretized form, this leads to
z = G.cells.centroids(:,3); % z-coordinate of grid cells
pressureEq = @(p, p0, dt) (1/dt) .* (pv(p).*rho(p) - pv(p0).*rho(p0)) ...
    - div( avg(rho(p) ./ mu) .* T .* grad(p - g*rho(p).*z) );

% Wellrates are given as Peaceman well-index times pressure drop
wc = W(1).cells; % perforation grid cells
WI = W(1).WI;    % well indices
dz = W(1).dZ;    % perforation depth relative to well reference depth
wellRates = ...
   @(p, bhp) WI .* (rho(p(wc)) ./ mu) .* (bhp - p(wc) + g*dz.*rho(p(wc)));

%% Define ADI variables
% The primary variables are grid-cell pressures, well bhp and surface
% rate.
[p_ad, bhp_ad, qS_ad] = initVariablesADI(p_init, p_init(wc(1)), 0);

% For convenience, create indices to variables when stacked:
pIx = 1:G.cells.num; bhpIx = G.cells.num + 1; qSIx = G.cells.num + 2;

%% Set up simulation parameters
%

numSteps = 52;
totTime  = 365*day;
dt       = totTime / numSteps;
% Tolerance and maximum number of iterations for the Newton solver.
tol      = 1e-5; 
maxits   = 10;

% Save output in array 'sol'
sol = repmat(struct('time', [], 'pressure', [], 'bhp', [], 'qS', []), ...
             [numSteps + 1, 1]);

sol(1).time     = 0;
sol(1).pressure = value(p_ad);
sol(1).bhp      = value(bhp_ad);
sol(1).qS       = value(qS_ad);

% Set up plot
set(0, 'CurrentFigure', f(2));
clf
set(f(2), 'Color', 'w')
subplot(2,1,1),
plotCellData(G, convertTo(p_init, barsa));
title('pressure [bar]', 'EdgeColor', 'w');
colorbar; view(3);
camproj perspective

subplot(2,1,2);
axis([0, convertTo(totTime,day), 0, 300]);
title('Surface volume rate [m^3/day]'); 
hold on


%% Main simulation
% We solve the equations implicitely. At each time step, the equations are assembled and
% the automatic differentiation framework takes care automatically of the computation of
% the Jacobian.

t = 0; 
step = 0;
nDigits = floor(log10(maxits)) + 1;
while t < totTime
    t = t + dt;
    step = step + 1;
    fprintf('\nTime step %d: Time %.2f -> %.2f days\n', ...
            step, convertTo(t - dt, day), convertTo(t, day));

    % Newton loop
    resNorm = 1e99;
    p0  = value(p_ad); % Previous step pressure
    nit = 0;
    while (resNorm > tol) && (nit < maxits)
        % Create equations:
        eqs = cell([3, 1]);

        eqs{1} = pressureEq(p_ad, p0, dt);
        % Add well contributions in perforated cells:
        eqs{1}(wc) = eqs{1}(wc) - wellRates(p_ad, bhp_ad);

        % Sum of wellrates should equal total rate:
        eqs{2} = qS_ad - sum(wellRates(p_ad, bhp_ad))/rhoS;

        % Final equation is prescribed bhp
        eqs{3} = bhp_ad - 100*barsa;

        % Concatenate equations and solve:
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
        fprintf('  Iteration %*d:  Res = %.4e\n', nDigits, nit, resNorm);
    end

    if nit > maxits
        error('Newton solves did not converge')
    else
        sol(step+1).time     = t;
        sol(step+1).pressure = value(p_ad);
        sol(step+1).bhp      = value(bhp_ad);
        sol(step+1).qS       = value(qS_ad);

        % Plot evolution
        set(0, 'CurrentFigure', f(2));
        subplot(2,1,1), cla, caxis([120, 205])
        plotCellData(G, convertTo(sol(step+1).pressure, barsa), 'EdgeColor', 'w');

        subplot(2,1,2)
        plot(convertTo(sol(step+1).time, day), ...
             convertTo(-sol(step+1).qS , meter^3/day), '*');

        drawnow
    end
end

% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
