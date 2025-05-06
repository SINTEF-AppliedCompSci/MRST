%% Conceptual three-phase compressible model
% The purpose of the example is to show how one can implement compressible
% three-phase flow models using automatic differentiation and the abstract
% differentiation operators in MRST. To this end, we set up a relatively
% simple model in which the upper 1/3 is initially filled with gas, the
% middle 1/3 is filled with oil, and the bottom 1/3 is filled with water.
% Water is injected in the lower southwest corner and fluids are produced
% from the northeast corner of the middle layer. To mimic the effect of
% injector and producer wells, we manipulate the flow equations to ensure
% reasonable inflow/outflow for the perforated cells.
%
% The code can be run with or without water injection. Relative
% permeabilities are sampled from the SPE3 benchmark. Physical properties
% are in generally chosen quite haphazardly for illustration
% purposes. The last plots use 'hold all' so that if you run the function
% twice with or without water injection, the plots of pressures and oil and
% gas production will be added to the same plots.
mrstModule add ad-core ad-props deckformat

if exist('pvi','var')~=1
    pvi = .5;
end
useInj = (pvi>0);

%% Set up model geometry
[nx,ny,nz] = deal(  11,  11,  3);
[Dx,Dy,Dz] = deal(200, 200, 50);
G = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
G = computeGeometry(G);

%% Define rock model
rock = makeRock(G, 30*milli*darcy, 0.2);
cr   = 1e-6/psia;
pR   = 200*barsa;
pvR  = poreVolume(G, rock);
pv   = @(p) pvR .* exp(cr * (p - pR) );

%% Fluid model
% Water is slightly compressible with quadratic relative permeability
muW    = 1*centi*poise;
cw     = 2e-6/psia;
rhoWR  = 1014*kilogram/meter^3;
rhoWS  = 1000*kilogram/meter^3;
rhoW   = @(p) rhoWR .* exp( cw * (p - pR) );

% Oil phase is lighter and has a cubic relative permeability
muO    = 0.5*centi*poise;
co     = 1e-5/psia;
rhoOR  = 850*kilogram/meter^3;
rhoOS  = 750*kilogram/meter^3;
rhoO   = @(p) rhoOR .* exp( co * (p - pR) );

% Gas
muG    = 0.015*centi*poise;
cg     = 1e-3/psia;
rhoGR  = 1.2;
rhoGS  = 1;
rhoG   = @(p) rhoGR .* exp( cg * (p - pR) );

% Get relative permeabilities from the SPE3 data set
pth  = getDatasetPath('spe3');
fn   = fullfile(pth, 'BENCH_SPE3.DATA');
deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);
f    = initDeckADIFluid(deck);
f    = assignRelPerm(f);

%% Extract grid information
nf = G.faces.num;                                 % Number of faces
nc = G.cells.num;                                 % Number of cells
N  = double(G.faces.neighbors);                   % Map: face -> cell
F  = G.cells.faces(:,1);                          % Map: cell -> face
intInx = all(N ~= 0, 2);                          % Interior faces
N  = N(intInx, :);                                % Interior neighbors

%% Compute transmissibilities
hT = computeTrans(G, rock);                       % Half-transmissibilities
T  = 1 ./ accumarray(F, 1 ./ hT, [nf, 1]);        % Harmonic average
T  = T(intInx);                                   % Restricted to interior
nf = size(N,1);

%% Define discrete operators
D     = sparse([(1:nf)'; (1:nf)'], N, ones(nf,1)*[-1 1], nf, nc);
grad  = @(x)  D*x;
div   = @(x)-D'*x;
avg   = @(x) 0.5 * (x(N(:,1)) + x(N(:,2)));
upw   = @(flag, x) flag.*x(N(:, 1)) + ~flag.*x(N(:, 2));
gradz = grad(G.cells.centroids(:,3));

%% Impose vertical equilibrium
gravity reset on,
g     = norm(gravity);
equil = ode23(@(z,p) g.* rhoO(p), [0, max(G.cells.centroids(:,3))], pR);
p0    = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil
sW0   = f.krPts.w(1)*ones(G.cells.num, 1);
sW0   = reshape(sW0,G.cartDims); sW0(:,:,nz)=1; sW0 = sW0(:);
sG0   = zeros(G.cells.num, 1);
sG0   = reshape(sG0,G.cartDims); sG0(:,:, 1)=1-f.krPts.w(1); sG0 = sG0(:);

ioip  = sum(pv(p0).*(1-sW0 - sG0).*rhoO(p0));
igip  = sum(pv(p0).*sG0.*rhoG(p0));

%% Schedule and injection/production
nstep = 108;
Tf = 1080*day;
dt = Tf/nstep*ones(1,nstep);
dt = [dt(1).*sort(repmat(2.^-[1:5 5],1,1)) dt(2:end)];
nstep = numel(dt);

[inRate,  inIx ] = deal(pvi*sum(pv(p0))/Tf, nx*ny*nz-nx*ny+1);
[outPres, outIx] = deal(100*barsa, 2*nx*ny);

%% Initialize for solution loop
[p, sW, sG]     = initVariablesADI(p0, sW0, sG0);
[pIx, wIx, gIx] = deal(1:nc, (nc+1):(2*nc), (2*nc+1):(3*nc));
[tol, maxits]   = deal(1e-5,15);
pargs = {};% {'EdgeColor','k','EdgeAlpha',.1};

sol = repmat(struct('time',[],'pressure',[], 'sW', [], 'sG', []),[nstep+1,1]);
sol(1) = struct('time', 0, 'pressure', value(p),'sW', value(sW), 'sG', value(sG));

%% Main loop
t = 0;
hwb = waitbar(t,'Simulation ..');
nits = nan(nstep,1);
i = 1;
for n=1:nstep

    t = t + dt(n);
    fprintf('\nTime step %d: Time %.2f -> %.2f days\n', ...
        n, convertTo(t - dt(n), day), convertTo(t, day));

    % Newton loop
    resNorm   = 1e99;
    [p0, sW0, sG0] = deal(value(p), value(sW),value(sG));
    nit = 0;
    while (resNorm > tol) && (nit <= maxits)

        % Densities and pore volumes
        [rW,rW0,rO,rO0,rG,rG0] = ....
            deal(rhoW(p), rhoW(p0), rhoO(p), rhoO(p0), rhoG(p), rhoG(p0));
        [vol, vol0]     = deal(pv(p), pv(p0));

        % Mobility: Relative permeability over constant viscosity
        [krW, krO, krG] = f.relPerm(sW,sG);
        mobW = krW./muW;
        mobO = krO./muO;
        mobG = krG./muG;

        % Pressure differences across each interface
        dp  = grad(p);
        dpW = dp-g*avg(rW).*gradz;
        dpO = dp-g*avg(rO).*gradz;
        dpG = dp-g*avg(rG).*gradz;

        % Phase fluxes:
        % Density and mobility is taken upwinded (value on interface is
        % defined as taken from the cell from which the phase flow is
        % currently coming from). This gives more physical solutions than
        % averaging or downwinding.
        vW = -upw(value(dpW) <= 0, rW.*mobW).*T.*dpW;
        vO = -upw(value(dpO) <= 0, rO.*mobO).*T.*dpO;
        vG = -upw(value(dpG) <= 0, rG.*mobG).*T.*dpG;

        % Conservation of water and oil
        water = (1/dt(n)).*(vol.*rW.*sW - vol0.*rW0.*sW0) + div(vW);
        oil   = (1/dt(n)).*(vol.*rO.*(1-sW-sG) - vol0.*rO0.*(1-sW0-sG0)) + div(vO);
        gas   = (1/dt(n)).*(vol.*rG.*sG - vol0.*rG0.*sG0) + div(vG);

        % Injector: volumetric source term multiplied by surface density
        if useInj
            water(inIx) = water(inIx) - inRate.*rhoWS;
        end

        % Producer: replace equations by new ones specifying fixed pressure
        % and zero water saturation
        oil(outIx) = p(outIx) - outPres;
        water(outIx) = sW(outIx);
        gas(outIx)   = sG(outIx);

        % Collect and concatenate all equations (i.e., assemble and
        % linearize system)
        eqs = {oil, water, gas};
        eq  = cat(eqs{:});

        % Compute Newton update and update variable
        res = eq.val;
        upd = -(eq.jac{1} \ res);
        p.val  = p.val   + upd(pIx);
        sW.val = sW.val + upd(wIx);
        %sW.val = max( min(sW.val, 1), 0);
        sG.val = sG.val + upd(gIx);
        %sG.val = max( min(sG.val, 1), 0);

        resNorm = norm(res);
        nit     = nit + 1;
        fprintf('  Iteration %3d:  Res = %.4e\n', nit, resNorm);
    end
    if nit > maxits
        error('Newton solves did not converge')
    else % plot
        nits(n) = nit;

        figure(1); clf
        subplot(2,1,1)
        plotCellData(G, value(p)/barsa, pargs{:});
        title('Pressure'), view(30, 40);

        subplot(2,1,2)
        sg = value(sG); sw = value(sW);
        plotCellData(G,[sg, 1-sg-sw, sw]/1.001+1e-4,pargs{:});             % Avoid color artifacts
        caxis([0, 1]), view(30, 40); title('Saturation')
        drawnow

        sol(n+1) = struct('time', t, ...
                          'pressure', value(p), ...
                          'sW', value(sW), ...
                          'sG', value(sG));
        waitbar(t/Tf,hwb);
        pause(.1);
    end
end
close(hwb);

%%
figure(2),
subplot(2,1,2-value(useInj));
bar(nits,1,'EdgeColor','r','FaceColor',[.6 .6 1]);
axis tight, set(gca,'YLim',[0 10]);

%%
figure(3); hold all
ipr = arrayfun(@(x) x.pressure(inIx), sol)/barsa;
t   = arrayfun(@(x) x.time, sol)/day;
fpr = arrayfun(@(x) sum(x.pressure)/G.cells.num, sol)/barsa;
plot(t,ipr,'-', t,fpr, '--');

%%
figure(4);
oip = arrayfun(@(x) sum(pv(x.pressure).*(1-x.sW-x.sG).*rhoO(x.pressure)), sol);
gip = arrayfun(@(x) sum(pv(x.pressure).*x.sG.*rhoG(x.pressure)), sol);
subplot(1,2,1), hold all
plot(t,(ioip-oip)./rhoOS)
subplot(1,2,2), hold all
plot(t,(igip-gip)./rhoGS)


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
