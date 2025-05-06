%% Conceptual two-phase compressible model
% The purpose of the example is to show how one can implement compressible
% multiphase flow models using automatic differentiation and the abstract
% differentiation operators in MRST. To this end, we set up a relatively
% simple model in which the upper 2/3 is initially filled with oil and the
% bottom 1/3 is filled with water. Water is injected in the lower southwest
% corner and fluids are produced from the upper northeast corner. To mimic
% the effect of injector and producer wells, we manipulate the flow
% equations to ensure reasonable inflow/outflow for the perforated cells.
%
% The code can be run with compressible or incompressible fluids and rock.
% (The compressibility of oil is somewhat exaggerated for illustration
% purposes.) The last plots use 'hold all' so that if you run the function
% twice with compressible/incompressible model, the plots of pressures,
% condition numbers, and oil production will be added to the same plots.

if exist('isCompr','var')~=1
    isCompr = true;
end

%% Set up model geometry
[nx,ny,nz] = deal(  11,  11,  3);
[Dx,Dy,Dz] = deal(200, 200, 50);
G = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
G = computeGeometry(G);

figure(1); clf
plotGrid(G); view(3); axis tight

%% Define rock model
rock = makeRock(G, 30*milli*darcy, 0.2);
cr   = 1e-6/psia*double(isCompr);
pR   = 200*barsa;
pvR  = poreVolume(G, rock);
pv   = @(p) pvR .* exp(cr * (p - pR) );

p = linspace(100*barsa,250*barsa,30)';
s = linspace(0,1,50)';
figure(1), clf
plot(p/barsa, pvR(1).*exp(cr*(p-pR)),'LineWidth',2);


%% Fluid model
% Water is slightly compressible with quadratic relative permeability
muW    = 1*centi*poise;
cw     = 2e-6/psia*double(isCompr);
rhoWR  = 1014*kilogram/meter^3;
rhoWS  = 1000*kilogram/meter^3;
rhoW   = @(p) rhoWR .* exp( cw * (p - pR) );
krW    = @(S) S.^2;

% Oil phase is lighter and has a cubic relative permeability
muO    = 5*centi*poise;
co     = 1e-5/psia*double(isCompr);
rhoOR  = 850*kilogram/meter^3;
rhoOS  = 750*kilogram/meter^3;
rhoO   = @(p) rhoOR .* exp( co * (p - pR) );
krO    = @(S) S.^3;

% Plot compressibility factors
figure(1), subplot(1,2,1)
plot(p/barsa, [rhoW(p), rhoO(p)],'LineWidth',2);
legend('Water density', 'Oil density','Location','best')

% Plot relative permeability
subplot(1,2,2);
plot(s, [krW(s), krO(s)],'LineWidth',2);
legend('krW', 'krO','Location','NorthWest');

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
spy(D)

%% Impose vertical equilibrium
gravity reset on,
g     = norm(gravity);
equil = ode23(@(z,p) g.* rhoO(p), [0, max(G.cells.centroids(:,3))], pR);
p0    = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil
sW0   = zeros(G.cells.num, 1);
sW0   = reshape(sW0,G.cartDims); sW0(:,:,nz)=1; sW0 = sW0(:);
ioip  = sum(pv(p0).*(1-sW0).*rhoO(p0));

%% Schedule and injection/production
nstep = 36;
Tf = 360*day;
dt = Tf/nstep*ones(1,nstep);
dt = [dt(1).*sort(repmat(2.^-[1:5 5],1,1)) dt(2:end)];
nstep = numel(dt);

[inRate,  inIx ] = deal(.1*sum(pv(p0))/Tf, nx*ny*nz-nx*ny+1);
[outPres, outIx] = deal(100*barsa, nx*ny);

%% Initialize for solution loop
[p, sW]    = initVariablesADI(p0, sW0);
[pIx, sIx] = deal(1:nc, (nc+1):(2*nc));
[tol, maxits] = deal(1e-5,15);
pargs = {};% {'EdgeColor','k','EdgeAlpha',.1};

sol = repmat(struct('time',[],'pressure',[], 's', []),[nstep+1,1]);
sol(1) = struct('time', 0, 'pressure', value(p),'s', value(sW));

%% Main loop
t = 0;
hwb = waitbar(t,'Simulation ..');
nits = nan(nstep,1);
i = 1; clear cnd;
for n=1:nstep

    t = t + dt(n);
    fprintf('\nTime step %d: Time %.2f -> %.2f days\n', ...
        n, convertTo(t - dt(n), day), convertTo(t, day));

    % Newton loop
    resNorm   = 1e99;
    [p0, sW0] = deal(value(p), value(sW));
    nit = 0;
    while (resNorm > tol) && (nit <= maxits)

        % Densities and pore volumes
        [rW,rW0,rO,rO0] = deal(rhoW(p), rhoW(p0), rhoO(p), rhoO(p0));
        [vol, vol0]     = deal(pv(p), pv(p0));

        % Mobility: Relative permeability over constant viscosity
        mobW = krW(sW)./muW;
        mobO = krO(1-sW)./muO;

        % Pressure differences across each interface
        dp  = grad(p);
        dpW = dp-g*avg(rW).*gradz;
        dpO = dp-g*avg(rO).*gradz;

        % Phase fluxes:
        % Density and mobility is taken upwinded (value on interface is
        % defined as taken from the cell from which the phase flow is
        % currently coming from). This gives more physical solutions than
        % averaging or downwinding.
        vW  = -upw(value(dpW) <= 0, rW.*mobW).*T.*dpW;
        vO  = -upw(value(dpO) <= 0, rO.*mobO).*T.*dpO;

        % Conservation of water and oil
        water = (1/dt(n)).*(vol.*rW.*sW - vol0.*rW0.*sW0) + div(vW);
        oil   = (1/dt(n)).*(vol.*rO.*(1-sW) - vol0.*rO0.*(1-sW0)) + div(vO);

        % Injector: volumetric source term multiplied by surface density
        water(inIx) = water(inIx) - inRate.*rhoWS;

        % Producer: replace equations by new ones specifying fixed pressure
        % and zero water saturation
        oil(outIx) = p(outIx) - outPres;
        water(outIx)   = sW(outIx);

        % Collect and concatenate all equations (i.e., assemble and
        % linearize system)
        eqs = {oil, water};
        eq  = combineEquations(eqs{:});

        % Measure condition number
        cnd(i) = condest(eq.jac{1}); i = i+1;                              %#ok<SAGROW>

        % Compute Newton update and update variable
        res = eq.val;
        % Alternatively, you can use a LU factorization which may be
        % required for some Matlab versions when isCompr = false
        % [L, U] = lu(eq.jac{1});
        % upd = -(U\(L \ res));
        upd = -(eq.jac{1} \ res);
        p.val  = p.val   + upd(pIx);
        sW.val = sW.val + upd(sIx);
        sW.val = max( min(sW.val, 1), 0);

        resNorm = norm(res);
        nit     = nit + 1;
        fprintf('  Iteration %3d:  Res = %.4e\n', nit, resNorm);
    end
    % Add line with NaN in cnd variable to signify end of time step
    cnd(i) = NaN; i=i+1;

    if nit > maxits
        error('Newton solves did not converge')
    else % plot
        nits(n) = nit;

        figure(1); clf
        subplot(2,1,1)
        plotCellData(G, value(p)/barsa, pargs{:});
        title('Pressure'), view(30, 40); caxis([100 250]);

        subplot(2,1,2)
        plotCellData(G, 1-value(sW), pargs{:});
        caxis([0, 1]), view(30, 40); title('Saturation')
        drawnow

        sol(n+1) = struct('time', t, ...
                          'pressure', value(p), ...
                          's', value(sW));
        waitbar(t/Tf,hwb);
        pause(.1);
    end
end
close(hwb);

%%
figure(2),
subplot(2,1,2-double(isCompr));
bar(nits,1,'EdgeColor','r','FaceColor',[.6 .6 1]);
axis tight, set(gca,'YLim',[0 10]);
hold on,
plotyy(nan,nan,1:numel(dt),dt/day,@plot, ...
    @(x,y) plot(x,y,'-o','MarkerFaceColor',[.6 .6 .6]));

%%
figure(3); hold all
ipr = arrayfun(@(x) x.pressure(inIx), sol)/barsa;
t   = arrayfun(@(x) x.time, sol)/day;
fpr = arrayfun(@(x) sum(x.pressure)/G.cells.num, sol)/barsa;
plot(t,ipr,'-', t,fpr, '--');
% chl = get(gca,'Children');
% col=flipud(lines(2));
% for i=1:numel(chl), set(chl(i),'Color',col(ceil(i/2),:),'LineWidth',2); end

%%
figure(4); hold all
oip = arrayfun(@(x) sum(pv(x.pressure).*(1-x.s).*rhoO(x.pressure)), sol);
plot(t,(ioip - oip)./rhoOS,'-o','MarkerFaceColor',[.6 .6 .6])

%%
figure(5); hold all
plot(cumsum(isnan(cnd))+1, 1./cnd,'-o');

%%
%{
g = cartGrid([1 1 1],[Dx Dy Dz]);
figure; colormap(winter);
for i=1:nstep+1
    s = sol(i).s;
    clf
    plotCellData(G,1-s,s>1e-3,'EdgeColor','k','EdgeAlpha',.1);
    plotGrid(g,'EdgeColor','k','FaceColor','none');
    view(3); plotGrid(G,[inIx; outIx],'FaceColor',[.6 .6 .6]);
    view(30,40); caxis([0 1]); axis tight off
    axis([-.5 Dx+.5 -.5 Dy+.5 -.5 Dz+.5]);
    h=colorbar('horiz');
    pos = get(gca,'Position');
    set(h,'Position',[.13 .15 .775 .03]); set(gca,'Position',pos);
    title(['Time: ' num2str(t(i)) ' days']);
    drawnow, pause(.5)
    print('-dpng','-r0',sprintf('twoPhaseAD-%02d.png',i-1));
end
%}


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
