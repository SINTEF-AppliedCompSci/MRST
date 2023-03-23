function sec105()
%% Single-phase compressible AD solver with thermal effects
% The purpose of the example is to give demonstrate rapid prototying in
% MRST using the automatic differentiation (AD) class. To this end, we
% extend the <singlePhaseAD.m flow simulator for compressible single-phase
% flow> to include temperature effects, which are modeled by incorporating
% a second conservation equation for energy. Except for this, the
% computational setup is the same with a single horizontal well draining
% fluids from a simple box-geometry reservoir.

mrstModule add incomp

%% Set up model: grid, perm and poro
[nx,ny,nz] = deal( 10,  10, 10);
[Dx,Dy,Dz] = deal(200, 200, 50);
G = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
G = computeGeometry(G);

%% Define rock model
rock.perm = repmat(30*milli*darcy, [G.cells.num, 1]);
rock.poro = repmat(0.3, [G.cells.num, 1]);
model.cr   = 1e-6/barsa;
model.p_r  = 200*barsa;
model.pv_r = poreVolume(G, rock);

%% Define model for compressible fluid
model.mu0   = 5*centi*poise;
model.cmup  = 2e-3/barsa;
model.cmut  = 1e-3;
model.T_r   = 300;

model.alpha = 5e-3;
model.beta  = 1e-3/barsa;
model.rho_r = 850*kilogram/meter^3;
model.rho_S = 750*kilogram/meter^3;
rho = @(p,T) ...
   model.rho_r .* (1+model.beta*(p-model.p_r)).*exp(-model.alpha*(T-model.T_r) );

%% Quantities for energy equation
model.Cp = 4e3;

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
equil  = ode23(@(z,p) g.* rho(p,model.T_r), [z_0, z_max], model.p_r);
p_init = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil
T_init = ones(G.cells.num,1)*model.T_r;

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
model.Tp = Tp;

kap = 4;                                      % Heat conduction of granite
tmp = struct('perm',kap*ones(G.cells.num,1)); % Temporary rock object
hT  = computeTrans(G, tmp);                   % Half-transmissibilities
Th  = 1 ./ accumarray(cf, 1 ./ hT, [nf, 1]);  % Harmonic average
Th  = Th(intInx);                             % Restricted to interior
model.Th = Th;

%% Define discrete operators
n    = size(N,1);
C    = sparse( [(1:n)'; (1:n)'], N, ones(n,1)*[-1 1], n, G.cells.num);
model.grad = @(x) C*x;
model.div  = @(x) -C'*x;
model.avg  = @(x) 0.5 * (x(N(:,1)) + x(N(:,2)));
model.upw  = @(x,flag) x(N(:,1)).*double(flag)+x(N(:,2)).*double(~flag);
model.gdz  = model.grad(G.cells.centroids)*gravity()';

%% Initialize for solution loop
state = vertcat(p_init,T_init, p_init(W.cells(1)), 0);
nc = G.cells.num;
[model.pIx, model.TIx, model.bhpIx, model.qSIx] = ....
   deal(1:nc, nc+1:2*nc, 2*nc+1, 2*nc+2);

numSteps = 78;                  % number of time-steps
totTime  = 365*day*1.5;         % total simulation time
model.dt = totTime / numSteps;  % constant time step
tol      = 1e-5;                % Newton tolerance
maxits   = 10;                  % max number of Newton its

sol = repmat(struct('time',[], 'pressure',[], 'bhp',[], 'qS',[], ...
    'T',[], 'qH',[]),[numSteps+1,1]);
sol(1)  = struct('time', 0, 'pressure', state(model.pIx), ...
    'bhp', state(model.bhpIx), 'qS', state(model.qSIx), ...
    'T', state(model.TIx),'qH', 0);

%% Main loop
t = 0; step = 0;
while t < totTime
   t = t + model.dt;
   step = step + 1;
   fprintf('\nTime step %d: Time %.2f -> %.2f days\n', ...
      step, convertTo(t - model.dt, day), convertTo(t, day));

   % Newton loop
   resNorm = 1e99;
   state0  = state;
   nit = 0;
   while (resNorm > tol) && (nit < maxits)

      eq = eqsfiThermal(state0, state, G, W, model);
      J   = eq.jac{1};  % Jacobian
      res = eq.val;     % residual
      upd = -(J \ res); % Newton update
      state = state + upd;
      resNorm = norm(res);
      nit     = nit + 1;
      fprintf('  Iteration %3d:  Res = %.4e\n', nit, resNorm);
   end

   if nit > maxits
      error('Newton solves did not converge')
   else % store solution
      sol(step+1)  = struct('time', t, 'pressure', state(model.pIx), ...
                            'bhp', state(model.bhpIx), ...
                            'qS', state(model.qSIx), ...
                            'T', state(model.TIx), 'qH', []);
   end
end
%% Plot production rate
clf, set(gca,'FontSize',20);
stairs([sol(2:end).time]/day,-[sol(2:end).qS]*day,'LineWidth',2);
xlabel('days');ylabel('m^3/day')

%% Plot pressure evolution
dt = model.dt;
figure; clf;
steps = [2 5 10 25];
p = vertcat(sol(:).pressure);
cax=[min(p) max(p)]./barsa;
for i=1:4
   subplot(2,2,i);
   plotCellData(G, sol(steps(i)).pressure/barsa, show, 'EdgeColor',.5*[1 1 1]);
   plotWell(G, W, 'FontSize',12);
   view(-125,20), camproj perspective
   caxis(cax);
   axis tight off; zoom(1.2), set(gca,'Clipping','off')
   text(200,170,-8,[num2str(round(steps(i)*model.dt/day)) ' days'],'FontSize',12);
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
   plotCellData(G, sol(steps(i)).T, show, 'EdgeColor',.5*[1 1 1]);
   plotWell(G, W, 'FontSize',12);
   view(-125,20), camproj perspective
   caxis(cax);
   axis tight off; zoom(1.2), set(gca,'Clipping','off')
   text(200,170,-8,[num2str(round(steps(i)*dt/day)) ' days'],'FontSize',12);
end
colorbar('horiz','Position',[.1 .05 .8 .025]);
colormap(jet(55));
%%
wc = W(1).cells;
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
figure; plot(t,Pm,t,PM,t,Pa,t,Pw,'.k','LineWidth',2);
legend('min(p)', 'max(p)', 'avg(p)', 'wells');
figure; plot(t,Tm,t,TM,t,Ta,t,Tw,'.k','LineWidth',2);
legend('min(T)', 'max(T)', 'avg(T)', 'wells');

end

function [eqs] = eqsfiThermal(s0, s, G, W, m)

[p,T,bhp,qS] = initVariablesADI(s(m.pIx), s(m.TIx), s(m.bhpIx), s(m.qSIx));
p0 = s0(m.pIx);
T0 = s0(m.TIx);

pv   = m.pv_r .* exp( m.cr * (p - m.p_r) );
pv0  = m.pv_r .* exp( m.cr * (p0 - m.p_r) );
sv   = G.cells.volumes - pv;
sv0  = G.cells.volumes - pv0;
pf   = m.avg(p);
Tf   = m.avg(T);
muf  = m.mu0*(1+m.cmup*(pf-m.p_r)).*exp(-m.cmut*(Tf-m.T_r));
rho  = m.rho_r .* (1+m.beta*(p -m.p_r)).*exp(-m.alpha*(T -m.T_r));
rho0 = m.rho_r .* (1+m.beta*(p0-m.p_r)).*exp(-m.alpha*(T0-m.T_r));
rhof = m.avg(rho);
Hf   =  m.Cp*T  + (1-m.T_r*m.alpha).*(p -m.p_r)./rho;
Hf0  =  m.Cp*T0 + (1-m.T_r*m.alpha).*(p0-m.p_r)./rho0;
Er   =  m.Cp*T;
Er0  =  m.Cp*T0;
Efrho  =  Hf.*rho  - p;
Efrho0 =  Hf0.*rho0 - p0;

v   = -(m.Tp./muf).*(m.grad(p) - rhof.*m.gdz);
pEq = (1/m.dt)*(pv.*rho - pv0.*rho0) + m.div(rhof.*v);
hEq = (1/m.dt)*(pv.*Efrho + sv.*Er - pv0.*Efrho0 - sv0.*Er0) ...
     + m.div( m.upw(Hf,v>0).*rhof.*v) + m.div( -m.Th.*m.grad(T));

wc = W(1).cells;
assert(numel(wc)==numel(unique(wc)));
WI  = W(1).WI;
dz  = W(1).dZ;
bhT = ones(size(wc))*200;
g = norm(gravity);

rhob   = m.rho_r .* (1+m.beta*(bhp -m.p_r)).*exp(-m.alpha*(bhT -m.T_r));
Hfb    = m.Cp*bhT + (1-m.T_r*m.alpha).*(bhp -m.p_r)./rhob;
muW    = m.mu0*(1+m.cmup*(p(wc)-m.p_r)).*exp(-m.cmut*(T(wc)-m.T_r));
p_conn = bhp + g*dz.*rhob;
q_conn = WI .* (rho(wc)./muW).*(p_conn-p(wc));

pEq(wc)  = pEq(wc) - q_conn;
qw       = q_conn;
hq       = Hfb.*qw;
hq(qw<0) = Hf(wc(qw<0)).*qw(qw<0);
hEq(wc)  = hEq(wc) - hq;

rateEq = qS - sum(q_conn)/m.rho_S;
ctrlEq = bhp-100*barsa;

eqs = {pEq, hEq, rateEq, ctrlEq};
eqs = cat(eqs{:});
end