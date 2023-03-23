%% Advection: classic schemes 
% Initial data: set up the initial data and remember to add one ghost cell
% at each end of the interval. These cells will be used to impose boundary
% conditions
FontSize = 12; pargs = {'MarkerSize',8};
dx = 1/100;
x  = -.5*dx:dx:1+.5*dx;
u0 = sin((x-.1)*pi/.3).^2.*double(x>=.1 & x<=.4);
u0((x<.9) & (x>.6)) = 1;

% Flux function
f  = @(u) u;
df = @(u) 0*u + 1;

% Run simulation
uu = upw(u0, .5, dx, 1, f, df, @periodic);
uf = lxf(u0, .5, dx, 1, f, df, @periodic);
uw = lxw(u0, .5, dx, 1, f, df, @periodic);

% Plot results
figure('Position',[293 539 1000 300]);
i = (x>=.1 & x<=.4);
subplot(1,4,1,'FontSize',FontSize)
plot([0 x(i) .6 .6 .9 .9 1], [0 u0(i) 0 1 1 0 0], '-k',...
    x(2:end-1), uu(2:end-1), '.', pargs{:});
axis([0 1 -.2 1.2]); title('Upwind');

subplot(1,4,2,'FontSize',FontSize)
plot([0 x(i) .6 .6 .9 .9 1], [0 u0(i) 0 1 1 0 0], '-k',...
    x(2:end-1), uf(2:end-1), '.', pargs{:});
axis([0 1 -.2 1.2]); title('Lax-Friedrichs');

subplot(1,4,3,'FontSize',FontSize)
plot([0 x(i) .6 .6 .9 .9 1], [0 u0(i) 0 1 1 0 0], '-k',...
    x(2:end-1), uw(2:end-1), '.', pargs{:});
axis([0 1 -.2 1.2]); title('Lax-Wendroff');


%% Advection: high-resolution schemes
lim = @(r) max(0,min(2*r,min(.5+.5*r,2)));
xh  = -1.5*dx:dx:1+1.5*dx;
u0 = sin((xh-.1)*pi/.3).^2.*double(xh>=.1 & xh<=.4);
u0((xh<.9) & (xh>.6)) = 1;

un = nt (u0, .25, dx, 1, f, df, @periodic, lim);
%uc = cuw(u0, .25, dx, 1, f, df, @periodic, lim);

subplot(1,4,4,'FontSize',FontSize)
i = (xh>=.1 & xh<=.4);
plot([0 x(i) .6 .6 .9 .9 1], [0 u0(i) 0 1 1 0 0], '-k', ...
    xh(3:end-2), un(3:end-2), '.', pargs{:});
title('Nessyahu-Tadmor');
axis([0 1 -.2 1.2]);

%% Burgers' equation
f = @(u) .5*u.^2;
df = @(u) u;

% Run simulations
dx = 1/20;
x  = -1-0.5*dx:dx:1+0.5*dx;
uf = sin(2*pi*x).*double(x>=-.5 & x<=.5);
uw = lxw (uf, .995, dx, .6, f, df, @outflow);
uf = lxf (uf, .995, dx, .6, f, df, @outflow);

xh = -1-1.5*dx:dx:1+1.5*dx;
uh = sin(2*pi*xh).*double(xh>=-.5 & xh<=.5);
uh = nt(uh, .495, dx, .6, f, df, @outflow, lim);

dx = 1/1000;
xr = -1-1.5*dx:dx:1+1.5*dx;
ur = sin(2*pi*xr).*double(xr>=-.5 & xr<=.5);
ur = cuw(ur, .495, dx, .6, f, df, @outflow);

% Plot results
figure('Position',[293 539 1000 300]);
subplot(1,3,1)
plot(xr(3:end-2),ur(3:end-2), '-', x(2:end-1), uf(2:end-1), 'o'); 
axis([-1 1 -1.1 1.1]); title('Lax--Friedrichs');
subplot(1,3,2)
plot(xr(3:end-2),ur(3:end-2), '-', x(2:end-1), uw(2:end-1), 'o');
axis([-1 1 -1.1 1.1]); title('Lax-Wendroff');
subplot(1,3,3)
plot(xr(3:end-2),ur(3:end-2), '-', xh(3:end-2), uh(3:end-2), 'o');
axis([-1 1 -1.1 1.1]); title('Nessyahu-Tadmor');


%% Buckley-Leverett problem: classical schemes vs high-resolution
% Flux function
f = @(u) u.^2./(u.^2 + (1-u).^2);
s = linspace(0,1,501);
dfv = max(diff(f(s))./diff(s));
df = @(u) 0*u + dfv;

% reference solution
T = 0.65;
dx = 1/1000;
xr = -.5*dx:dx:1+.5*dx;
u0 = 0*xr; u0(xr<.1)=1.0;
ur = upw(u0, .995, dx, T, f, df, @outflow);

% Solutions on coarser grids
N  = 50;
dx = 1/N;
x  = -.5*dx:dx:1+.5*dx;
u0 = 0*x; u0(x<.1)=1.0;

figure('Position',[293 539 1000 300]);
subplot(1,4,1,'FontSize',FontSize)
uu = upw(u0, .995, dx, T, f, df, @outflow);
plot(xr(2:end-1),ur(2:end-1), '-k', x(2:end-1), uu(2:end-1), '.', pargs{:}); 
axis([0 1 -.05 1.05]); title('Upwind');

subplot(1,4,2,'FontSize',FontSize)
uf = lxf(u0, .995, dx, T, f, df, @outflow);
plot(xr(2:end-1),ur(2:end-1), '-k', x(2:end-1), uf(2:end-1), '.', pargs{:});
axis([0 1 -.05 1.05]); title('Lax-Friedrichs');

subplot(1,4,3,'FontSize',FontSize)
uw = lxw(u0, .995, dx, T, f, df, @outflow);
plot(xr(2:end-1),ur(2:end-1), '-k', x(2:end-1), uw(2:end-1), '.', pargs{:});
axis([0 1 -.05 1.05]); title('Lax-Wendroff');

% High-resolution schemes
xh  = -1.5*dx:dx:1+1.5*dx;
u0 = 0*xh; u0(xh<.1)=1.0;
un = nt(u0, .495, dx, T, f, df, @inflow);
uc = cuw(u0, .495, dx, T, f, df, @inflow);
subplot(1,4,4,'FontSize',FontSize)
plot(xh(3:end-2), un(3:end-2), '.', xh(3:end-2), uc(3:end-2), '.r', ...
    xr(2:end-1),ur(2:end-1), 'k-', pargs{:});
axis([0 1 -.05 1.05]); title('High-resolution');

%% Buckley-Leverett problem: using 'incomp' module
mrstModule add incomp
G = computeGeometry(cartGrid([100,1],[1 1]));
rock = makeRock(G, ones(G.cells.num,1), ones(G.cells.num,1));

fluid = initSimpleFluid('mu' , [   1,    1] .* centi*poise     , ...
                        'rho', [1000, 1000] .* kilogram/meter^3, ...
                        'n'  , [   2,    2]);
bc = fluxside([], G, 'Left',   1, 'sat', [1 0]);
bc = fluxside(bc, G, 'Right', -1, 'sat', [0 1]);

hT = computeTrans(G, rock);
rSol = initState(G, [], 0, [0 1]);

rSol = incompTPFA(rSol, G, hT, fluid, 'bc', bc);
rSole = explicitTransport(rSol, G, .65, rock, fluid, 'bc', bc);
rSoli = rSol; n = 13;
for i=1:n
    rSoli = implicitTransport(rSoli, G, .65/n, rock, fluid, 'bc', bc);
end
rSolt = rSol; n = 130;
for i=1:n
    rSolt = implicitTransport(rSolt, G, .65/n, rock, fluid, 'bc', bc);
end

dx = 1/1000;
xr = -.5*dx:dx:1+.5*dx;
u0 = 0*xr; u0(1)=1.0;
ur = upw(u0, .995, dx, .65, f, df, @outflow);

figure;
plot(G.cells.centroids(:,1), rSole.s(:,1),'o', ...
    G.cells.centroids(:,1), rSolt.s(:,1), 's', ...
    G.cells.centroids(:,1), rSoli.s(:,1), '*', ...
    xr(2:end-1),ur(2:end-1),'-k');
legend('Explicit','Implicit, CFL=1','Implicit, CFL=10','Location','SouthWest');
