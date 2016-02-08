%% Advection 
% Initial data: set up the initial data and remember to add one ghost cell
% at each end of the interval. These cells will be used to impose boundary
% conditions
N  = 100;
dx = 1/N;
x  = -.5*dx:dx:1+.5*dx;
u0 = sin((x-.1)*pi/.3).^2.*double(x>=.1 & x<=.4);
u0((x<.9) & (x>.6)) = 1;

% Flux function
f = @(u) u;

% Run simulation
dt = 0.995*dx/1.0;
uu = upw(u0, dt, dx, 20, f, @periodic);
uf = lxf(u0, dt, dx, 20, f, @periodic);
uw = lxw(u0, dt, dx, 20, f, @periodic);

% Plot results
figure('Position',[293 539 1000 300]);
i = (x>=.1 & x<=.4);
subplot(1,3,1)
plot([0 x(i) .6 .6 .9 .9 1], [0 u0(i) 0 1 1 0 0], '-',...
    x(2:end-1), uu(2:end-1), 'o');
axis([0 1 -.2 1.2]); title('Upwind');

subplot(1,3,2)
plot([0 x(i) .6 .6 .9 .9 1], [0 u0(i) 0 1 1 0 0], '-',...
    x(2:end-1), uf(2:end-1), 's');
axis([0 1 -.2 1.2]); title('Lax-Friedrichs');

subplot(1,3,3)
plot([0 x(i) .6 .6 .9 .9 1], [0 u0(i) 0 1 1 0 0], '-',...
    x(2:end-1), uw(2:end-1), '*');
axis([0 1 -.2 1.2]); title('Lax-Wendroff');

%% Burgers' equation
N  = 50;
dx = 1/N;
x  = -1-.5*dx:dx:1+.5*dx;
u0 = sin(2*pi*x).*double(x>=-.5 & x<=.5);

% Flux function
f = @(u) .5*u.^2;

% Run simulation
dt = 0.995*dx/1.0;
[uf,uw] = deal(u0);
figure;
for i=1:12
    uf = lxf(uf, dt, dx, .05, f, @outflow);
    uw = lxw(uw, dt, dx, .05, f, @outflow);
    plot(x(2:end-1), uf(2:end-1), 's', ...
        x(2:end-1), uw(2:end-1), '*');
    axis([-1 1 -1.1 1.1]);
    drawnow; pause(.1);
end

%% Buckley--Leverett problem
% Flux function
f = @(u) u.^2./(u.^2 + (1-u).^2);
s = linspace(0,1,501);
df = max(diff(f(s))./diff(s));

% reference solution
N  = 1000;
dx = 1/N;
xr = -.5*dx:dx:1+.5*dx;
u0 = 0*xr; u0(1)=1.0;
ur = upw(u0, .995*dx/df, dx, .65, f, @inflow);

% Solutions on coarser grids
N  = 100;
dx = 1/N;
x  = -.5*dx:dx:1+.5*dx;
u0 = 0*x; u0(1)=1.0;
dt = .995*dx/df;

uu = upw(u0, dt, dx, .65, f, @inflow);
uf = lxf(u0, dt, dx, .65, f, @inflow);
uw = lxw(u0, dt, dx, .65, f, @inflow);

% Plot results
figure('Position',[293 539 1000 300]);
subplot(1,3,1)
plot(xr(2:end-1),ur(2:end-1), '-', x(2:end-1), uu(2:end-1), 'o'); 
axis([0 1 -.1 1.1]); title('Upwind');
subplot(1,3,2)
plot(xr(2:end-1),ur(2:end-1), '-', x(2:end-1), uf(2:end-1), 's');
axis([0 1 -.1 1.1]); title('Lax-Friedrichs');
subplot(1,3,3)
plot(xr(2:end-1),ur(2:end-1), '-', x(2:end-1), uw(2:end-1), '*');
axis([0 1 -.1 1.1]); title('Lax-Wendroff');


%% Buckley--Leverett problem using 'incomp' module
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
figure;
plot(G.cells.centroids(:,1), rSole.s(:,1),'bo', ...
    G.cells.centroids(:,1), rSoli.s(:,1), 'ks', ...
    G.cells.centroids(:,1), rSolt.s(:,1), 'r*', ...
    xr(2:end-1),ur(2:end-1),'-');
legend('Explicit','Implicit, CFL=1','Implicit, CFL=10',3);
axis([0 1 -.1 1.1]);
    