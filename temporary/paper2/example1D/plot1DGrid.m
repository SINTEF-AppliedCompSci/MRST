%G,Gt
%%
depth=1300;

[nx,ny,nz] = deal(100*n_fac, 1, 1);     % Cells in Cartsian grid
[Lx,Ly,H]  = deal(30e3,10e3,50); % Physical dimensions of reservoir
total_time = 5*year;             % Total simulation time
nsteps     = 40;                 % Number of time steps in simulation
dt         = total_time/nsteps;  % Time step length
perm       = 1000;                % Permeability in milli darcies
phi        = 0.03;                % Porosity
         % Initial depth


%% Create input deck and construct grid
% Create an input deck that can be used together with the fully-implicit
% solver from the 'ad-fi' module. Since the grid is constructed as part of
% setting up the input deck, we obtain it directly. 
G=cartGrid([nx,ny,nz],[Lx,Ly,H]);
x=G.nodes.coords(:,1);
LL=Lx*2/3;


G.nodes.coords(:,3)=G.nodes.coords(:,3)+depth-LL*sin(x/LL)*tan(phi)+2*sin(2*pi*x/0.3e3);
G=computeGeometry(G);
Gt=topSurfaceGrid(G);

clf,hold on
cartDims=Gt.cartDims;
nx=cartDims(1);
ind=1:nx+1;
x=Gt.nodes.coords(ind,1)/1e3;
z=Gt.nodes.z(ind);
x=[x;x(end:-1:1)];
z=[z;z(end:-1:1)+50];
patch(x,z,'y');set(gca,'YDir','reverse')
%set(gca,'FontSize',19)
%set(gca,'Ydir','reverse')
xlabel('X km','FontSize',16)
ylabel('Depth m','FontSize',16)
%set(gca,'FontSize',16)
h1 = gca;
pos= get(h1,'Position')
h3 = axes('Position',get(h1,'Position')*0.4+[0.5,0.13,0, 0]);
patch(x,z,'y');set(gca,'YDir','reverse');axis([5 6 1120 1205])
set(h3,'XTick',[5 6]);set(h3,'YTick',[1150 1175 1200])

axes(h1)
plot([5 6 6 5 5]',[1120 1120 1205 1205 1120],'k')

h4 = axes('Position',get(h1,'Position')*0.35+[0.11,0.56,0, 0]);
patch(x,z,'y');set(gca,'YDir','reverse');axis([25 26 720 733])
set(h4,'YAxisLocation','right') 
set(h4,'XTick',[25 26])

axes(h1)
plot([25 26 26 25 25]',[720 720 733 733 720],'k')
axes(h3)
axes(h4)
print('-depsc2','grid_1D example')