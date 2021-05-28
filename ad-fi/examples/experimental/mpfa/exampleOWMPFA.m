%% VE simulation in a standard black-oil solver
%  In this example we show how to set up a standard format black-oil
%  model that can be used to simulate a VE model. For the actual
%  simulation,  we use the fully-implicit solver in MRST from the 'ad-fi'
%  module, which is based on automatic differentiation. 
use_mpfa=true;
try
   require deckformat ad-fi
catch %#ok<CTCH>
   mrstModule add deckformat ad-fi
end

%% Parameters for the simulation
gravity off
[nx,ny,nz] = deal(100, 1, 1);     % Cells in Cartsian grid
[Lx,Ly,H]  = deal(400,10,15); % Physical dimensions of reservoir
total_time = year/25;             % Total simulation time
nsteps     = 25;                 % Number of time steps in simulation
dt         = total_time/nsteps;  % Time step length
perm       = 100;                % Permeability in milli darcies
phi        = 0.1;                % Porosity
depth      = 1000;               % Initial depth
ipress     = 200;                % Initial pressure

%% Create input deck and construct grid
% Create an input deck that can be used together with the fully-implicit
% solver from the 'ad-fi' module. Since the grid is constructed as part of
% setting up the input deck, we obtain it directly. 
fac=2;%inner_prod='ip_simple';
nx = fac*20+1; ny = fac*10; nz = 1;
nx_c=ceil(nx/2);
shift=floor(nx_c*2/3);
Lx=2000;Ly=1000;
G = cartGrid([nx, ny, nz],[2000 1000, 10]);
nodes_x=reshape(G.nodes.coords(:,1),G.cartDims+1);
nodes_y=reshape(G.nodes.coords(:,2),G.cartDims+1);
nodes_x_old=nodes_x;
ind_x=nx_c-shift +1:nx_c+shift -1;
ind_x=2:nx;
ind_y=1:ny-1;
nodes_x(ind_x,ind_y,:)=nodes_x_old(ind_x,ind_y,:)+(Ly-nodes_y(ind_x,ind_y,:))/20.*5.4.*sin(1*3.14*nodes_x(ind_x,ind_y,:)/Lx);
corr_cells=[nx_c,ny,nx_c,ny;
            nx_c-shift,1,nx_c-shift-fac*1,1;
            nx_c+shift,1,nx_c+shift-fac*2,1]
        
nw= size(corr_cells,1);       
wcells_ind=nan(nw,2);        
for j=1:size(corr_cells,1)
   wcells_new(j)=sub2ind(G.cartDims, corr_cells(j,3),corr_cells(j,4));
   wcells_ind(j,:)=[corr_cells(j,3),corr_cells(j,4)]
   nodes_x(corr_cells(j,3),corr_cells(j,4),:)=nodes_x_old(corr_cells(j,1),corr_cells(j,2),:)
   nodes_x(corr_cells(j,3)+1,corr_cells(j,4),:)=nodes_x_old(corr_cells(j,1)+1,corr_cells(j,2),:)
   nodes_x(corr_cells(j,3)+1,corr_cells(j,4)+1,:)=nodes_x_old(corr_cells(j,1)+1,corr_cells(j,2)+1,:)
   nodes_x(corr_cells(j,3),corr_cells(j,4)+1,:)=nodes_x_old(corr_cells(j,1),corr_cells(j,2)+1,:)    
end

G.nodes.coords(:,1)=nodes_x(:);
G = computeGeometry(G);
rock.perm  = repmat(100*milli*darcy, [G.cells.num, 1]);
rock.poro  = repmat(0.3, [G.cells.num, 1]);
clf,plotGrid(G)


%% Introduce wells
% <html>
% We will include two wells, one rate-controlled vertical well and one
% horizontal well controlled by bottom-hole pressure. Wells are described
% using a Peacemann model, giving an extra set of equations that need to be
% assembled, see <a href="../../1ph/html/simpleWellExample.html#3"> "Using
% Peacemann well models"</a> for more details.
% </html>

inner_prod='ip_tpf'
W = verticalWell([], G, rock, wcells_ind(1,1) , wcells_ind(1,2),1:G.cartDims(3),...          ...
            'Type', 'rate', 'Val', 1e4/day(), ...
            'Radius', 0.1, 'Comp_i', [1, 0, 0],'InnerProduct',inner_prod);
W = verticalWell(W, G, rock, wcells_ind(2,1) , wcells_ind(2,2), 1:G.cartDims(3),...
            'Type', 'bhp' , 'Val', 1.0e5, ...
            'Radius', 0.1, 'Dir', 'y', 'Comp_i', [0, 1, 0],'InnerProduct',inner_prod,'Sign'        , -1);
W = verticalWell(W, G, rock, wcells_ind(3,1) , wcells_ind(3,2), 1:G.cartDims(3),...
            'Type', 'bhp' , 'Val', 1.0e5, ...
            'Radius', 0.1, 'Dir', 'y', 'Comp_i', [0, 1, 0],'InnerProduct',inner_prod, 'Sign'        , -1);
Wext=W;
figure(30)
plotGrid(G)
plotWell(G,W)
plotGrid(G,vertcat(W.cells),'FaceColor','r');

%figure(, plotGrid(G),view([0 -1 0]), box on
                                      
%% Initialize data structures
% First, we convert the input deck to SI units, which is the unit system
% used by MRST. Second, we initialize the rock parameters from the deck;
% the resulting data structure may have to be post-processed to remove
% inactive cells. Then we set up the fluid object and tell the ad-fi solver
% that that we are working with an oil-gas system.
fluid = initSimpleADIFluid('mu' , [   1,  10, 1] .* centi*poise     , ...
                        'rho', [1000, 700, 700] .* kilogram/meter^3, ...
                        'n'  , [   2,   2, 2]);
%[krW, krO] = f.relPerm(sW);
fluid.relPerm =@(sw) deal(fluid.krW(sw),fluid.krO(1-sw));

systemOW  = initADISystem({'Oil', 'Water'}, G, rock, fluid);
s_org=systemOW.s;
s=s_org;
if(use_mpfa)
   %% TPF
    s.grad=@(x) s_org.T.*s_org.grad(x);
    s.T=ones(size(s.T));
    systemOW.s=s;
else
   %% MPFA
    mrstModule add mpfa
    nc=G.cells.num;
    % matrix from cells+boundary to hface, where all value on a face is
    % equal (??)
    TM=computeMultiPointTrans(G, rock);
    T=TM.T;
    cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
    b  = any(G.faces.neighbors==0, 2);
    I1 = [(1:nc)'; G.cells.num + find(b)];
    C =sparse(1:nc,1:nc,1,nc,nc+sum(b));
    C= C+sparse(double(sum(G.faces.neighbors(b,:),2)), [nc+1:nc+sum(b)]', 1,nc,nc+sum(b));
    ph=G.faces.neighbors(G.cells.faces(:,1),1)==cellNo;
    THF=T(:,I1)*C';
    clear TF;
    TF(G.cells.faces(ph),:)=THF(ph,:);
    % grad fpr internal faces
    grad=TF(~b,:);
    %
    s.T=ones(size(s.T));
    s.grad=@(x) grad*x;
    systemOW.s=s;
end


% calculate conducivity for fock
systemOW.getEquations =@ eqsfiOWExplicitWells;
systemOW.nonlinear.linesearch=true;

%% Run the schedule setup in the file
% Before we can run the schedule, we make sure that we have an initial
% hydrostatic pressure distribution. Then we pick the schedule from the
% input deck and start the simulator.
p0=ones(G.cells.num,1)*300*barsa;
s0=repmat([0 1],G.cells.num);
x0= struct('pressure',p0,'s',s0);
x0.wellSol=initWellSolLocalOrig(W, x0);
x0.wellSol(1).bhp=300;
%x0 = initEclipseState(G, deck, initEclipseFluid(deck));

%Wext=processWellsLocal(G, rock, deck.SCHEDULE.control(1));
%for i=1:numel(Wext)
%  Wext(i).T=273;
%  
%  Wext(i).I=ones(1,size(x0.I,2));
%end

T      = 400*day();
dT     = T/20;
dTplot = 100*day();  % plot only every 100th day
N      = fix(T/dTplot);
N_step = ceil(T/dT); 
mrst_schedule.step=struct('control',ones(N_step,1),'val',dT*ones(N_step,1));
mrst_schedule.W={W};

[wellSols, states] = runMrstADI(x0, G, systemOW, mrst_schedule,'Wext',Wext);

%%
TT=cumsum(ones(N_step,1)*dt)/day;
qOs=[];qWs=[];
for i=1:numel(wellSols)
 qOs=[qOs;wellSols{i}.qOs];
 qWs=[qWs;wellSols{i}.qWs];
end
wcut=qWs./(qOs+qWs);
figure(31)
plot(TT,wcut)
xlabel('days')
ylabel('water cut')

%%
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
