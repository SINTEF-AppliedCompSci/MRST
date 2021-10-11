%% LOAD NECESSARY MODULES
mrstModule add hfm;             % hybrid fracture module
mrstModule add mrst-gui;        % plotting routines
mrstModule add ad-blackoil;     % AD blackoil solver
mrstModule add ad-core;         % AD core module
mrstModule add ad-props;        % AD properties

%% SET UP A STRUCTURED MATRIX GRID
physdim = [350, 200, 100]; % 350m x 200m x 100m domain
celldim = [35, 20, 10]; % 10m x 10m x 10m grid cell sizes
G = cartGrid(celldim, physdim);
G = computeGeometry(G);
G.rock=makeRock(G,100*milli*darcy,0.3); % km=100mD, matrix porosity = 0.3

%% SET UP FRACTURE 1 (RED)
fracplanes(1).points = [40 100 0;
    90 160 0;
    90 160 100;
    40 100 100]; % Vertices
fracplanes(1).aperture = 1/25;
fracplanes(1).poro=0.8;
fracplanes(1).perm=10000*darcy;

figure;
plotfracongrid(cartGrid([1 1 1],physdim),fracplanes,'label',false,...
    'fracfacealpha',1,'fracfacecolor','r','fracturelist',1);
view(30,45);
axis equal tight

%% SET UP FRACTURE 2 (GREEN)
points = [80 160 0;
    290 40 0;
    290 40 100;
    80 160 100]; %vertices

f2normal = getnormal(points);
points([1,2],:)=points([1,2],:)-f2normal*15; % displace top points
points([3,4],:)=points([3,4],:)+f2normal*15; % displace bottom points

fracplanes(2).points = points;
fracplanes(2).aperture = 1/25;
fracplanes(2).poro=0.8;
fracplanes(2).perm=10000*darcy;

figure;
plotfracongrid(cartGrid([1 1 1],physdim),fracplanes,'label',false,...
    'fracfacealpha',1,'fracfacecolor','g','fracturelist',2);
view(30,45);
axis equal tight

%% SET UP FRACTURE 3 (YELLOW)
fracplanes(3).points = [200 70 0;
    280 160 0;
    280 160 100;
    200 70 100]; % Vertices
fracplanes(3).aperture = 1/25;
fracplanes(3).poro=0.8;
fracplanes(3).perm=10000*darcy;

figure;
plotfracongrid(cartGrid([1 1 1],physdim),fracplanes,'label',false,...
    'fracfacealpha',1,'fracfacecolor','y','fracturelist',3);
view(30,45);
axis equal tight

%% VISUALIZE FRACTURES
figure;
plotfracongrid(G,fracplanes,'label',false);
view(30,45);
axis equal tight

%% GENERATE FRACTURE GRID
tol=1e-5;
Nm = G.cells.num;

for i=1:length(fracplanes)
    points = fracplanes(i).points;
    aperture = fracplanes(i).aperture;
    
    % Calculate plane unit normal
    diffp = diff(points,1);
    planenormal = cross(diffp(1,:), diffp(2,:));
    planenormal = planenormal/norm(planenormal);
    
    % Instantiate data types to hold fracture grid information
    fraccellpoints=cell(Nm,1); % vertices of each fracture grid cell
    area=-1*ones(Nm,1); % area of each fracture grid cell
    
    % Calculate intersection of fracture with each matrix grid cell
    for j=1:Nm
        [cn,cpos]=gridCellNodes(G,j);
        cellnodecoords = G.nodes.coords(cn,:);
        [~,area(j),~,~,fraccellpoints{j}]=...
            pebiAABBintersect(points,cellnodecoords,tol);
    end
    
    % Consolidate intersection data
    intersected=~cellfun('isempty',fraccellpoints);
    fraccellpoints=fraccellpoints(intersected);
    mcells=find(intersected);
    area=area(intersected);
    
    % Generate grid using intersection data
    V=vertcat(fraccellpoints{:}); % Coordinates of all vertices
    C=cellfun(@(c) 1:size(c,1),...
        fraccellpoints,'UniformOutput',false); % Indirection map
    for j=2:size(C,1)
        addTo=C{j-1}(end);
        C{j}=C{j}+addTo;
    end
    Fgrid(i).grid=fractureplanegeneralgrid(V,C,points,...
        planenormal,aperture,tol);
    Fgrid(i).matrix_connection.cells=mcells;
    Fgrid(i).matrix_connection.area=area;
end

%% PLOT FRACTURES
figure;
hold on;
colors = ['r','g','y'];
for i=1:3
    plotGrid(Fgrid(i).grid,'FaceColor',colors(i));
end
axis equal tight
view(30,45);
xlim([0 physdim(1)]);
ylim([0 physdim(2)]);

%% PLOT INTERSECTED MATRIX GRIDS
figure;
hold on;
colors = ['r','g','y'];
for i=1:3
    plotGrid(G,Fgrid(i).matrix_connection.cells,...
        'FaceAlpha', 0.5, 'FaceColor', colors(i));
end
axis equal tight
view(30,45);
xlim([0 physdim(1)]);
ylim([0 physdim(2)]);


%% ASSEMBLE GLOBAL GRID
% Initiate starting indices for cells, faces and nodes
cstart = G.cells.num+1;
fstart = G.faces.num+1;
nstart = G.nodes.num+1;

% Append Fgrid to G
for i=1:length(fracplanes)
    fieldname=['Frac',num2str(i)];
    
    % Add fracture grid to G.FracGrid
    G.FracGrid.(fieldname)=Fgrid(i).grid;
    
    % Save global starting indices and compute next one
    G.FracGrid.(fieldname).cells.start = cstart;
    G.FracGrid.(fieldname).faces.start = fstart;
    G.FracGrid.(fieldname).nodes.start = nstart;
    cstart = cstart + G.FracGrid.(fieldname).cells.num;
    fstart = fstart + G.FracGrid.(fieldname).faces.num;
    nstart = nstart + G.FracGrid.(fieldname).nodes.num;
    
    % Append poroperm data
    G.FracGrid.(fieldname).rock.perm = ones(G.FracGrid.(fieldname).cells.num,1)*fracplanes(i).perm;
    G.FracGrid.(fieldname).rock.poro = ones(G.FracGrid.(fieldname).cells.num,1)*fracplanes(i).poro;
end
G.nnc=[]; % Instantiate an empty list of NNCs
G=assembleGlobalGrid(G); % Create a global grid

%% GENERATE FRACTURE MATRIX NNC
G.nnc.cells=[];
G.nnc.T=[];
G.nnc.area=[];
G.nnc.type=[];
tol=1e-6;

for i=1:length(fracplanes)
    points = fracplanes(i).points;
    diffp = diff(points,1);
    planenormal = cross(diffp(1,:), diffp(2,:));
    planenormal = planenormal/norm(planenormal);
    fieldname=['Frac',num2str(i)];
    
    % Generate list of NNC pairs
    fcellstart=G.FracGrid.(fieldname).cells.start;
    Nf=Fgrid(i).grid.cells.num;
    mcells=Fgrid(i).matrix_connection.cells;
    N_nnc = length(mcells);
    nncpairs=[mcells,repmat(((1:Nf)+fcellstart-1)',N_nnc/Nf,1)];
    
    CI=ones(N_nnc,1); % instantiate list of CI
    
    for j=1:N_nnc
        [cn,cpos]=gridCellNodes(G,Fgrid(i).matrix_connection.cells(j));
        cellnodecoords = G.nodes.coords(cn,:);
        
        % Calculate average distance
        davg=calcdavg(cellnodecoords,planenormal,points(1,:),tol);
        
        % Calculate CI
        CI(i)=calcfracmatCI(cellnodecoords,...
            Fgrid(i).matrix_connection.area(j),...
            planenormal,points(1,:),davg,tol);
    end
    
    % Calculate cell to cell transmissibility
    pv = poreVolume(G,G.rock);
    w1 = pv(nncpairs(:,1))./G.rock.perm(nncpairs(:,1));
    w2 = pv(nncpairs(:,2))./G.rock.perm(nncpairs(:,2));
    wt = pv(nncpairs(:,1))+pv(nncpairs(:,2));
    Trans = 2*CI.*(wt./(w1+w2));
    
    % Append NNC data
    G.nnc.cells=[G.nnc.cells;nncpairs];
    G.nnc.T=[G.nnc.T;Trans];
    G.nnc.area=[G.nnc.area;Fgrid(i).matrix_connection.area];
end

%% GENERATE FRACTURE-FRACTURE NNC
tol=1e-6;
nnc_pairs=[];
nnc_T=[];

% Indices of matrix cells which contain fracture cells
mcells = unique(G.nnc.cells(:,1));
N_nnc_max = length(mcells);

for i=1:N_nnc_max
    mcelli=mcells(i);
    ind=ismember(G.nnc.cells(:,1),mcelli);
    
    if sum(ind)==2        
        fraccells = G.nnc.cells(ind,2)';
        fraccellareas = G.nnc.area(ind);
        
        cellnodes1 = G.nodes.coords(gridCellNodes(G,fraccells(1)),:);
        cellnodes2 = G.nodes.coords(gridCellNodes(G,fraccells(2)),:);
        
        % Intersection calculation
        [intersected,~,xlength,df]=...
            PEBIPEBIintersection3D(cellnodes1,cellnodes2,tol);
        
        if intersected
            nnc_pairs=[nnc_pairs; sort(fraccells)];
            
            % Calculate transmissibility
            aperture1 = G.cells.volumes(fraccells(1))/fraccellareas(1);
            aperture2 = G.cells.volumes(fraccells(2))/fraccellareas(2);
            Trans1=G.rock.perm(fraccells(1))*aperture1*xlength/df(1);
            Trans2=G.rock.perm(fraccells(2))*aperture2*xlength/df(2);
            Trans=(Trans1*Trans2)/(Trans1+Trans2);
            nnc_T=[nnc_T; Trans];
        end
    end  
end

G.nnc.cells=[G.nnc.cells;nnc_pairs];
G.nnc.T=[G.nnc.T;nnc_T];

%% PLOT INTERSECTING PAIRS
figure;
hold on
axis equal tight
xlim([0 physdim(1)])
ylim([0 physdim(2)])
colors=['r','g','y'];
allnnccells=nnc_pairs(:);
for i=1:3
    fieldname=['Frac',num2str(i)];
    cstart=G.FracGrid.(fieldname).cells.start;
    cend=cstart-1+G.FracGrid.(fieldname).cells.num;
    
    plotGrid(G,cstart:cend,'FaceAlpha',0);
    plotGrid(G,allnnccells(ismember(allnnccells',cstart:cend)),'FaceColor',colors(i));
    
end
view(30,45)

%% PLOT FRACTURE GRID
figure;
plotGrid(cartGrid([1 1 1],physdim),'facealpha',0);
hold on;
plotGrid(G,G.FracGrid.Frac1.cells.start-1+(1:G.FracGrid.Frac1.cells.num),'facecolor','r');
plotGrid(G,G.FracGrid.Frac2.cells.start-1+(1:G.FracGrid.Frac2.cells.num),'facecolor','g');
plotGrid(G,G.FracGrid.Frac3.cells.start-1+(1:G.FracGrid.Frac3.cells.num),'facecolor','y');
axis tight equal
title('Fracture Grid')
view(30,45);

%% FLUID PROPERTIES
% Define a three-phase fluid model without capillarity. Properties are
% listed in the order 'Water-Oil-Gas'.
pRef = 100*barsa;
fluid = initSimpleADIFluid('phases' , 'WOG', ...
    'mu' , [   1,  5, 0.2] * centi*poise     , ...
    'rho', [1000, 700, 250] * kilogram/meter^3, ...
    'c',   [1e-8, 1e-5, 1e-3] / barsa, ...
    'n'  , [   2,   2, 2], ...
    'pRef' , pRef);

%% DEFINE BLACK OIL MODEL
% Define a black oil model without dissolved gas or vaporized oil. Gravity
% is disabled.
gravity off
model = ThreePhaseBlackOilModel(G, G.rock, fluid, ...
    'disgas', false, 'vapoil', false);
model.operators = setupEDFMOperatorsTPFA(G, G.rock, tol);

%% ADD INJECTOR
totTime = 5*year;
tpv = sum(model.operators.pv);
wellRadius = 0.1;
[nx, ny, nz] = deal(G.cartDims(1), G.cartDims(2), G.cartDims(3));
cellinj = 1:nx*ny:(1+(nz-1)*nx*ny);
W   = addWell([], G, G.rock, cellinj, 'Type', 'rate', ...
    'Val', tpv/totTime, 'Radius', wellRadius, ...
    'Comp_i', [1, 0, 0], 'Name', 'Injector');

%% ADD PRODUCER
cellprod = nx*ny : nx*ny : nz*nx*ny;
W   = addWell(W, G, G.rock, cellprod, 'Type', 'bhp', ...
    'Val', 50*barsa, 'Radius', wellRadius, ...
    'Comp_i', [1, 1, 0], 'Name', 'Producer');

%% INITIALIZATION
% At time zero, the model is saturated only with oil. The initial pressure
% is set to the reference pressure.
s0 = [0, 1, 0];
state  = initResSol(G, pRef, s0);

%% SET UP SCHEDULE
% Time step is set to 30 days, with an initial ramp up to the designated
% time step.
dt = rampupTimesteps(totTime, 30*day, 10);
schedule = simpleSchedule(dt, 'W', W);

%% LAUNCH SIMULATION
[ws, states, report] = simulateScheduleAD(state, model, schedule);

%% PLOT RESULTS
figure;
plotToolbar(G,states);
view(30,45);
axis equal tight
grid on