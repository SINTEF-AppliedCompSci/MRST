%% Inspect the Pliocenesand Formation
% This formation has a very low percentage (0.02%) of trapping compared to
% the overall volume of the whole model
try
    require coarsegrid deckformat mex
catch %#ok<CTCH>
    mrstModule add coarsegrid deckformat mex
end

grdecl = getAtlasGrid('Pliocenesand');
G      = processGRDECL(grdecl{1});
G      = computeGeometry(G(1));
Gt     = topSurfaceGrid(G);
ta     = trapAnalysis(Gt, false);



%%
% Next, we will perform a VE simulation. This could, of course, have been
% launched from inside the interactive viewer, but to make the example as
% reproducible as possible, we launch it manually from the outside.

%%
% Rerun for a longer time
petrodata.avgperm = 1.2*darcy;
petrodata.avgporo = 0.25;

%%
% Adding depth to the plio example
% It is to shallow for a real storage sight
depth=1200;


% cut grid to avoid calculation on not relevant domain
wpos=Gt.parent.cells.centroids(5280,1:2);
wpos(:,1)=4.85e5;
G=Gt.parent;
rm_cells=abs(Gt.cells.centroids(:,2)-wpos(:,2))>2.5e4;
G=removeCells(G,rm_cells);
G.nodes.coords(:,3)=G.nodes.coords(:,3)+depth;
G=computeGeometry(G);
Gt=topSurfaceGrid(G);
%clear G;

% set the permeability and porosity
rock.poro = repmat(petrodata.avgporo, Gt.cells.num, 1);
rock.perm = repmat(petrodata.avgperm, Gt.cells.num, 1);
rock2D    = averageRock(rock, Gt);
%% Find pressure boundary
% Setting boundary conditions is unfortunately a manual process and may
% require some fiddling with indices, as shown in the code below. Here, we
% need to find all outer vertical faces
i = any(Gt.faces.neighbors==0, 2);  % find all outer faces
I = i(Gt.cells.faces(:,1));         % vector of all faces of all cells, true if outer
j = false(6,1);                     % mask, cells can at most have 6 faces,
%j(1:4)=true;
j(1)=true;
%   extract east, west, north, south
J = j(Gt.cells.faces(:,2));         % vector of faces per cell, true if E,W,N,S
bcIxVE = Gt.cells.faces(I & J, 1);

%% Set time and fluid parameters

% Fluid data are taken from paper SPE 134891
gravity on

Ti  =   50*year;
%dTi =  5*year;
dTi =  2*year;
istep = linspace(0.1*year, dTi, 10)';
istep = [istep; ones(floor((Ti-sum(istep))/dTi), 1)*dTi];
istep = [istep; Ti-sum(istep)];

%Tm  = 1600*year;
Tm  = 3000*year;
dTm = 25*year;
%dTm = 10*year;
mstep = linspace(0.5*year, dTm, 5)';
mstep = [mstep; ones(floor((Tm-sum(mstep))/dTm),1)*dTm];
mstep = [mstep; Tm-sum(mstep)];

%%
amount = 5; %Mt/year
% convert mass flux at surface to surface reference volume/"surface volume
rhoc=760;rhow=1100;
rate = amount*1e9*kilogram./(year*rhoc*kilogram*meter^3);
% find well position
dist=sqrt(sum(bsxfun(@minus,Gt.cells.centroids(:,1:2),wpos).^2,2));
[dd,cellnum]=min(dist);
[ix,iy]=ind2sub(Gt.cartDims,Gt.cells.indexMap(cellnum));
wellIx = double([ix iy]);
% make well
assert(G.cartDims(3)==1)
W = createSampleWell([],Gt.parent, rock, cellnum ,     ...
                     'Type', 'rate', 'Val', rate, ...
                     'Radius', 0.125, 'Name', 'I','Comp_i',[0 0 1]);
%{
W      = verticalWell([], Gt.parent, rock, wellIx(1), wellIx(2), ...
    1, 'Type', 'rate', 'Val', rate, ...
    'Radius', 0.1, 'comp_i', [1,0], 'name', 'I','InnerProduct','ip_tpf');
                     %}
%WVE = convertwellsVE(W, Gt.parent, Gt, rock2D);
control = struct('W',[],'step',struct('val',[],'control',[]));
control.W = {W, []};
control.step.val = [istep; mstep];
control.step.control = [ones(size(istep));ones(size(mstep))*2];

bc = addBC([], bcIxVE, 'pressure', ...
    Gt.faces.z(bcIxVE)*rhow*norm(gravity));
bc.sat = zeros(size(bc.face));

aquifer.G  = G;
aquifer.Gt = Gt;
aquifer.rock   = rock;
aquifer.rock2D = rock2D;
aquifer.W  = W;
aquifer.bc = bc;

%% start loop over cases
s = setupSimCompVe(aquifer.Gt, aquifer.rock2D);
residual=true;
for k=1:3
    switch k
        case 1
            dissolution=false;
            fluid = makeFluidModel(aquifer, 'residual', residual, ...
                'dissolution', dissolution, 'fluidType', 'sharp interface');
            system = initADISystemVE({'Oil','Gas'}, aquifer.Gt, aquifer.rock2D, ...
                fluid, 'simComponents', s, 'VE', true);
        case 2
            dissolution=true;
            fluid = makeFluidModel(aquifer, 'residual', residual, ...
                'dissolution', dissolution, 'fluidType', 'sharp interface');
            fluid=rmfield(fluid,'dis_rate');
            system = initADISystemVE({'Oil', 'Gas','DisGas'}, aquifer.Gt, aquifer.rock2D,...
                fluid,'simComponents',s,'VE',true);
        case 3
            % deside which fluid to use
            dissolution=true;
            fluid = makeFluidModel(aquifer, 'residual', residual, ...
                'dissolution', dissolution, 'fluidType', 'sharp interface');
            system = initADISystemVE({'Oil', 'Gas','DisGas'}, aquifer.Gt, aquifer.rock2D,...
                fluid,'simComponents',s,'VE',true);
        otherwise
            error('No such case')
    end
    assert(fluid.rhoOS==rhow);
    assert(fluid.rhoGS==rhoc);
    
    system.nonlinear.linesearch    = false;
    system.nonlinear.maxIterations = 10;
    system.nonlinear.tol           = 1e-6;
    
    % Well in 2D model
    z  = aquifer.Gt.cells.z;
    clear x0;
    
    
    x0.pressure = z(:)*norm(gravity)*fluid.rhoOS;
    x0.s(:,1)   = ones(G.cells.num,1);
    x0.s(:,2)   = zeros(G.cells.num,1);
    x0.rs       = ones(G.cells.num,1)*0.0;
    x0.smax     = x0.s;
    x0.smin     = x0.s;
    x0.sGmax    = x0.s(:,2);
    
    
    if(dissolution)
        systemOG = initADISystemVE({'Oil', 'Gas','DisGas'}, aquifer.Gt, aquifer.rock2D,...
            fluid,'simComponents',s,'VE',true);
    else
        systemOG = initADISystemVE({'Oil','Gas'}, aquifer.Gt, aquifer.rock2D, ...
            fluid, 'simComponents', s, 'VE', true);
    end
    systemOG.nonlinear.linesearch    = false;
    systemOG.nonlinear.maxIterations = 10;
    systemOG.nonlinear.tol           = 1e-6;
    
    t2=tic;
    [wellSols, states] = runMrstADI(x0, Gt, system, control, ...
        'force_step', false, 'dt_min', 0.5*year, 'report_all', false, 'bc', bc);
    t2=toc(t2);
    save(['data/secondPlioExample_',num2str(depth),'_',num2str(k),'.mat'],'t2','states','wellSols','control')
end




