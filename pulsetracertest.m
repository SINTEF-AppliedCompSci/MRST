% Setup pulse tracer case since NTPFA performed surprisingly bad on
% the notes/slides2/scripts/plotPulseAtPos using rate wells, homogeneous Neumann bc
% datacases   = {'small-perm', 'small-perm2', 'small-perm3', 'small-perm4'};
% datacase = datacases{1}

%{
summary
cart grid: tpfa=ntpfa, mfd a bit different
cart grid + twist: same as w/o twister
tet grid: tpfa zero, mfd ok. ntpfa so-so
tet grid + twist: tpfa zero, ntpfa pretty bad (conv problems on refine=4), mfd great but not equal to cart grid case (but it probably converges)
tet grid + twist low tol: same bad ntpfa
%}

clear all
close all

mrstModule add incomp mimetic  ad-blackoil ad-core postprocessing ...
    ad-props mpfa ntpfa-kobaisi ntpfa mrst-gui wellpaths blackoil-sequential

% Grid (twist later)
nx = 7;
ny = 5;
nz = 3;
N = [nx, ny, nz];
L = [10, 20, 30];
refine = 2;
G = cartGrid(N*refine, L);
%G = tetrahedralGrid(G.nodes.coords);
G = computeGeometry(G);
xyzmax = max(G.nodes.coords);
xyzmin = min(G.nodes.coords);
h = (xyzmax - xyzmin)./N;

% Rock has one layer with high perm
dim = 1;
zmean = mean(G.nodes.coords(:,dim));
dz = xyzmax(:,dim) - xyzmin(:,dim);
cells = (1:G.cells.num)';
zcells_ii = cells(abs(G.cells.centroids(:,dim) - zmean) < 0.2*dz);
perm_low = 1*milli*darcy;
perm_high = 100*milli*darcy;
rock = makeRock(G, perm_low, 1);
rock.perm(zcells_ii) = perm_high;
pv = sum(poreVolume(G,rock));

% Wells
W = [];
W = addWell(W, G, rock, well_cells(G, '1'), 'Type', 'rate', ...
            'comp_i', 1, 'radius', 0.01, ...
            'Val', 10, 'Name', 'I1');
W = addWell(W, G, rock, well_cells(G, '5'), 'Type', 'rate', ...
            'comp_i', 1, 'radius', 0.01, ...
            'Val', 1, 'Name', 'I2');
% bhp_val = 10.0;
% W = addWell(W, G, rock, well_cells(G, '1'), 'Type', 'bhp', ...
%             'comp_i', 1, ...
%             'Val', bhp_val, 'Name', 'I1');
% W = addWell(W, G, rock, well_cells(G, '5'), 'Type', 'bhp', ...
%             'comp_i', 1, ...
%             'Val', -bhp_val, 'Name', 'I2');

% Optional twist
%G = twister(G, 0.1);

figure
plotCellData(G,rock.perm)
hold on
plotWell(G,W)
view(3)

% Simulation setup
p0 = 1;
s0 = 1;
state0 = initState(G, W, p0, s0);
%bc_std = bc_hom_neumann(G);
bc_std = bc_dirichlet_side_pairs(G, 1);
bc_flow_ntpfa = convertBC2FlowNTPFA(G, bc_std);
%schedule = simpleSchedule(dt, 'W', W, 'bc', bc_std);
mu_value = 1;
rho_value = 1;
fluid = initSingleFluid('mu',mu_value,'rho',rho_value);
gravity off

refpt = [xyzmax(:,1:2), 0.5*(xyzmax(:,3)+xyzmin(:,3))] + h.*[-1.5,-1.5,-0.1];
%tol = 1e-3*h;
%refpt = [max(G.nodes.coords(:,1)), max(G.nodes.coords(:,2)), 0.5*(max(G.coodrdsmin(G.nodes.coords(:,3))] + tol.*[-1,-1,1];
%refpt = max(G.nodes.coords) - 2.5*h;
%refpt = mean(G.nodes.coords) - 1e-3*h;
refcell = findEnclosingCell(G, refpt);
assert(refcell > 0);

plotconc = false;
tracer = true;
results = {};

%% TPFA mrst version
results{end+1}.name = 'tpfa mrst';
disp(results{end}.name);
T = computeTrans(G,rock);
state = incompTPFA(state0, G, T, fluid, 'wells', W, 'bc', bc_std);
%figure,plotToolbar(G,state);plotWell(G,W);view(3)
%title(['Pressure ', results{end}.name]);h=colorbar;h.Label.String='Pressure [Pa]';clear h;
%figure,plotFaceData(G,cells_xdir, state.flux);colorbar,view(3)
%title(['flux ', results{end}.name])
if plotconc
    Qinj=sum(state.wellSol(1).flux);tt=pv/Qinj;nstep=100;dt=tt/nstep;
    state = tracerTransport_implicit(G,rock,W,state,dt,nstep);
    figure, plotCellData(G,state.cres(:,end));plotWell(G,W);
    view(3);colorbar;title(['Concentration ', results{end}.name]);
drawnow
end
results{end}.state = state;


%% ntpfa mrst
W = [];
W = addWell(W, G, rock, well_cells(G, '1'), 'Type', 'rate', ...
            'comp_i', [1,0], 'radius', 0.01, ...
            'Val', 10, 'Name', 'I1');
W = addWell(W, G, rock, well_cells(G, '5'), 'Type', 'rate', ...
            'comp_i', [1,0], 'radius', 0.01, ...
            'Val', 1, 'Name', 'I2');
fluidAD = initSimpleADIFluid();
state0 = initResSol(G, 0, [1, 0]);
results{end+1}.name = 'ntpfa mrst';
disp(results{end}.name);
model = PressureOilWaterModelNTPFA(G, rock, fluidAD);
schedule = simpleSchedule(1, 'W', W, 'bc', bc_std);
[~, states] = simulateScheduleAD(state0, model, schedule);
results{end}.state = states{1};


%% nonlinear TPFA
results{end+1}.name = 'ntpfa kobaisi';
mrstverbosestatus=mrstVerbose;
mrstVerbose on
disp(results{end}.name);
interpFace=findHAP(G,rock,bc_flow_ntpfa);
disp(['fraction of faces with centroids outside convex hull ', num2str(interpFace.percentage)]);
interpFace=correctHAP(G,interpFace);
OSflux=findOSflux(G,rock,bc_flow_ntpfa,interpFace);
%state=FlowNTPFA(G,bc_flow_ntpfa,fluid,W,OSflux,p0*ones(G.cells.num, ...
%                                                  1),1e-7,1000);
state=FlowNTPFA(G,bc_flow_ntpfa,fluid,W,OSflux,results{1}.state.pressure,1e-7,1000);
%figure,plotCellData(G,state.pressure/1e6);view(3)
%figure,plotToolbar(G,state);plotWell(G,W);view(3)
%title(['Pressure ', results{end}.name]);h=colorbar;h.Label.String='Pressure [Pa]';clear h;
%figure,plotFaceData(G,cells_xdir, state.flux);colorbar,view(3)
%title(['flux ', results{end}.name])
if plotconc
    Qinj=sum(state.wellSol(1).flux);tt=pv/Qinj;nstep=100;dt=tt/nstep;
    state=tracerTransport_implicit(G,rock,W,state,dt,nstep);
    figure, plotCellData(G,state.cres(:,end));plotWell(G,W);
    view(3);colorbar;title(['Concentration ', results{end}.name]);
    drawnow
end
results{end}.state = state;
mrstVerbose(mrstverbosestatus);


%% MFD (Must be at the end)
results{end+1}.name = 'mfd';
disp(results{end}.name);
state=incompMimetic(state0,G,computeMimeticIP(G,rock),fluid,'wells',W,'bc',bc_std);
%figure, plotCellData(G,state.pressure/1e6);view(3);
%figure,plotToolbar(G,state);plotWell(G,W);view(3)
%title(['Pressure ', results{end}.name]);h=colorbar;h.Label.String='Pressure [Pa]';clear h;
%figure,plotFaceData(G,cells_xdir, state.flux);colorbar,view(3)
%title(['flux ', results{end}.name])
if plotconc
    Qinj=sum(state.wellSol(1).flux);tt=pv/Qinj;nstep=100;dt=tt/nstep;
    state=tracerTransport_implicit(G,rock,W,state,dt,nstep);
    figure, plotCellData(G,state.cres(:,end));plotWell(G,W);
    view(3);colorbar;title(['Concentration ', results{end}.name]);
    drawnow
end
results{end}.state = state;

if tracer
    refsols = cell(numel(results), 1);
    Ws = cell(1,2);
    Ws{1} = W;
    Ws{1}(1).tracer = 1;
    Ws{1}(2).tracer = 0;
    Ws{2} = W;
    Ws{2}(1).tracer = 0;
    Ws{2}(2).tracer = 0;
    
    time  = 1e4;
    dt = time/100;
    dtvec = rampupTimesteps(time, dt);
    schedule = simpleSchedule(dtvec, 'W', W);
    schedule.step.control(2 : end) = 2;
    control(1) = schedule.control;
    control(2) = control(1);
    control(1).W = Ws{1};
    control(2).W = Ws{2};    
    schedule.control = control;
    tracermodel = TracerModel(G, rock, fluid, 'tracerNames', {'tracer'});

    for i = 1:numel(results)
        disp(results{i}.name)
        [wellsols, states] = simulateScheduleAD(results{i}.state, tracermodel, schedule);
        figure
        plotToolbar(G,states);colorbar,view(3)
        hold on
        plotGrid(G, refcell, 'facecolor', 'k');
        title(results{i}.name)
       
        refsol = zeros(numel(dtvec),1);
        for j = 1:numel(dtvec)
            refsol(j) = states{j}.tracer(refcell);
        end
        refsols{i} = refsol;
    end

    figure, hold on
    names = cell(numel(results), 1);
    for i = 1:numel(results)
        plot(cumsum(dtvec),refsols{i},'-')
        names{i} = results{i}.name;
    end
    legend(names)
end
