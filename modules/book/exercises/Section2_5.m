%% Exercise 2.5.1
newplot
mrstModule add spe10
rock = getSPE10rock(); p=rock.poro; K=rock.perm(:,1);
plot(log10(K),p,'.','MarkerSize',4);

%% Exercise 2.5.2
newplot
filename = fullfile(getDatasetPath('CaseB4'),'pillar_36x48.grdecl');
grdecl = readGRDECL(filename);
grdecl = convertInputUnits(grdecl, getUnitSystem('METRIC'));
G = processGRDECL(grdecl);

% Homogeneous data
rock1 = makeRock(G, 200*milli*darcy, 0.2);
subplot(1,2,1), plotCellData(G, rock1.poro,'EdgeAlpha',.1);
view(130,30); axis tight off

% Heterogeneous data
p = gaussianField(G.cartDims, [0.2 0.4], [11 3 3], 2.5);
K = p.^3.*(1e-5)^2./(0.81*72*(1-p).^2);
rock2 = makeRock(G, K(:), p(:));
subplot(1,2,2), plotCellData(G, rock2.poro,'EdgeAlpha',.1), 
view(130,30); axis tight off

%% Exercise 2.5.3
newplot
G = computeGeometry(cartGrid([58 48 10]));
load rock1.mat;
plotToolbar(G, K);

% And in case you did not find the structure:
% newplot
% plotGrid(cartGrid([1 1 1], [58 48 10]),'FaceColor','none');
% view(70,60), axis tight
% [i,j,k]=gridLogicalIndices(G);
% h = plotCellData(G, log10(K), k==7, 'EdgeAlpha',.1);

%% Exercise 2.5.4
fn = {'mortarTestModel', 'non-periodicTilted', 'periodicTilted', ...
    'simpleLayeredUniform', 'testModel1', 'testModel2'};
pth = getDatasetPath('BedModels1');
for i=1:numel(fn)
    grdecl = readGRDECL(fullfile(pth,[fn{i},'.grdecl']));
    grdecl = convertInputUnits(grdecl, getUnitSystem('METRIC'));
    clc; disp(grdecl)
    
    G = processGRDECL(grdecl);
    rock = grdecl2Rock(grdecl, G);
    figure(1); clf
    subplot(1,2,2), 
    plotCellData(G,log10(rock.perm(:,1)),'EdgeAlpha',.05); 
    view(3), axis tight off
    subplot(1,2,1)
    plotCellData(cartGrid(G.cartDims), log10(grdecl.PERMX), 'EdgeAlpha',.05)
    view(3), axis tight off
    
    if isfield(grdecl, 'SATNUM')
        figure(2), clf, hold all
        sn = grdecl.SATNUM(G.cells.indexMap);
        sval = unique(sn);
        for n=1:numel(sval)
            hist(log10(rock.perm(sn==sval(n),1)),101);
        end
        col = colorcube(numel(sval));
        h = get(gca,'Children');
        for n=1:numel(h)
            set(h(n), 'FaceColor','none','EdgeColor',col(n,:),'LineWidth',3);
        end
        set(gca,'XLim',log10([.01*milli*darcy 10*darcy]));
    end
    
    drawnow; pause
end

%% Exercise 2.5.5
mrstModule add incomp
gravity reset on

pth = getDatasetPath('BedModels1');
grdecl = readGRDECL(fullfile(pth,'mortarTestModel.grdecl'));
grdecl = convertInputUnits(grdecl, getUnitSystem('METRIC'));
G      = computeGeometry(processGRDECL(grdecl));
rock   = grdecl2Rock(grdecl, G);


fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);
bc    = pside([], G, 'TOP', 100.*barsa());
T     = simpleComputeTrans(G, rock);
sol   = simpleIncompTPFA(initResSol(G, 0.0), G, T, fluid, 'bc', bc);

clf
plotFaces(G, 1:G.faces.num, convertTo(sol.facePressure, barsa()));
set(gca, 'ZDir', 'reverse'), title('Pressure [bar]')
view(3), colorbar
set(gca,'DataAspect',[1 1 .2]); axis tight

