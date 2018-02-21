%% Bed Models 1
pth = getDatasetPath('BedModels1');
fn = dir(fullfile(pth,'*.grdecl'));
for i = 1:numel(fn)
    grdecl = readGRDECL(fullfile(pth, fn(i).name));
    G = processGRDECL(grdecl);
    rock = grdecl2Rock(grdecl,G.cells.indexMap);
    clf
    if isfield(grdecl,'SATNUM')
        q = grdecl.SATNUM(G.cells.indexMap);
    else
        q = rock.poro;
    end
    plotCellData(G, q, 'EdgeAlpha',.1, 'EdgeColor','k');
    view(3); axis tight off, drawnow
    print('-dpng',sprintf('bedmodels1-%d.png', i));
end

%% Bed Model 2
grdecl = readGRDECL(fullfile(getDatasetPath('BedModel2'),'BedModel2.grdecl'));
G      = processGRDECL(grdecl);
rock   = grdecl2Rock(grdecl,G.cells.indexMap);
clf, 
plotCellData(G,grdecl.SATNUM(G.cells.indexMap),'EdgeColor','none');
view(3); axis tight off

