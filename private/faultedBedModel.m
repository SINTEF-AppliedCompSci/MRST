pth = getDatasetPath('bedmodel2');
grdecl = readGRDECL(fullfile(pth, 'BedModel2.grdecl'));

usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

z = grdecl.ZCORN;
[ii, jj, kk, ll] = ind2sub(2*grdecl.cartDims, 1:numel(z));
fault1 = ii > 30;
fault2 = jj > 30;

z(fault1) = z(fault1) + 1.25;
z(fault2) = z(fault2) + 1.75;
grdecl.ZCORN = z(:);
%%
G = processGRDECL(grdecl);
G = computeGeometry(G);
rock = grdecl2Rock(grdecl, G.cells.indexMap);
%%
mrstModule add mrst-gui
figure;
plotToolbar(G, rock);
axis tight off
view(30, 30);

figure;
plotGrid(G, 'facecolor', 'none');
plotFaces(G, G.faces.tag > 0);
view(30, 30);
