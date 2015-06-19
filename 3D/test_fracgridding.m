%test_gridfracture

close all;
clearvars -except METISPATH mrstVerbose screenSize verbose GlobTri
% dbstop if error
dbclear if error
% nx = 20; ny = 20; nz = 20;
% lx = 30; ly = 30; lz = 5;
% G = cartGrid([nx ny nz], [lx ly lz]);
% G = computeGeometry(G);

pth = getDatasetPath('bedmodel2');
grdecl = readGRDECL(fullfile(pth, 'BedModel2.grdecl'));
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

G = processGRDECL(grdecl);
G = computeGeometry(G);


% Gscaled = G;
% Gscaled.nodes.coords = [G.nodes.coords(:,1)/dims(1) , ...
%                         G.nodes.coords(:,2)/dims(2) , ...
%                         G.nodes.coords(:,3)/dims(3) ];
% Gscaled = computeGeometry(Gscaled);
% minfA = min(Gscaled.faces.areas);
%%
fracplanes = planes_input();
checkIfCoplanar(fracplanes)
%%
dims = max(G.nodes.coords);
fracScaled = fracplanes; 
Ar = zeros(numel(fracScaled),1);
Ar_scaled = zeros(numel(fracScaled),1);
for i = 1:numel(fracScaled)
    fracScaled(i).points = [fracScaled(i).points(:,1)./dims(1), ...
                            fracScaled(i).points(:,2)./dims(2), ...
                            fracScaled(i).points(:,3)./dims(3)];
    Ar(i) = polyArea3D(fracplanes(i).points);
    Ar_scaled(i) = polyArea3D(fracScaled(i).points);
end
scale = 500*max(G.cells.volumes./dims(1)./dims(2)./dims(3));
esize_scaled = scale/max(Ar_scaled);
%
fracplanes = getPlaneNormals(fracplanes);
fracScaled = getPlaneNormals(fracScaled);
rectangular = findRectangularFractures(fracplanes);

%
for i = 1:numel(fracplanes)
isrectangular = any(rectangular==i);
Gf = gridFractureDistmesh(G, fracplanes(i), fracScaled(i), esize_scaled, ...
     'type', 'pebi','rectangular',isrectangular, 'scale', max(dims(3)/dims(1),dims(1)/dims(3)));
plotGrid(G,'FaceAlpha',0,'EdgeAlpha',0.05); hold on; plotGrid(Gf); figure;
end