%test_markcells

% close all;
clearvars -except METISPATH mrstVerbose screenSize GlobTri
dbclear if error
% 
pth = getDatasetPath('bedmodel2');

grdecl = readGRDECL(fullfile(pth, 'BedModel2.grdecl'));
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

G = processGRDECL(grdecl);
G = computeGeometry(G);


rock = grdecl2Rock(grdecl, G.cells.indexMap);


% nx = 20; ny = 20; nz = 20;
% lx = 30; ly = 30; lz = 5;
% G = cartGrid([nx ny nz], [lx ly lz]);
% G = computeGeometry(G);

% global GlobTri
% if isempty(GlobTri)
%     GlobTri = globalTriangulation(G,GlobTri);
% elseif size(GlobTri.Tri.Points,1) ~= G.nodes.num
%     GlobTri = globalTriangulation(G,GlobTri);
% end

load('GlobTri_BedModel');
%%
fracplanes = planes_input();
checkIfCoplanar(fracplanes)
%
fracplanes = getPlaneNormals(fracplanes);
Sign = establishSign(G, fracplanes);

%
[fraCells,remove] = markcells(G, fracplanes, 'Sign', Sign, 'GlobTri', GlobTri);
if ~isempty(remove)
    fracplanes = fracplanes(setdiff(1:numel(fracplanes),remove));
end

for i = 1:numel(fracplanes)
    figure; plotGrid(G,'FaceColor','none','EdgeAlpha',0.01); hold on; 
    plotPlane(fracplanes)
    plotGrid(G,fraCells{i,1},'FaceColor','r','FaceAlpha',0.1);
end
