%% Outline of Formations in the CO2 Storage Atlas: North Sea
% In this example, show the outline of all North Sea formations that can be
% constructed along with a map of Norway and point plots of all production
% wells in the Norwegian continental shelf.
%
% The well data comes from the Norwegian Petroleum Directorate and can be
% found in more detail at http://factpages.npd.no/factpages/.
%
% The map of Norway comes from The Norwegian Mapping and Cadastre Authority
% and can be found at http://www.kartverket.no/. Note that the map is only
% provided for scale and rough positioning - no claims are made regarding
% the accuracy in relation the subsea reservoirs.
grdecls = getAtlasGrid(getNorthSeaNames(),'coarsening', 3);
ng = numel(grdecls);

grids = cell(ng,1);
for i = 1:ng
    gd = processGRDECL(grdecls{i});
    grids{i} = computeGeometry(gd(1));
end

%%
clf
mymap = [
         0         0    1.0000
         0    1.0000         0
    1.0000         0         0
         0         0    0.5000
         0    0.5000         0
    0.7500         0         0         
         0    1.0000    1.0000
    1.0000         0    1.0000
    1.0000    1.0000         0
    0.8500    0.8500    0.8500
    0.4000    0.4000    0.4000    
         0         0         0
    0.2500    0.5000    1.0000
    0.5000    0.2500    0.7500
    0.7500    0.2500    0.5000
    0.7500    0.5000    0.2500
    0.2500    1.0000    0.5000
    0.5000    1.0000    0.2500
         ];

for i=1:ng;
    G = grids{i};
    bf = boundaryFaces(G);
    ind = G.faces.normals(bf,3)==0;
    plotFaces(G,bf(ind), 'FaceColor', 'none', 'EdgeColor', mymap(i,:), 'LineWidth',2);
end
legend(cellfun(@(x) x.name, grdecls, 'UniformOutput', false), 'Location', 'EastOutside')
box on
view(2)
set(gcf,'Color','none');set(gca,'Color',[1 1 1]);
set(gca,'XColor',[0,0,0])
set(gca,'YColor',[0,0,0]);

ax = axis();
load(fullfile(getDatasetPath('co2atlas'), 'norway.mat'));
for k=1:length(norway), 
    line(norway{k}(:,1) + 6.8e5, norway{k}(:,2)); 
end;
axis(ax)

hold on
load(fullfile(getDatasetPath('co2atlas'), 'welldata.mat'));
plot(welldata(:,2), welldata(:,1), '.k', 'MarkerSize', 5)
hold off
