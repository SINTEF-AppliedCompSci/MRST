%% Bjarmeland Platform (BP) in Barents Sea.
% Simulation studies have been undertaken in this formation by NPD.
% The following specifications were used:
% - injection period of 50 years
% - simulation continues for 1000 years

% - prospect A details found on page 139 (chp 6) of Atlas.
% - well 7125/1-1
% - closed structure (but table says open).
% - residual oil present
% - hydrostatic pressure
% - max pressure increase a) 75 bars, b) 30 bars
% - reported storage capacity: 176 Mt

% - prospect B details found on page 140 (chp 6) of Atlas.
% - well 7124/4-1S
% - half-open structure, west bdry closed
% - reported storage capacity: 19 Mt



%% Hammerfest Basin Aquifer (HB)
% - half-open aquifer, with prospects C-H.
%
% - prospect C details found on page 141 (chp 6) of Atlas.
% - well 7122/4-1





%% PLOT all ATLAS FORMATIONS:

% map of Norway might be in Euref89(wgs84) sone 33
% Euref89 - 
% wgs84 - world geodetic system 1984, has no zones

% taken from CAGEO-75:
%% Outline of Formations in the CO2 Storage Atlas
% In this example, show the outline of all formations that can be
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

% look at maps in slides by NDP, for better location referencing.

% Note: map of Norway is likely in Euref89(wgs84) sone 33
% Note: formation maps are in ED50 UTM32 (i.e., zone 32)

%grdecls = getAtlasGrid('coarsening', 3);
grdecls = getAtlasGrid({'Utsirafm', 'Huginfm', 'Johansenfm', 'Pliocenesand', ...
     'Skadefm', 'Statfjordfm', ...
     'Brynefm', ...
     'Bjarmeland', 'Stofm', 'Nordmela', 'Tubaenfm', 'Knurrfm', ...
     'Tiljefm', 'Arefm', 'Rorfm', 'Notfm', 'Ilefm', 'Garnfm'}, 'coarsening', 3);
%grdecls = getAtlasGrid('coarsening',3);


ng = numel(grdecls);

grids = cell(ng,1);
for i = 1:ng
    gd = mprocessGRDECL(grdecls{i});
    grids{i} = mcomputeGeometry(gd(1)); % use mcomputeGeometry for large grid sizes (i.e., Garn)
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
     0         0         0
      0         0         0
       0         0         0
        0         0         0
         0         0         0
          0         0         0
           0         0         0
            0         0         0
             0         0         0
              0         0         0
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
%set(gcf,'Color','none'); set(gca,'Color',[1 1 1]);
%set(gca,'XColor',[0,0,0])
%set(gca,'YColor',[0,0,0]);
ax = axis();

load(fullfile(getDatasetPath('co2atlas'), 'norway.mat'));
a = 0.65e6 %0.65e6; %0.3954e6;
b = 0.5 %0.4; %0.4; %0.78364;
for k=1:length(norway),
    %line(norway{k}(:,1), norway{k}(:,2)); 
    line(a+b*(norway{k}(:,1)), norway{k}(:,2)); 
    %line(norway{k}(:,1) + 6.8e5, norway{k}(:,2)); 
end;
axis equal
xlim([0 14e5])
ylim([6200000 8400000])

axis(ax)

hold on
load(fullfile(getDatasetPath('co2atlas'), 'welldata.mat'));
plot(welldata(:,2), welldata(:,1), '.k', 'MarkerSize', 5)
hold off