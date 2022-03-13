%% Trapping on Johansen grids
% The Johansen formation has been proposed as a potential injection site
% for CO2, in particular when it was planned to capture CO2 from the
% gas-power plant at Mongstad. A simplified representation of the Johansen
% formation was also used as a benchmark case in the Stuttgart code
% comparison study based on a paper by Eigestad et al. Herein, we consider
% the same injection point as suggested by Eigestad et al. and evaluate the
% potential for structural trapping that can be estimated from two
% different models: (i) a sector model developed by the Norwegian Petroleum
% Directorate, which was used to produce the Stuttgart benchmark test, and
% (ii) a model of a somewhat larger region derived from the CO2 Storage
% Atlas for the Norwegian North Sea.

mrstModule add co2lab libgeometry deckformat coarsegrid matlab_bgl

%% Load NPD data: sector model
try
   jdir = fullfile(mrstPath('co2lab'), 'data', 'johansen');
   sector = 'NPD5';
   sector = fullfile(jdir, sector);
   grdecl = readGRDECL([sector '.grdecl']);
catch me
   disp(' -> Download data from: http://www.sintef.no/Projectweb/MatMoRA/')
   disp(['    Putting data in ', jdir]);
   unzip('https://www.sintef.no/project/MatMoRA/Johansen/NPD5.zip', jdir);
   grdecl = readGRDECL([sector '.grdecl']);
end

% Extract the part that represents the Johansen formation
grdecl = cutGrdecl(grdecl,[1 grdecl.cartDims(1); 1 grdecl.cartDims(2);  6 11]);
try
   Gs = processgrid(grdecl);
   Gs = mcomputeGeometry(Gs);
catch
   Gs = processGRDECL(grdecl);
   Gs = computeGeometry(Gs);
end

Gts = topSurfaceGrid(Gs);

% Get the position of the well (data given from 'Sector5_Well.txt');
wi    = find(Gs.cells.indexMap==sub2ind(Gs.cartDims, 48, 48, 1));
wcent = Gs.cells.centroids(wi,:);
d = sqrt(sum(bsxfun(@minus, Gts.cells.centroids, wcent(1:2)).^2, 2));
[~,wi_s] = min(d);

%% Load NPD data: full-field model
try
   jdir = fullfile(mrstPath('co2lab'), 'data', 'johansen');
   sector = 'FULLFIELD_IMAXJMAX';
   sector = fullfile(jdir, sector);
   grdecl = readGRDECL([sector '.GRDECL']);
catch me
   disp(' -> Download data from: http://www.sintef.no/Projectweb/MatMoRA/')
   disp(['    Putting data in ', jdir]);
   unzip('https://www.sintef.no/project/MatMoRA/Johansen/FULLFIELD_Eclipse.zip', jdir);
   grdecl = readGRDECL([sector '.GRDECL']);
end

% Extract the part that represents the Johansen formation
grdecl = cutGrdecl(grdecl,[1 grdecl.cartDims(1); 1 grdecl.cartDims(2);  10 14]);
try
   Gf = processgrid(grdecl);
   Gf = mcomputeGeometry(Gf);
catch
   Gf = processGRDECL(grdecl);
   Gf = computeGeometry(Gf);
end
Gtf = topSurfaceGrid(Gf);

% Get the position of the well
d = sqrt(sum(bsxfun(@minus, Gtf.cells.centroids, wcent(1:2)).^2, 2));
[~,wi_f] = min(d);

%% Load atlas data
grdecl = getAtlasGrid('Johansenfm');
try
   Ga = processgrid(grdecl{1});
   Ga = mcomputeGeometry(Ga);
catch
   Ga = processGRDECL(grdecl{1});
   Ga = computeGeometry(Ga);
end
Gta = topSurfaceGrid(Ga);

% Get the position of the well
d = sqrt(sum(bsxfun(@minus, Gta.cells.centroids, wcent(1:2)).^2, 2));
[~,wi_a] = min(d);

%% Plot the three data sets
clf
zm = min(Ga.nodes.coords(:,3));
zM = max(Ga.nodes.coords(:,3));

plotCellData(Ga,Ga.cells.centroids(:,3),'FaceAlpha',.95,'EdgeColor','none');
plotGrid(Gf,'FaceColor','none','EdgeAlpha',.2,'EdgeColor','r');
plotGrid(Gs,'FaceColor','none','EdgeAlpha',.4,'EdgeColor','k');
axis tight, view(-62,60);
light, lighting phong,
light('Position',[max(Gta.cells.centroids) 4*zM],'Style','infinite');

hold on; plot3(wcent([1 1],1),wcent([1 1],2),[zm zM],'b','LineWidth',2);

legend('Atlas','Full field','Sector','Inj.pt',...
   'Location','SouthOutside','Orientation','horizontal');

%% Interactive trapping: atlas grid
interactiveTrapping(Gta, 'method', 'cell', 'light', true, ...
   'spillregions', true, 'colorpath', false, 'injpt', wi_a);
view(-80,64);

%% Interactive trapping: NPD sector grid
interactiveTrapping(Gts, 'method', 'cell', 'light', true, ...
   'spillregions', true, 'colorpath', false, 'injpt', wi_s);
view(-80,64);


%% Interactive trapping: NPD full-field grid
interactiveTrapping(Gtf, 'method', 'cell', 'light', true, ...
   'spillregions', true, 'colorpath', false, 'injpt', wi_f);
view(-80,64);

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
