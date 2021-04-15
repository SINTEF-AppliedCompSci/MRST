%% Model 2 of the 10th SPE Comparative Solution Project
% The data set was originally posed as a benchmark for upscaling methods.
% The 3-D geological model consists of 60x220x85 grid cells, each of size
% 20ftx10ftx2ft. The model is a geostatistical realization from the
% Jurassic Upper Brent formations, in which one can find the giant North
% Sea fields of Statfjord, Gullfaks, Oseberg, and Snorre. In this specific
% model, the top 70 ft (35 layers) represent the shallow-marine Tarbert
% formation and the lower 100 ft (50 layers) the fluvial Ness formation.
% The data can be obtained from the SPE website:
% http://www.spe.org/web/csp/
% The data set is used extensively in the literature, and for this reason,
% MRST has a special module that will download and provide easy access to
% the data set.
%
% In this example, we will inspect the SPE10 model in more detail.

mrstModule add spe10

%% Load the model
% The first time you access the model, the MRST dataset system will prompt
% you to download the model from the official website of the comparative
% solution project. Depending upon your internet connection this may take
% quite a while--even several minutes--so please be patient. Notice that
% |getSPE10rock| returns permeabilities in strict SI units (i.e., m^2), and
% the petrophysical data of this model may therefore be used directly in
% simulations in MRST.
rock = getSPE10rock();
p = reshape(rock.poro,60,220,85);
Kx = reshape(log10(rock.perm(:,1)),60,220,85) - log10(milli*darcy);
Ky = reshape(log10(rock.perm(:,2)),60,220,85) - log10(milli*darcy);
Kz = reshape(log10(rock.perm(:,3)),60,220,85) - log10(milli*darcy);

%% show p
clf
slice( p, [1 220], 60, [1 85]);
shading flat, axis equal, set(gca,'zdir','reverse'), box on;
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
colorbarHist(p(:),[-.01 max(p(:))],'South',100);

%% show Kx/Ky
% These are identical so we only plot one of them
clf
slice( Kx, [1 220], 60, [1 85]);
shading flat, axis equal, set(gca,'zdir','reverse'), box on;
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
h=colorbarHist(Kx(:), caxis, 'South', 100);
set(h,'XTickLabel',10.^get(h,'XTick'));
set(h,'YTick',mean(get(h,'YLim')),'YTickLabel','mD');

%% show Kz
clf
slice( Kz, [1 220], 60, [1 85]);
shading flat, axis equal, set(gca,'zdir','reverse'), box on;
set(gca,'XTick',[],'YTick',[],'ZTick',[]);
h=colorbarHist(Kz(:), caxis, 'South', 100);
set(h,'XTickLabel',10.^get(h,'XTick'));
set(h,'YTick',mean(get(h,'YLim')),'YTickLabel','mD');

%% show horizontal permeability distributions
clf
hist(Kx(220*60*35+1:end),101);
hold on
hist(Kx(1:220*60*35),101)
hold off
h=get(gca,'Children');
set(h(1),'EdgeColor',[0 0 0.4],'FaceColor','none')
set(h(2),'EdgeColor',[0.7 0 0],'FaceColor','none')
h=legend('Ness','Tarbert'); set(h,'FontSize',16);

%% show vertical permeability distributions
clf;
hist(Kz(1:220*60*35),101)
hold on
hist(Kz(220*60*35+1:end),101);
hold off
h=get(gca,'Children');
set(h(2),'EdgeColor',[0 0 0.4],'FaceColor','none')
set(h(1),'EdgeColor',[0.7 0 0],'FaceColor','none')
h=legend('Tarbert','Ness');set(h,'FontSize',16);

%% show porosity distributions
pT = p(1:220*60*35); pN = p(220*60*35+1:end);
pb = linspace(0,0.5,101);
N=histc(pN,pb); bar(pb,N,'histc');
hold on
N=histc(pT,pb); bar(pb,N,'histc');
hold off;
h=get(gca,'Children');
set(h(1),'EdgeColor',[0 0 0.4],'FaceColor','none')
set(h(2),'EdgeColor',[0.7 0 0],'FaceColor','none')
h=legend('Ness','Tarbert');set(h,'FontSize',16);
set(gca,'XLim',[0 0.5]);

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
