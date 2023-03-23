%% Modle 2 of the 10th SPE Comparative Solution Project
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
doprint = false;

%% Load the model
% The first time you access the model, it will be downloaded from the
% official websiet of the comparative solution project. Depending upon your
% internet connection, this may take quite a while, so please be patient.
rock = getSPE10rock();
p = rock.poro; 
K = convertTo(rock.perm,milli*darcy);

%% show p
slice( reshape(p,60,220,85), [1 220], 60, [1 85]); 
shading flat, axis equal off, set(gca,'zdir','reverse'), box on;
colorbar('horiz');
set(gcf,'renderer','painters');
if doprint
   print -depsc2 rock-spe10-poro.eps;                           %#ok<UNRCH>
   print -dpng rock-spe10-poro.png;
end

%% show Kx
slice( reshape(log10(K(:,1)),60,220,85), [1 220], 1, [1 85]); 
shading flat, axis equal off, set(gca,'zdir','reverse'), box on;
h=colorbar('horiz');
set(h,'XTickLabel',10.^get(h,'XTick'));
set(h,'YTick',mean(get(h,'YLim')),'YTickLabel','mD');
set(gcf,'renderer','painters');
if doprint
   print -depsc2 rock-spe10-Kx.eps;                             %#ok<UNRCH>
   print -dpng rock-spe10-Kx.png;
end

%% show Ky
slice( reshape(log10(K(:,2)),60,220,85), [1 220], 60, [1 85]); 
shading flat, axis equal off, set(gca,'zdir','reverse'), box on;
h=colorbar('horiz');
set(h,'XTickLabel',10.^get(h,'XTick'));
set(h,'YTick',mean(get(h,'YLim')),'YTickLabel','mD');
set(gcf,'renderer','painters');
if doprint
   print -depsc2 rock-spe10-Ky.eps;                             %#ok<UNRCH>
   print -dpng rock-spe10-Ky.png;
end

%% show Kz
slice( reshape(log10(K(:,3)),60,220,85), [1 220], 60, [1 85]); 
shading flat, axis equal off, set(gca,'zdir','reverse'), box on;
h=colorbar('horiz');
set(h,'XTickLabel',10.^get(h,'XTick'));
set(h,'YTick',mean(get(h,'YLim')),'YTickLabel','mD');
set(gcf,'renderer','painters');
if doprint
   print -depsc2 rock-spe10-Kz.eps;                             %#ok<UNRCH>
   print -dpng rock-spe10-Kz.png;
end

%% show horizontal permeability distributions
hist(log10(K(220*60*35+1:end,1)),101);
hold on
hist(log10(K(1:220*60*35,1)),101)
hold off
h=get(gca,'Children'); 
set(h(1),'EdgeColor',[0 0 0.4],'FaceColor','none')
set(h(2),'EdgeColor',[0.7 0 0],'FaceColor','none')
h=legend('Ness','Tarbert'); set(h,'FontSize',16);
if doprint
   print -depsc2 rock-spe10-Kx-hist.eps;                        %#ok<UNRCH>
end

%% show vertical permeability distributions
clf;
hist(log10(K(1:220*60*35,3)),101)
hold on
hist(log10(K(220*60*35+1:end,3)),101);
hold off
h=get(gca,'Children'); 
set(h(2),'EdgeColor',[0 0 0.4],'FaceColor','none')
set(h(1),'EdgeColor',[0.7 0 0],'FaceColor','none')
h=legend('Tarbert','Ness');set(h,'FontSize',16);
if doprint
   print -depsc2 rock-spe10-Kz-hist.eps;                        %#ok<UNRCH>
end

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
if doprint
   print -depsc2 rock-spe10-poro-hist.eps;                      %#ok<UNRCH>
end

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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
