%% Show fluid models from SPE 1 and SPE 9
% The purpose of this example is to demonstrate how to read an Eclipse
% input file, create a fluid model, and inspect it. To this end, we will
% consider the SPE1 and SPE9 data sets. The example assumes that you have
% downloaded them already. If not, you should uncomment and run the
% following line, to bring up the dataset manager that enables you to
% download the data sets

% mrstDatasetGUI

%% Read data
% If you have used the mrstDatasetGUI functionality to download the SPE 1
% and SPE 9 data set, it should be possible to find them on a standard
% path. We read the data and convert to SI units.
mrstModule add deckformat

pth  = getDatasetPath('spe1');
fn   = fullfile(pth, 'BENCH_SPE1.DATA');
%pth  = getDatasetPath('spe9');
%fn   = fullfile(pth, 'BENCH_SPE9.DATA');
deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);

%% Create fluid objects
mrstModule add ad-props ad-core
f = initDeckADIFluid(deck);

%% Plot two-phase relative permeability curves
% For a three-phase model we have four relative permeability curves. One
% for both gas and water and two curves for the oil phase. The oil relative
% permeability is tabulated for both water-oil and oil-gas systems.
s = linspace(0,1,51)';

figure;
plot(s, f.krW(s), 'linewidth', 2)
grid on
xlabel('Water saturation');
title('Water relative permeability curve')
ylabel('k_r')

figure;
plot(s, [f.krOW(s), f.krOG(s)], 'linewidth', 2)
grid on
xlabel('Oil saturation');
legend('Oil-Water system', 'Oil-Gas system', 'location', 'northwest')
title('Oil relative permeability curves')
ylabel('k_r')

figure;
plot(s, f.krG(s), 'linewidth', 2)
grid on
xlabel('Gas saturation');
title('Gas relative permeability curve')
ylabel('k_r')

%% Plot three-phase relative permeability for oil
% When all three phases are present simultaneously in a single cell, we
% need to use some functional relationship to combine the two-phase curves
% in a reasonable manner, resulting in a two-dimensional relative
% permeability model. Herein, we use a simple linear interpolation, which
% is also the default in Eclipse
%

[x, y] = meshgrid(linspace(0,1,201));
[~, krO] = f.relPerm(x(:),y(:));
krO = reshape(krO,size(x));
krO(x+y>1)=nan;
z=1-x-y;
figure;
[mapx,mapy] = ternaryAxis('names',{'S_w';'S_g';'S_o'});
contourf(mapx(x,y,z), mapy(x,y,z), krO, 20)
set(gca,'YLim',[0 sqrt(3)/2]);
text(mapx(.45,.45,.1),mapy(.45,.45,.1),'Oil relative permeability', ...
    'BackgroundColor','w','HorizontalAlignment','center'); 