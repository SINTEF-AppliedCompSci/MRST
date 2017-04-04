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
mrstModule add deckformat ad-props ad-core

nr   = [1 3 9];
f    = cell(3,1);
col  = lines(3);
name = cell(3,1);
for i=1:3
    name{i} = sprintf('spe%d',nr(i));
    pth  = getDatasetPath(name{i});
    fn   = fullfile(pth, sprintf('BENCH_SPE%d.DATA',nr(i)));
    deck = readEclipseDeck(fn);
    deck = convertDeckUnits(deck);
    f{i} = initDeckADIFluid(deck);
end

%% Plot two-phase relative permeability curves
% For a three-phase model we have four relative permeability curves. One
% for both gas and water and two curves for the oil phase. The oil relative
% permeability is tabulated for both water-oil and oil-gas systems.
s = linspace(0,1,51)';
for i=1:numel(name), name{i} = upper(name{i}); end

figure; hold all, set(gca,'FontSize',14)
for i=1:3
    plot(s, f{i}.krW(s), 'linewidth', 2, 'Color',col(i,:))
end
grid on
legend(name{:},2);
xlabel('Water saturation');
title('Water relative permeability curve')

figure; hold all, set(gca,'FontSize',14)
for i=1:3
    krOW = f{i}.krOW(s); krOW(s>1-f{i}.sWcon)=nan;
    krOG = f{i}.krOG(s); krOG(s>1-f{i}.sWcon)=nan;
    plot(s, krOW, '-','Color',col(i,:), 'linewidth',2);
    plot(s, krOG, '--', 'linewidth', 4, 'Color',col(i,:))
end
grid on
h = findobj(gca,'Type','line');
legend(h([6 5 6 4 2]),'oil/water system', 'oil/gas system', name{:},2);
xlabel('Oil saturation');
title('Oil relative permeability curves')

figure; hold all, set(gca,'FontSize',14)
for i=1:3
    krG = f{i}.krG(s); krG(s>1-f{i}.sWcon)=nan;
    plot(s, krG, 'linewidth', 2,'Color',col(i,:))
end
grid on
legend(name{:},4);
xlabel('Gas saturation');
title('Gas relative permeability curve')


%% Plot three-phase relative permeability for oil
% When all three phases are present simultaneously in a single cell, we
% need to use some functional relationship to combine the two-phase curves
% in a reasonable manner, resulting in a two-dimensional relative
% permeability model. Herein, we use a simple linear interpolation, which
% is also the default in Eclipse
%
[sw, sg] = meshgrid(linspace(0,1,201));
so=1-sw-sg;
for i=1:3
    [~, krO] = f{i}.relPerm(sw(:),sg(:));
    krO = reshape(krO,size(sw)); krO(sw+sg>1)=nan;
    
    figure; set(gca,'FontSize',14);
    [mapx,mapy] = ternaryAxis('names',{'S_w';'S_g';'S_o'},'FontSize',14);
    contourf(mapx(sw,sg,so), mapy(sw,sg,so), krO, 20)
    
    set(gca,'YLim',[0 sqrt(3)/2]);
    text(mapx(.45,.45,.1),mapy(.45,.45,.1),'Oil relative permeability', ...
        'BackgroundColor','w','HorizontalAlignment','center','FontSize',14);
    caxis([0 1]);
end