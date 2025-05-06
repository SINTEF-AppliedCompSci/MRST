%% Show fluid models from SPE 1, SPE 3 and SPE 9
% The purpose of this example is to demonstrate how to read an Eclipse
% input file, create a fluid model, and inspect it. To this end, we will
% consider the SPE1, SPE 3, and SPE9 data sets. The example assumes that
% you have downloaded them already. If not, you should uncomment and run
% the following line, to bring up the dataset manager that enables you to
% download the data sets

% mrstDatasetGUI

%% Read data
% If you have used the mrstDatasetGUI functionality to download the SPE 1
% and SPE 9 data set, it should be possible to find them on a standard
% path. We read the data and convert to SI units.
mrstModule add deckformat ad-props ad-core

nr    = [1 3 9];
fo    = cell(3,1);
col   = lines(3);
name  = cell(3,1);
decks = cell(3,1);
for i=1:3
    name{i}  = sprintf('spe%d',nr(i));
    pth      = getDatasetPath(name{i});
    fn       = fullfile(pth, sprintf('BENCH_SPE%d.DATA',nr(i)));
    decks{i} = readEclipseDeck(fn);
    decks{i} = convertDeckUnits(decks{i});
    fo{i}    = initDeckADIFluid(decks{i});
    if ~isfield(fo{i},'sWcon')
        fo{i}.sWcon = 0.0;
    end
    fo{i}    = assignRelPerm(fo{i});
end

%% Plot two-phase relative permeability curves
% For a three-phase model we have four relative permeability curves. One
% for both gas and water and two curves for the oil phase. The oil relative
% permeability is tabulated for both water-oil and oil-gas systems.
s = linspace(0,1,51)';
for i=1:numel(name), name{i} = upper(name{i}); end

figure; hold all, set(gca,'FontSize',14)
for i=1:3
    plot(s, fo{i}.krW(s), 'linewidth', 2, 'Color',col(i,:))
end
grid on
legend(name{:},'Location','NorthWest');
xlabel('Water saturation');
title('Water relative permeability curve')

figure; hold all, set(gca,'FontSize',14)
for i=1:3
    krOW = fo{i}.krOW(s); krOW(s>1-fo{i}.sWcon)=nan;
    krOG = fo{i}.krOG(s); krOG(s>1-fo{i}.sWcon)=nan;
    plot(s, krOW, '-','Color',col(i,:), 'linewidth',2);
    plot(s, krOG, '--', 'linewidth', 4, 'Color',col(i,:))
end
grid on
h = findobj(gca,'Type','line');
legend(h([6 5 6 4 2]),'oil/water system', 'oil/gas system', name{:}, ...
       'Location', 'NorthWest');
xlabel('Oil saturation');
title('Oil relative permeability curves')

figure; hold all, set(gca,'FontSize',14)
for i=1:3
    krG = fo{i}.krG(s); krG(s>1-fo{i}.sWcon)=nan;
    plot(s, krG, 'linewidth', 2,'Color',col(i,:))
end
grid on
legend(name{:},'Location','SouthEast');
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
    [~, krO] = fo{i}.relPerm(sw(:),sg(:));
    krO = reshape(krO,size(sw)); krO(sw+sg>1)=nan;
    
    figure; set(gca,'FontSize',14);
    [mapx,mapy] = ternaryAxis('names',{'S_w';'S_g';'S_o'},'FontSize',14);
    contourf(mapx(sw,sg,so), mapy(sw,sg,so), krO, 20)
    
    set(gca,'YLim',[0 sqrt(3)/2]);
    text(mapx(.45,.45,.1),mapy(.45,.45,.1),'Oil relative permeability', ...
        'BackgroundColor','w','HorizontalAlignment','center','FontSize',14);
    caxis([0 1]);
end

% =========================================================================

%% Phase behavior of live oil for the SPE 1 or SPE 9 model
% In the following, we will briefly show the phase behavior for gas and oil
% in the model. That is, we will show both the input data points and the
% interpolated data obtained from the functions in the fluid object

% Set input deck and fluid object
[deck, f] = deal(decks{1}, fo{1});

% Extract data from the input deck
pvto = deck.PROPS.PVTO{1};
rsd  = pvto.key([1:end end]);
pbp  = pvto.data([pvto.pos(1:end-1)' end],1);
Bod  = pvto.data([pvto.pos(1:end-1)' end],2);
muOd = pvto.data([pvto.pos(1:end-1)' end],3);

%% Plot Rs
figure, set(gca,'FontSize',14)
pargs = {'LineWidth', 2, 'MarkerSize',7,'MarkerFaceColor',[.5 .5 .5]};
plot(pbp/barsa,rsd,'-o',pargs{:}); xlabel('Pressure [bar]');

%% Plot oil formation-volume factor and viscosity for oil
% For Bo, we plot both the tabulated data as well as interpolated data at
% various combinations of dissolved gas and reservoir pressures. For each
% data point, we must first determine whether the state is saturated or not
% by comparing the actual amount of dissolved gas against the maximum
% possible amount of solved gas for this pressure.
[M, N]    = deal(11,51);
[RsMax,pMax] = deal(max(rsd), max(pbp));
[rs,p]    = meshgrid(linspace(10,RsMax-10,M), linspace(0,pMax,N));
Rv        = reshape(f.rsSat(p(:)), N, M);
isSat     = rs >= Rv;
rs(isSat) = Rv(isSat);
Bo        = 1./reshape(f.bO(p(:), rs(:), isSat(:)),N,M);
muO       = reshape(f.muO(p(:), rs(:), isSat(:)),N,M);

% Formation-volume factor
figure, hold on, set(gca,'FontSize',14)
for j=1:M
    i = isSat(:,j);
    plot(p(i,j)/barsa, Bo(i,j),'b-',p(~i,j)/barsa,Bo(~i,j),'-r');
end
plot(pbp/barsa,Bod,'-bo',pargs{:});
hold off, axis tight, xlabel('Pressure [bar]');
title('Oil formation-volume factor [-]');

% Viscosity
figure, hold on, set(gca,'FontSize',14)
for j=1:M
    i = isSat(:,j);
    plot(p(i,j)/barsa,  convertTo(muO(i,j),  centi*poise),'b-',...
         p(~i,j)/barsa, convertTo(muO(~i,j), centi*poise),'-r');
end
plot(pbp/barsa, convertTo(muOd, centi*poise),'-bo',pargs{:});
hold off, axis tight, xlabel('Pressure [bar]');
title('Oil viscosity [cP]');

%% Phase behavior for the dry gas (SPE 1 or SPE9)
pvdg = deck.PROPS.PVDG{1};

figure, set(gca,'FontSize',14);
plot(pvdg(:,1)/barsa,pvdg(:,2),'-o',pargs{:})
set(gca,'YScale','log')
xlabel('Pressure [bar]'); 
title('Gas formation-volume factor [-]');

figure, set(gca,'FontSize',14);
plot(pvdg(:,1)/barsa,convertTo(pvdg(:,3), centi*poise), '-o',pargs{:})
xlabel('Pressure [bar]');
title('Gas viscosity [cP]');


%% Show inconsistencies for large pressures (only SPE1)
% For sufficintly high pressure/gas-oil ratios, the interpolated shrinkage
% factors will cross zero, which causes blowup and negative Bo values. In
% turn, this leads to negative partial volumes.
[M, N]    = deal(101,201);
[rs,p]    = meshgrid(linspace(f.rsSat(100*barsa),f.rsSat(1300*barsa),M), ...
              linspace(100,1300,N)*barsa);
Rv        = reshape(f.rsSat(p(:)), N, M);
isSat     = rs >= Rv;
rs(isSat) = Rv(isSat);
bo        = reshape(f.bO (p(:), rs(:), isSat(:)),N,M);
Bo        = 1./bo;
muO       = reshape(f.muO(p(:), rs(:), isSat(:)),N,M);
bg        = reshape(f.bG (p(:)), N, M);
Bg        = 1./bg;

figure
contourf(p/barsa,rs,bo,[-inf linspace(-.2,.9,12) inf]); axis tight
set(gca,'FontSize',20); xlabel('Pressure [bar]'); ylabel('Gas-oil ratio');
cb = colorbar; set(cb,'FontSize',20);
colormap(repmat([linspace(.2,.3,4) linspace(.6,1,20)]',1,3));

figure
contourf(p/barsa,rs,1./bo - rs./bg,[-inf linspace(-10,10,41) inf]); axis tight
set(gca,'FontSize',20); xlabel('Pressure [bar]');ylabel('Gas-oil ratio');
caxis([-10,10]); cb = colorbar; set(cb,'FontSize',20);
colormap(repmat([linspace(.25,.5,20) linspace(.7,.95,20)]',1,3));

% =========================================================================
%% Dead oil and gas with vaporized oil (SPE3)
% The SPE 3 benchmark was designed to study gas cycling in a richâ€?gas
% reservoir with retrograde condensation. The original setup was for
% compositional simulation. The current input file describes a black-oil
% translation of the compositional model.
% into a 
[deck, f] = deal(decks{2}, fo{2});
pvdo = deck.PROPS.PVDO{1};

figure, set(gca,'FontSize',14);
pargs = {'LineWidth', 2, 'MarkerSize',7,'MarkerFaceColor',[.5 .5 .5]};
plot(pvdo(:,1)/barsa, pvdo(:,2),'-o',pargs{:})
xlabel('Pressure [bar]'); 
title('Oil formation-volume factor [-]');

figure, set(gca,'FontSize',14);
plot(pvdo(:,1)/barsa, convertTo(pvdo(:,3), centi*poise), '-o',pargs{:})
xlabel('Pressure [bar]');
title('Oil viscosity [cP]');

%% Rich gas with retrograde condensation
% Extract data from the input deck
pvtg = deck.PROPS.PVTG{1};
pbp  = convertTo(pvtg.key,barsa);
rvd  = pvtg.data([1 pvtg.pos(2:end-1)'],1);
Bgd  = pvtg.data(:,2);
muGd = convertTo(pvtg.data(:,3), centi*poise);
p    = pvtg.key(rldecode((1:numel(pvtg.key))',diff(pvtg.pos),1))/barsa;

%% Plot Rv
figure, set(gca,'FontSize',14)
hold on
for i=1:numel(pbp)
    ind = pvtg.pos(i):pvtg.pos(i+1)-1;
    plot(p(ind),pvtg.data(ind,1),'-ok',...
        'MarkerSize',6,'MarkerFaceColor',[.5 .5 .5]);
end
plot(pbp,rvd,'-o',pargs{:}); xlabel('Pressure [bar]');
title('Maximum vaporized oil-gas ratio [-]');
hold off

%% Plot Bg
figure, set(gca,'FontSize',14)
plot(pbp,Bgd([pvtg.pos(2:end-1)'-1 end]),'-r', ...
    pbp, Bgd(pvtg.pos(1:end-1)), '-bo', pargs{:});
xlabel('Pressure [bar]');
title('Gas formation-volume factor [-]');

% Make inset
axes('position',[.48 .45 .4 .4]);
hold on
for i=1:numel(pbp)
    ind = pvtg.pos(i):pvtg.pos(i+1)-1;
    plot(p(ind),pvtg.data(ind,2),'-ok',...
        'MarkerSize',6,'MarkerFaceColor',[.5 .5 .5]);
end
plot(pbp,Bgd([pvtg.pos(2:end-1)'-1 end]),'-r', ...
    pbp, Bgd(pvtg.pos(1:end-1)), '-bo', pargs{:});
hold off
axis([205 240 4.2e-3 4.8e-3]);

%% Plot gas density
bG = f.bG(p*barsa, pvtg.data(:,1), false(size(p)));
rho = bG*f.rhoGS + bG.*pvtg.data(:,1)*f.rhoOS;
figure, hold on, set(gca,'FontSize', 14)
for i=1:numel(pbp),
    ind = pvtg.pos(i):pvtg.pos(i+1)-1;
    plot(p(ind),rho(ind),'-ok',...
        'MarkerSize',6,'MarkerFaceColor',[.5 .5 .5]);
end
plot(pbp,rho([pvtg.pos(2:end-1)'-1 end]),'-or', ...
    pbp, rho(pvtg.pos(1:end-1)), '-bo', pargs{:});
hold off
xlabel('Pressure [bar]'); title('Density of gaseous phase [kg/m3]');

%% Plot muG
figure, hold on, set(gca,'FontSize',14)
for i=1:numel(pbp)
    ind = pvtg.pos(i):pvtg.pos(i+1)-1;
    plot(p(ind), muGd(ind), '-ok',...
        'MarkerSize',6,'MarkerFaceColor',[.5 .5 .5] );
end
plot(pbp,muGd([pvtg.pos(2:end-1)'-1 end]),'-or', ...
    pbp, muGd(pvtg.pos(1:end-1)), '-ob',pargs{:});
hold off
xlabel('Pressure [bar]');
title('Gas viscosity [cP]');


%% Interpolation/extrapolation of viscosity
% For sufficintly high pressure/gas-oil ratios, the interpolated viscosity
% values exhibit unphysical behavior 
[M, N]    = deal(51,51);
[rs,pp]    = meshgrid(linspace(0,f.rvSat(240*barsa),M), ...
              linspace(200,240,N)*barsa);
Rv        = reshape(f.rvSat(pp(:)), N, M);
isSat     = rs >= Rv;
rs(isSat) = Rv(isSat);
muG       = reshape(f.muG(pp(:), rs(:), isSat(:)),N,M);

figure;
contourf(pp/barsa,rs,convertTo(muG,centi*poise),21);
hold on
p    = pvtg.key(rldecode((1:numel(pvtg.key))',diff(pvtg.pos),1))/barsa;
for i=1:numel(pbp)
    ind = pvtg.pos(i):pvtg.pos(i+1)-1;
    plot(p(ind),pvtg.data(ind,1),'-ok',...
        'MarkerSize',6,'MarkerFaceColor',[.5 .5 .5]);
end
plot(pbp,rvd,'-o',pargs{:});
hold off
cb = colorbar; set([gca;cb],'FontSize',20); xlabel('Pressure [bar]');

%%
% Plot viscosity as function of pressure only for a larger pressure
% interval
[M, N]    = deal(11,21);
[rs,pp]    = meshgrid(linspace(0,f.rvSat(280*barsa),M), ...
              linspace(160,280,N)*barsa);
Rv        = reshape(f.rvSat(pp(:)), N, M);
isSat     = rs >= Rv;
rs(isSat) = Rv(isSat);
muG       = convertTo(reshape(f.muG(pp(:), rs(:), isSat(:)),N,M),centi*poise);

figure
y = muG; y(isSat)=NaN;
z = muG; z(~isSat) = NaN;
plot(pp/barsa, y,'k--');
hold on
plot(pp/barsa, z, '--k','LineWidth',3);
hold on, set(gca,'FontSize',20)
for i=1:numel(pbp)
    ind = pvtg.pos(i):pvtg.pos(i+1)-1;
    plot(p(ind), muGd(ind), '-ok',...
        'MarkerSize',6,'MarkerFaceColor',[.5 .5 .5] );
end
plot(pbp,muGd([pvtg.pos(2:end-1)'-1 end]),'-ok', ...
    pbp, muGd(pvtg.pos(1:end-1)), '-ok',pargs{:});
hold off
xlabel('Pressure [bar]');
axis([160 280 0.015 0.08]);