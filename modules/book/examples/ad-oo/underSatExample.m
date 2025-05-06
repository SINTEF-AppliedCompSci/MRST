%% Plotting for Sector Model Example
% This script contains plotting of results from blackoilSectorModelExample
% from the ad-blackoil module
mrstModule add ad-blackoil
blackoilSectorModelExample;

t = convertTo(cumsum(schedule.step.val),year);
dT = dt/day;
lw = {'LineWidth', 2};
ms = {'MarkerSize', 5};

%% Bottom-hole pressure and field-average pressure
clf, set(gcf,'Position',[800 200 640 420],'PaperPositionMode','auto');
set(gca,'FontSize',20);
bhp   = convertTo(cellfun(@(x) x.bhp, ws  ), barsa);
bhpc  = convertTo(cellfun(@(x) x.bhp, ws_c), barsa);
pavg  = convertTo(cellfun(@(x) mean(x.pressure), states  ), barsa);
pavgc = convertTo(cellfun(@(x) mean(x.pressure), states_c), barsa);
plot(t, bhp, '-', t, bhpc,'-o',lw{:}, ms{:});
hold on
plot(t, pavg, '--', t, pavgc, '--',lw{:});
plot([870 870]*day/year,[100 310],':k');
hold off; set(gca,'YLim',[90 310]);
xlabel('Time [years]');
title('Bottom-hole pressure [bar]');
legend('Pressure','Closed','Location','Best');
box off
% print -depsc2 sector-bhp.eps;

%% Surface oil rate
clf
set(gca,'FontSize',20);
qOs   = -convertTo(cellfun(@(x) x.qOs, ws  ), 1000*meter.^3/day);
qOsc  = -convertTo(cellfun(@(x) x.qOs, ws_c), 1000*meter.^3/day);
plot(t, qOs, '-', t, qOsc,'-o',lw{:}, ms{:});
hold on
plot([870 870]*day/year,[0 2.1],':k');
hold off
set(gca,'YLim', [0 2.2]);
xlabel('Time [years]');
title('Surface oil rate [10^3 m^3/day]');
legend('Pressure','Closed','Location','Best');
box off
% print -depsc2 sector-qOs.eps;

%% Cumulative surface rates
field = {'qOs','qWs','qGs'};
name = {'oil','water','gas'};
for i=1:numel(field)
    clf,
    set(gca,'FontSize',20);
    d  = convertTo(cumsum(-cellfun(@(x) x.(field{i}), ws  ).*dt),1e6*meter^3);
    dc = convertTo(cumsum(-cellfun(@(x) x.(field{i}), ws_c).*dt),1e6*meter^3);
    plot(t, d, '-', t, dc,'-o',lw{:}, ms{:});
    hold on
    plot([870 870]*day/year,[0 max(d)],':k');
    hold off
    axis tight
    xlabel('Time [years]');
    title(['Cumulative ', name{i}, ' production [10^6 m^3]']);
    legend('Pressure','Closed','Location','Best');
    box off
    disp(' ---- Press any key ----');
    pause
%    print('-depsc2',['sector-',field{i},'-cum.eps']);
end

%% Instantaneous reservoir rates
field = {'qOr','qWr','qGr'};
name = {'Oil','Water','Gas'};
for i=1:numel(field)
    clf,
    set(gca,'FontSize',20);
    d  = convertTo(-cellfun(@(x) x.(field{i}), ws  ),1e3*meter^3/day);
    dc = convertTo(-cellfun(@(x) x.(field{i}), ws_c),1e3*meter^3/day);
    plot(t, d, '-', t, dc,'-o',lw{:}, ms{:});
    hold on
    plot([870 870]*day/year,[0 max(d)],':k');
    hold off
    axis tight
    xlabel('Time [years]');
    title([name{i}, ' reservoir rate [10^3 m^3/day]']);
    legend('Pressure','Closed','Location','Best');
    box off
    disp(' ---- Press any key ----');
    pause
%    print('-depsc2',['sector-',field{i},'.eps']);
end

%% Plot liberation of free gas
clf, set(gcf,'Position',[200 360 1200 260],'PaperPositionMode','auto');
g = cartGrid([1 1 1], [1000, 1000, 100]);
for i=1:3
    subplot(1,3,i)
    plotGrid(g,'FaceColor','none','LineWidth',1); plotWell(G,W); view(3)
    sg = states_c{10+i}.s(:,3);
    plotCellData(G, sg, sg>5e-3,'EdgeColor','none');
    title([num2str(t(10+i)*year/day), ' days'],'FontSize',14);
    axis equal off
    axis([-5 1005 -5 1005 -5 105]);
    caxis([0 0.025]);
end
pos = get(gca,'Position'); 
cb=colorbar('Location','EastOutside');
set(cb,'Position',[0.92 0.25 0.02 0.6],'FontSize',14);
set(gca,'Position',pos)
for i=2:-1:1
    subplot(1,3,i),
    set(gca,'Position',get(gca,'Position')+[.1/i 0 0 0]);
end