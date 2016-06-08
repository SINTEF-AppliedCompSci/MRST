clc; clear all; close all;

addpath('../VEM2D/')
addpath('../')

load('basisGrid');

figure;
fill(P(:,1), P(:,2), 'y', 'facealpha', .2)
hold on
Kc = [.5,.45];
Ec = (P + P([2:end,1],:))/2;
nN = size(P,1);
nE = nN;

hold on;
for i = 1:nN
    h1 = plot(P(:,1), P(:,2), 'ok', 'MarkerFaceColor', [0 0.4470 0.7410], 'markersize', 6);
end
for i = 1:nE
    h2 = plot(Ec(:,1), Ec(:,2), 'sk', 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'markersize', 6);
end

h3 = plot(Kc(:,1), Kc(:,2), 'dk', 'MarkerFaceColor', 'w', 'markersize', 6);


h = legend([h1, h2, h3], '$\mathcal{V}^K$', '$\mathcal{E}^K$',...
                         '$\mathcal{P}^K$');
set(h, 'interpreter', 'latex')
set(h, 'fontsize', 14);
axis equal off;

%%

w = 0;
h = 2;
ps = get(gcf, 'Position');
ratio = 1;
paperWidth = 10;
paperHeight = paperWidth*ratio;
set(gcf, 'paperunits', 'centimeters');
set(gcf, 'papersize', [paperWidth-w paperHeight-h]);
set(gcf, 'PaperPosition', [-w    -h   paperWidth+w paperHeight+h]);
print(gcf, '-dpdf', '../../tex/thesis/fig/BasisElementDofs2D.pdf');