%{
G    = model.G;
rock = model.rock;
pos = [541 476 832 330];
wcol = [0.8, 0.8, 0.8];
%%
paras = {rock.perm(:,1)/(milli*darcy), state0.s};
%%
ii = 1;
figure, hold on
plotCellData(G, paras{ii})
view(0, 0)
axis on
colorbar
set(gcf, 'position', pos)
colormap(jet)


% 2500 2525 2550 2575 2600 2625
xx = [100, 100];
yy = [-2, -2];
zz = [2580, 2620];
plot3(xx, yy,zz, 'linewidth', 5, 'color', wcol)

zz = [2490, 2580];
plot3(xx, yy,zz, 'linewidth', 1, 'color', wcol)
text(xx(1)+50, yy(1), zz(1), 'I', 'color', 'r', 'FontName', 'Times','fontsize', 16)


xx = [3900, 3900];
yy = [-2, -2];
zz = [2505, 2545];
plot3(xx, yy,zz, 'linewidth', 5, 'color', wcol)

zz = [2490, 2505];
plot3(xx, yy,zz, 'linewidth', 1, 'color', wcol)
text(xx(1)+50, yy(1), zz(1), 'P', 'color', 'r', 'FontName', 'Times','fontsize', 16)

axis off
%%
figure
plotCellData(G, state0.s)
view(0, 0)
axis off
colorbar
set(gcf, 'position', pos)
%}
%%
G    = model.G;
%%
report = {reportW, reportP, reportS, reportsSP, reportB};
wellSols = {wellSolsW, wellSolsP, wellSolsS, wellSolsSP, wellSolsB};

T = cellfun(@(x)x.ReservoirTime / day, report, 'UniformOutput', false);
%%
casenanme = {'SURFACTANTPOLYMER2D_NOSP',  ...
    'SURFACTANTPOLYMER2D_P', 'SURFACTANTPOLYMER2D_S', 'SURFACTANTPOLYMER2D'};

prefixs = cellfun(@(x)['/home/xinsun/Desktop/Eclipse SP results/', x], casenanme, 'UniformOutput', false);

states_Ecl = cellfun(@(x)convertRestartToStates(x, G), prefixs, 'UniformOutput', false);

T_Ecl = cellfun(@(y)cellfun(@(x)x.time/day, y(2:end)), states_Ecl, 'UniformOutput', false);
wellSols_Ecl = cellfun(@(y)cellfun(@(x)x.wellSol, y(2:end), 'UniformOutput', false), states_Ecl, 'UniformOutput', false);
%%
% w = 1;
% field = 'bhp';
% unit = (mega*Pascal);

w = 2;
field = 'qOs';
unit = -(meter^3/day);
%%
para = cellfun(@(y) cellfun(@(x)x(w).(field)/ unit, y),  wellSols, 'UniformOutput', false);
para_Ecl = cellfun(@(y) cellfun(@(x)x(w).(field)/ unit, y),  wellSols_Ecl, 'UniformOutput', false);
%%
% para = cellfun(@(y) cellfun(@(x)x(w).qWs / (x(w).qOs + x(w).qWs) , y),  wellSols, 'UniformOutput', false);
% para_Ecl = cellfun(@(y) cellfun(@(x)x(w).qWs / (x(w).qOs + x(w).qWs), y),  wellSols_Ecl, 'UniformOutput', false);
%%
figure, hold on
cols = {'r', 'g', 'b', 'k'};
xx = [1,2,3,4];
arrayfun(@(x)plot(T{x}, para{x}, '-', 'color', cols{x}), xx)
arrayfun(@(x)plot(T_Ecl{x}, para_Ecl{x}, '--', 'color', cols{x}), xx)
%%
figure, hold on
plot(T{1}, para{1}, 'r-')
plot(T{5}, para{5}, 'g-')
plot(T_Ecl{1}, para_Ecl{1}, 'b-')


%%
figure, hold on
s_Ecl = states_Ecl(4);
plot(statesSP{228,1}.s(:,1), 'r-')
plot(s_Ecl{1,1}{229,1}.s(:,1), 'ro-')
plot(statesSP{228,1}.s(:,2), 'g-')
plot(s_Ecl{1,1}{229,1}.s(:,2), 'go-')
plot(statesSP{228,1}.s(:,3), 'b-')
plot(s_Ecl{1,1}{229,1}.s(:,3), 'bo-')