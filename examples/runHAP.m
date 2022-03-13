% Test the effect of myRatio on the accuracy of NFVMs
clear;
close all;
mrstModule add ad-core mpfa ad-blackoil compositional ad-props mrst-gui nfvm

% Set mrstVerbose on to display number of HAPs corrected
mrstVerbose on

%%
Lx = 2;
Ly = 2;
iangle = 7 * pi / 16;
slope = tan(iangle);
tol = 1e-7;
maxIter = 500;
permconfig = 2;
switch permconfig
  case 1
    k1 = 1;
    k2x = 10;
    k2y = 10;
    k3 = 100;
    ia = 0.5645;
    a = [1.0, -4.8591, -0.9664];
    b = [-12.0414, -6.0699, -0.2837];
  case 2
    k1 = 1;
    k2x = 10;
    k2y = 1e3;
    k3 = 100;
    ia = 0.6142;
    a = [1.0, -0.4275, -0.7604];
    b = [-1.0546, 0.2142, -0.6495];
  case 3
    k1 = 1;
    k2x = 10;
    k2y = 1e5;
    k3 = 100;
    ia = 0.8866;
    a = [1.0, -0.0144, 0.7544];
    b = [-0.3706, 0.0022, -0.6564];
end
xdomain = [-0.5 * Lx, 0.5 * Lx, 0.5 * Lx, -0.5 * Lx, -0.5 * Lx];
ydomain = [-0.5 * Ly, -0.5 * Ly, 0.5 * Ly, 0.5 * Ly, -0.5 * Ly];
set(figure, 'color', 'w');
plot(xdomain, ydomain, 'k');
hold on
plot([0.5 * Ly / slope, -0.5 * Ly / slope], [0.5 * Ly, -0.5 * Ly], 'k')
plot([0, 0.5 * Lx], [0, 0], 'k'); axis image off
text([0.25 * Lx, -0.25 * Lx, 0.25 * Lx, 0.01 * Lx], [0.25 * Ly, 0, -0.25 * Ly, 0.05 * Ly], ...
     {'1', '2', '3', ' \it\beta=\rm7\pi/16'}, 'fontsize', 15); clear xdomain ydomain
%--------------------------------------------------------------------------
n = [8, 16, 32, 64, 128]';
h = zeros(size(n));
num_solvers = 5;
ep = zeros(length(n), num_solvers);
ef = zeros(length(n), num_solvers);

gravity reset off;
fluid = initSimpleADIFluid('cR', 1e-8/barsa, 'rho', [1, 1000, 100]);

for i_con = 1:length(n)
    % Grid
    nx = n(i_con);
    ny = n(i_con);
    [i_con, nx, ny]
    G = Grid_Incline(Lx, Ly, nx, ny, iangle);
    h(i_con) = mean(G.faces.areas);
    figure, plotGrid(G, 'facecolor', 'none');
    axis image;

    % Rock
    x = G.cells.centroids(:, 1);
    y = G.cells.centroids(:, 2);
    th = car2pol(x, y);
    ind1 = th >= 0 & th <= iangle;
    ind2 = th >= iangle & th <= pi + iangle;
    ind3 = th > pi + iangle;
    rock = makeRock(G, 0, 0.2);
    rock.perm = zeros(G.cells.num, 3);
    rock.perm(ind1, :) = repmat([k1, 0, k1], [sum(ind1), 1]);
    rock.perm(ind2, :) = repmat([k2x, 0, k2y], [sum(ind2), 1]);
    rock.perm(ind3, :) = repmat([k3, 0, k3], [sum(ind3), 1]);
    clear x y th ind1 ind2 ind3

    % BC
    % bc.face = boundaryFaces(G);
    % bc.type = repmat({'pressure'}, [numel(bc.face), 1]);
    % bc.value = cell(numel(bc.face), 1);
    % x = G.faces.centroids(bc.face, 1);
    % y = G.faces.centroids(bc.face, 2);
    % th = car2pol(x, y);
    % ind1 = th >= 0 & th <= iangle;
    % ind2 = th >= iangle & th <= pi + iangle;
    % ind3 = th >= pi + iangle;
    % bc.value(ind1) = repmat({@(x)(sqrt(x(1)^2 + x(2)^2))^ia * ...
    %                     (a(1) * cos(ia * car2pol(x(1), x(2))) + b(1) * ...
    %                      sin(ia * car2pol(x(1), x(2)))) + 10}, [sum(ind1), 1]);
    % bc.value(ind2) = repmat({@(x)(sqrt(x(1)^2 + x(2)^2))^ia * ...
    %                     (a(2) * cos(ia * car2pol(x(1), x(2))) + b(2) * ...
    %                      sin(ia * car2pol(x(1), x(2)))) + 10}, [sum(ind2), 1]);
    % bc.value(ind3) = repmat({@(x)(sqrt(x(1)^2 + x(2)^2))^ia * ...
    %                     (a(3) * cos(ia * car2pol(x(1), x(2))) + b(3) * ...
    %                      sin(ia * car2pol(x(1), x(2)))) + 10}, [sum(ind3), 1]);

    % MRST BC
    bf = boundaryFaces(G);
    x = G.faces.centroids(bf, 1);
    y = G.faces.centroids(bf, 2);
    th = car2pol(x, y);
    ind1 = th >= 0 & th <= iangle;
    ind2 = th >= iangle & th <= pi + iangle;
    ind3 = th >= pi + iangle;
    fcn1 = @(x)(sqrt(x(:,1).^2 + x(:,2).^2)).^ia ...
           .* (a(1) * cos(ia * car2pol(x(:,1), x(:,2))) ...
               + b(1) * sin(ia * car2pol(x(:,1), x(:,2)))) + 10;
    fcn2 = @(x)(sqrt(x(:,1).^2 + x(:,2).^2)).^ia ...
           .* (a(2) * cos(ia * car2pol(x(:,1), x(:,2))) ...
               + b(2) * sin(ia * car2pol(x(:,1), x(:,2)))) + 10;
    fcn3 = @(x)(sqrt(x(:,1).^2 + x(:,2).^2)).^ia ...
           .* (a(3) * cos(ia * car2pol(x(:,1), x(:,2))) ...
               + b(3) * sin(ia * car2pol(x(:,1), x(:,2)))) + 10;
    val1 = fcn1(G.faces.centroids(bf(ind1), :));
    val2 = fcn2(G.faces.centroids(bf(ind2), :));
    val3 = fcn3(G.faces.centroids(bf(ind3), :));
    bc = [];
    bc = addBC(bc, bf(ind1), 'pressure', val1);%, 'sat');%, [1 0]);
    bc = addBC(bc, bf(ind2), 'pressure', val2);%, 'sat');%, [1 0]);
    bc = addBC(bc, bf(ind3), 'pressure', val3);%, 'sat');%, [1 0]);
    clear x y th ind1 ind2 ind3

    % Source
    x = G.faces.centroids(:, 1);
    y = G.faces.centroids(:, 2);
    [th, r] = car2pol(x, y);
    ind1 = (th >= 0 & th <= iangle + 1e-3);
    ind2 = th > iangle & th <= pi + iangle + 1e-3;
    ind3 = th > pi + iangle;
    ur = zeros(G.faces.num, 1);
    uth = zeros(G.faces.num, 1);
    ur(ind1) = ia .* r(ind1).^(ia - 1) .* (a(1) .* (cos(ia .* th(ind1))) + b(1) .* (sin(ia .* th(ind1))));
    uth(ind1) = ia .* r(ind1).^ia .* (-a(1) .* sin(ia .* th(ind1)) + b(1) .* (cos(ia .* th(ind1))));
    ur(ind2) = ia .* r(ind2).^(ia - 1) .* (a(2) .* (cos(ia .* th(ind2))) + b(2) .* (sin(ia .* th(ind2))));
    uth(ind2) = ia .* r(ind2).^ia .* (-a(2) .* sin(ia .* th(ind2)) + b(2) .* (cos(ia .* th(ind2))));
    ur(ind3) = ia .* r(ind3).^(ia - 1) .* (a(3) .* (cos(ia .* th(ind3))) + b(3) .* (sin(ia .* th(ind3))));
    uth(ind3) = ia .* r(ind3).^ia .* (-a(3) .* sin(ia .* th(ind3)) + b(3) .* (cos(ia .* th(ind3))));
    vx = cos(th) .* ur - 1 ./ r .* sin(th) .* uth;
    vy = sin(th) .* ur + 1 ./ r .* cos(th) .* uth;
    vx(ind1) = -k1 .* vx(ind1);
    vy(ind1) = -k1 .* vy(ind1);
    vx(ind2) = -k2x .* vx(ind2);
    vy(ind2) = -k2y .* vy(ind2);
    vx(ind3) = -k3 .* vx(ind3);
    vy(ind3) = -k3 .* vy(ind3);
    fa = dot(G.faces.normals, [vx, vy], 2);
    hf_flux = fa(G.cells.faces(:, 1));
    hf2cn = gridCellNo(G);
    sgn = 2 * (G.faces.neighbors(G.cells.faces(:, 1), 1) == hf2cn) - 1;
    hf_flux = sgn .* hf_flux;
    q = accumarray(hf2cn, hf_flux);
    src = addSource([], 1:G.cells.num, q, 'sat', 1);
    clear x y th r ind1 ind2 ind3 vx vy ur uth hf_flux hf2cn sgn q

    % Physical model
    model = GenericBlackOilModel(G, rock, fluid, 'water', true, 'oil', false, 'gas', false);
    u0 = 10 * ones(G.cells.num, 1);
    s0 = 1;
    state0 = initResSol(G, u0, s0);
    % model = GenericBlackOilModel(G, rock, fluid, 'water', true, 'oil', true, 'gas', false);
    % state0 = initResSol(G, 1*barsa, [0, 1]);

    % Schedule
    dt = 1;
    schedule = simpleSchedule(dt, 'bc', bc, 'src', src);

    %----------------------------------------------------------------------
    [~, s0] = simulateScheduleAD(state0, model, schedule);
    s0 = s0{1};

    model_ntpfa = setNTPFADiscretization(model, 'myRatio', []);
    [~, sntp1] = simulateScheduleAD(state0, model_ntpfa, schedule);
    sntp1 = sntp1{1};

    model_ntpfa = setNTPFADiscretization(model, 'myRatio', 3);
    [~, sntp2] = simulateScheduleAD(state0, model_ntpfa, schedule);
    sntp2 = sntp2{1};

    model_ntpfa = setNTPFADiscretization(model, 'myRatio', 2);
    [~, sntp3] = simulateScheduleAD(state0, model_ntpfa, schedule);
    sntp3 = sntp3{1};

    model_ntpfa = setNTPFADiscretization(model, 'myRatio', 3);
    [~, sntp4] = simulateScheduleAD(state0, model_ntpfa, schedule);
    sntp4 = sntp4{1};

    % s0 = OnePhaseIncompMPFA(G, TransMPFAO(G, rock, bc), src);
    % u0 = 10 * ones(G.cells.num, 1);
    % interpFace = findHAP(G, rock, bc);
    % OSflux = findOSflux(G, rock, bc, interpFace);
    % sntp1 = PicardNTPFA(G, bc, src, OSflux, u0, tol, maxIter);
    % snmp1 = PicardNMPFA(G, bc, src, OSflux, u0, tol, maxIter);

    % interpFace = correctHAP(G, interpFace, 3);
    % OSflux = findOSflux(G, rock, bc, interpFace);
    % sntp2 = PicardNTPFA(G, bc, src, OSflux, u0, tol, maxIter);
    % snmp2 = PicardNMPFA(G, bc, src, OSflux, u0, tol, maxIter);

    % interpFace = correctHAP(G, interpFace, 2);
    % OSflux = findOSflux(G, rock, bc, interpFace);
    % sntp3 = PicardNTPFA(G, bc, src, OSflux, u0, tol, maxIter);
    % snmp3 = PicardNMPFA(G, bc, src, OSflux, u0, tol, maxIter);

    % interpFace = correctHAP(G, interpFace, 1);
    % OSflux = findOSflux(G, rock, bc, interpFace);
    % sntp4 = PicardNTPFA(G, bc, src, OSflux, u0, tol, maxIter);
    % snmp4 = PicardNMPFA(G, bc, src, OSflux, u0, tol, maxIter);
    % clear bc src OSflux u0 interpFace coin
    %----------------------------------------------------------------------
    ua = zeros(G.cells.num, 1);
    x = G.cells.centroids(:, 1);
    y = G.cells.centroids(:, 2);
    [th, r] = car2pol(x, y);
    ind1 = th >= 0 & th <= iangle;
    ind2 = th >= iangle & th <= pi + iangle;
    ind3 = th > pi + iangle;
    ua(ind1) = r(ind1).^ia .* (a(1) * cos(ia * th(ind1)) + b(1) * sin(ia * th(ind1)));
    ua(ind2) = r(ind2).^ia .* (a(2) * cos(ia * th(ind2)) + b(2) * sin(ia * th(ind2)));
    ua(ind3) = r(ind3).^ia .* (a(3) * cos(ia * th(ind3)) + b(3) * sin(ia * th(ind3)));
    ua = ua + 10; clear x y th r ind1 ind2 ind3
    [ep(i_con, 1), ef(i_con, 1)] = L2_Norm(G, ua, fa, s0.pressure, s0.flux);
    [ep(i_con, 2), ef(i_con, 2)] = L2_Norm(G, ua, fa, sntp1.pressure, sntp1.flux);
    [ep(i_con, 3), ef(i_con, 3)] = L2_Norm(G, ua, fa, sntp2.pressure, sntp2.flux);
    [ep(i_con, 4), ef(i_con, 4)] = L2_Norm(G, ua, fa, sntp3.pressure, sntp3.flux);
    [ep(i_con, 5), ef(i_con, 5)] = L2_Norm(G, ua, fa, sntp4.pressure, sntp4.flux);
    % [ep(i_con, 6), ef(i_con, 6)] = L2_Norm(G, ua, fa, snmp1.pressure, snmp1.flux);
    % [ep(i_con, 7), ef(i_con, 7)] = L2_Norm(G, ua, fa, snmp2.pressure, snmp2.flux);
    % [ep(i_con, 8), ef(i_con, 8)] = L2_Norm(G, ua, fa, snmp3.pressure, snmp3.flux);
    % [ep(i_con, 9), ef(i_con, 9)] = L2_Norm(G, ua, fa, snmp4.pressure, snmp4.flux);
end
clear a b i_con ia iangle k1 k2x k2y k3 Lx Ly maxIter n nx ny permconfig slope tol

%%
mylegend = {'TPFA', 'NTPFA(\itR^{''}=inf)', 'NTPFA(\itR^{''}=3R)', 'NTPFA(\itR^{''}=2R)', 'NTPFA(\itR^{''}=R)'};
plotConvergence(ep(:, 1:5), ef(:, 1:5), h, mylegend)
%mylegend = {'TPFA', 'NMPFA(\itR^{''}=inf)', 'NMPFA(\itR^{''}=3R)', 'NMPFA(\itR^{''}=2R)', 'NMPFA(\itR^{''}=R)'};
%plotConvergence(ep(:, [1, 6:end]), ef(:, [1, 6:end]), h, mylegend)
%clear mylegend


function [th, r] = car2pol(x, y)
    [th, r] = cart2pol(x,y);
    th = mod(th, 2*pi);
end

%% Copyright Notice
%
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
