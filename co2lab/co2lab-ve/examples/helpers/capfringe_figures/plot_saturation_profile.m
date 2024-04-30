function plot_saturation_profile(s, seff, smax, G)
%Undocumented Utility Function

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

    % define colors
    red = [1 0 0]; blue = [0 0 1]; white = [1 1 1]; black = [0 0 0];
    mix = @(a, b, c) (a * c + b * (100 - c)) / 100;
    trapped_color = mix(mix(red, white, 40), black, 80);;
    plume_color = mix(red, white, 40);
    residual_color = mix(mix(red, blue, 50), white, 20);
    brine_color = mix(blue, white, 30);
    
    zvals = G.cells.centroids(:,3);
    hold on;
    % current saturation, colored field
    fill([s; s * 0], [zvals; flipud(zvals)], plume_color); 
    % immobilized saturation (current - effective saturation)
    fill([s - seff; 0*s], [zvals; flipud(zvals)], residual_color);
    % brine zone, colored field
    fill([s; s * 0 + 1], [zvals; flipud(zvals)], brine_color); 
    
    % historically maximal saturation
    last_ix = find(smax>0, true, 'last');
    last_ix = min(last_ix+1, numel(value(smax)));
    plot(smax(1:last_ix), zvals(1:last_ix), 'linewidth', 2, 'color', 'r');
    set(gca, 'ydir', 'reverse')

    xlabel('Saturation');
    ylabel('Vertial depth');
    
    
    %% in-graph annotations
    ylim = get(gca, 'ylim');
    % smax 
    z1 = zvals(last_ix);
    z2 = zvals(find(seff == 0, true, 'first'));
    zmid = (z1+z2)/2;
    xend = interpTable(zvals, smax, zmid);
    [X, Y] = dsxy2figxy_new([0.0, xend], ylim(2) - [zmid, zmid]);
    annotation('doublearrow', X, Y);
    text(xend/2 - 0.03, zmid + 1.1, 's_{max}', 'fontsize', 12);
    
    % s reidual
    zmid = 0.25 * z1 + 0.75 * z2;
    xend = interpTable(zvals, s, zmid);
    [X, Y] = dsxy2figxy_new([0.0, xend], ylim(2) - [zmid, zmid]);
    annotation('doublearrow', X, Y);
    text(xend/2 - 0.03, zmid + 1.1, 's_{n, r}', 'fontsize', 12);
    
    % s effective
    zmid = 0.5 * z2;
    xend = interpTable(zvals, s, zmid);
    xstart = interpTable(zvals, s - seff, zmid);
    [X, Y] = dsxy2figxy_new([xstart, xend], ylim(2) - [zmid, zmid]);
    annotation('doublearrow', X, Y);
    text(xend/2 - 0.03, zmid + 1.1, 's_{eff}', 'fontsize', 12);
end
