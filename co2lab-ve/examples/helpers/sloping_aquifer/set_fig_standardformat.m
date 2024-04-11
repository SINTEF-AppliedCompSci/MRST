function set_fig_standardformat(fig, fig_title)
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

if ~iscell(fig_title)
    fig_title = {fig_title};
end

plots = get(fig, 'children');
set(gcf, 'color', 'white');
N = numel(plots);
for i = 1:N
    axis tight;
    pix = N-i+1;
    set(plots(pix), 'fontsize', 14);
    
    set(get(plots(pix), 'title'), 'String',  fig_title{i});
    set(get(plots(pix), 'xlabel'), 'String', 'meter');
    set(get(plots(pix), 'ylabel'), 'String', 'meter');
    set(get(plots(pix), 'zlabel'), 'String', 'meter');
end


if N == 1
    set(fig, 'position', [420, 70, 1200, 360]);
elseif N == 2
    set(fig, 'position', [300, 300, 1200, 400*numel(plots)])
else
    % multiple plots, we need larger figure
    set(fig, 'position', [300, 300, 2000, 800])
end
