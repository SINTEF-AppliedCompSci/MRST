%% IGEMS Data Set
% The IGEMS project studied how top surface morphology influences the CO2
% storage capacity. Alternative top-surface morphologies were created
% stochastically by combining different stratigraphic scenarios with
% different structural scenarios. In this example, we will use one of the
% top surfaces developed in the project to demonstrate capacity estimates
% and spill-point analysis on a model with a huge number of cells. For
% depositional features, two scenarios were chosen for which it was
% considered likely that a depositional/erosional topography could be
% preserved under a thick regional seal; the latter commonly formed by
% marine shale. The two scenarios reflect situations where sand deposition
% is suceeded by deposition of fines as a result of marine transgression:
% 
% * Offshore sand ridges covered by thick marine shale (OSS)
% * Preserved beach ridges under marine shale (FMM)
%
% In this example, we will consider an OSS-type surface that consists of
% sand dunes with amplitude up to 20 meters, width 2-4 km, and length 10-60
% km, and spacing 2-4 km. In addition, there is a conceptual structural
% scenario that consists of a single fault system with random 20-150 m
% displacement, random 300-6000 m length, and 90 degrees strike.
%
% More information about IGEMS
% Link: <http://www.nr.no/nb/IGEMS>
% Data description: <http://files.nr.no/igems/data.pdf>

mrstModule add co2lab
%% Read and prepare grid structure
% The project developed both full 3D grid saved in the ECLIPSE format and
% surfaces saved in the IRAP formate.  The ECLIPSE files are huge (588 MB)
% and reading and processing them typically requires a computer with at
% least 12 GB memory. Here, we will therefore only use the surface grid.
Gt = topSurfaceGrid( readIGEMSIRAP('OSSNP1', 1) );

%% First study: geometric analyisis of caprock (spill point analysis)
% The following command carries out a structural trapping analysis based on
% the geometry of the top surface.  It identifies _all traps_ (local
% pockets), along with their _depth_ (at which they 'spill over'),
% _connections_ (i.e., when a trap spills over, which trap(s) does the CO2
% flow into next), as well as the _rivers along which CO2 flows between
% traps_.  In addition, each cell in the top surface grid is associated
% with an unique _spill region_ associated with a trap (or with the
% exterior).  A spill region of a given trap consists of all the cells in
% the grid that 'leads into' the trap; in other words, it consists of all
% the cells from which a quantity of CO2 would eventually flow into the
% trap.

ts = trapAnalysis(Gt, false); % 'true' runs a cell-based algorithm for analysis;
                              % 'false' a node-based one (but result still
                              % presented for cells, not nodes).

%%
% The resulting structure contains the following fields:
ts %#ok

%% 
% |ts.traps| associates an integer to each cell of the top grid.  For cells
% located _within_ a trap, the integer will be the index of that trap.  For
% other cells, the integer will be a zero.  Traps are indexed from 1 and
% upwards. The total number of traps would thus be given by:
num_traps = max(ts.traps) %#ok

%%
% To have an idea about the location and distribution of traps, we can plot
% the trap cells on the grid by generating a scalar field over the grid
% cells, with two possible values ('trap-cell' and 'non-trap-cell') and
% then call the 'plotCellData' command.  On the resulting plot, we can see
% several long-narrow traps aligned with the crests of the surface, as well
% as a large number of scattered, small pockets.
figure; p = get(gcf,'Position'); set(gcf,'Position', p + [0 -300 0 300]);
trap_field = zeros(size(ts.traps));
trap_field(ts.traps>0) = 2;
plot_opts = {'EdgeColor', 'k', 'EdgeAlpha', 0.1};
plotCellData(Gt, trap_field, plot_opts{:});
view(-65,25); axis tight, colormap('jet');

%%
% Likewise, the 'rivers' exiting each trap (and then either entering
% another trap or exiting the domain) are stored in |ts.cell_lines|.  We
% can visualise them together with the traps by making a separate scalar
% field for rivers:

river_field = zeros(size(ts.traps));
for r = [ts.cell_lines{:}]'
    for c = 1:numel(r);
        river_field(r{c}) = 1;
    end
end
 
clf;
plot_opts = {'EdgeColor','none'};
plotCellData(Gt, max(trap_field, river_field), plot_opts{:});
view(-90,90); axis equal tight


%% 
% The vector |ts.trap_z| contains the _spill point depth_ for each trap,
% i.e., the depth at which the trap 'spills over'.  The trapping capacity
% of a trap is thus limited to the volume between the top of the grid
% within that trap and the plane defined by z equal to the spill point
% depth of that trap.  (In addition, one would have to multiply by the rock
% porosity, since the storage volume is limited to the pore volume.)
% 
% With |ts.traps|, |ts.trap_z| and the porosity values from our averaged
% rock structure, we can now easily compute the storage volumes for each
% trap.

poro = .2*ones(Gt.cells.num,1);
trap_volumes = volumesOfTraps(Gt, ts, 1:num_traps, 'poro', poro);

total_trapping_capacity = sum(trap_volumes);
fprintf('Total trapping capacity is: %6.2e\n\n', total_trapping_capacity);

%%
% We can get some information on the 10 largest traps of the grid:
[sorted_vols, sorted_ix] = sort(trap_volumes, 'descend');

fprintf('trap ix   | trap vol(m3)  | cells in trap\n');
fprintf('----------+---------------+--------------\n');
for i = 1:10
   fprintf('%7d   |  %10.3e   | %5d\n', ...
       sorted_ix(i), sorted_vols(i), numel(find(ts.traps == sorted_ix(i))));
end
fprintf(['\nTogether, these traps cover %6.2e m3, which represents %3.1f ' ...
   'percent of\nthe total trapping capacity of this grid.\n'], ...
   sum(sorted_vols(1:10)), sum(sorted_vols(1:10)) / total_trapping_capacity * 100)

%%
% Now, we visualise the location of these 10 traps (plotted in red, against
% the remaining traps in yellow):
largest_traps_field = zeros(size(ts.traps));
largest_traps_field(ismember(ts.traps, sorted_ix(1:10))) = 3;
clf; 
plotCellData(Gt, max(trap_field, largest_traps_field), plot_opts{:});
view(-90,90);axis equal tight

%%
% We can also color code each trap according to the total volume it holds:
for i = 1:num_traps
    trap_field(ts.traps == i) = trap_volumes(i);
end
clf; 
plotCellData(Gt, trap_field, plot_opts{:});
axis equal tight; view(-90,90); colorbar;


%%
% Each trap has an associated _accumulation region_, which consists of all
% cells 'spilling into' that trap.  In other words, injection of a quantity
% of CO2 into one of these cells would, assuming a purely gravity-driven
% flow, lead it to flow into the trap.  If we also count all cells leading
% _out_ of the grid as a separate accumulation region, then any grid cell
% belongs to exactly one accumulation region (either into a specific trap
% or out of the domain.
%
% The accumulation region a cell belongs to is given by |ts.trap_regions|.
% There is one integer value per cell, either indicating the index of the
% corresponding trap, or '0' if it belongs to the accumulation region that
% leads out of the grid domain.  We can visualize these regions by plotting
% them as a field on the grid, and use a colormap with sharp variations in
% color:

clf; plotCellData(Gt, ts.trap_regions, trap_field==0,plot_opts{:});
plotGrid(Gt, trap_field>0, 'FaceColor', 'k', 'EdgeColor','none');
nreg = max(ts.trap_regions);
colormap((colorcube(nreg+1)+2*ones(nreg+1,3))/3); 
view(-90,90); axis equal tight

%%
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
