function model = spe10(domain)
%Define a model structure that represents a subset of the SPE10 base case.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

try
   load spe10_rock
catch ME %#ok
   make_spe10_data;
   rock      = SPE10_rock();
end

%%
di = domain(2)-domain(1)+1;
dj = domain(4)-domain(3)+1;
dk = domain(6)-domain(5)+1;

G              = cartGrid([di, dj, dk], [di,dj,dk].*[20, 10, 2]*ft);
offset         = (domain([1,3,5])-1).*[20,10,2]*ft;
G.nodes.coords = bsxfun(@plus, G.nodes.coords, offset);
G              = computeGeometry(G);
clear di dj dk offset

rock.perm      = convertFrom(rock.perm(:,:), milli*darcy);

[cartIx{1:3}]  = ndgrid(domain(1) : domain(2), ...
                        domain(3) : domain(4), ...
                        domain(5) : domain(6));
c              = sub2ind([60, 220, 85], ...
                          cartIx{1}(:), cartIx{2}(:), cartIx{3}(:));

rock.perm      = rock.perm(c, :);
rock.poro      = max(rock.poro(c), 1e-3);


tpv       = sum(rock.poro.*G.cells.volumes);

W = [];
dims = G.cartDims;
W = verticalWell(W, G, rock,  1,   1, (1:dims(3)),     ...
                     'Type', 'rate', 'Val', tpv/(50*year), ...
                     'Radius', 0.125, 'Name', 'I1','Comp_i',[1,0], ...
                     'InnerProduct', 'ip_tpf');
W = verticalWell(W, G, rock,  dims(1),   dims(2), (1:dims(3)),     ...
                     'Type', 'bhp', 'Val', 200*barsa, ...
                     'Radius', 0.125, 'Name', 'P1','Comp_i',[1,0], ...
                     'InnerProduct', 'ip_tpf');

T  = computeTrans     (G, rock);
S  = computeMimeticIP (G, rock);

%%
s        = zeros(G.cells.num, 1);
state    = initState (G, W, 0, [s, 1-s]);

model = struct ('grid',  G,    ...
                'rock',  rock, ...
                'T',     T,    ...
                'S',     S,    ...
                'W',     W,    ...
                'state', state);
