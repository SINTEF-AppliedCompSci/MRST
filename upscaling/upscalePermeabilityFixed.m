function perm = upscalePermeabilityFixed(G,dp_scale,psolver,fluid_pure,rock,L)
% function for doing permeability upscaling on grids with fixed pressure conditions (pressure drop in one direction and noflow elsewhere)
%
% perm =
% upscalePermeabilityFixed(G,dp_scale,psolver,fluid_pure,rock,L)
%
%  G       - period grid
%  bcp      - boundaryConditions
%  dp_scale - scale of pressure drop
%  Transp   - transmissibilities on periodic grid
%  fluid_pure - one phase fluid
%{

Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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


% find boundary face combinations
if(G.griddim>1)
   bcl{1}=pside([],G,'LEFT',0);
   bcr{1}=pside([],G,'RIGHT',0);
   bcl{2}=pside([],G,'BACK',0);
   bcr{2}=pside([],G,'FRONT',0);
end
if(G.griddim>2)
   bcl{3}=pside([],G,'TOP',0);
   bcr{3}=pside([],G,'BOTTOM',0);
end
if(G.griddim==1)
   bcl{1}=struct('face',1,'value',dp_scale);
   bcr{1}=struct('face',G.faces.num,'value',dp_scale);
end

state = initResSol(G, 100*barsa, 1);

vmat = zeros(G.griddim, 1);
dp_mat = vmat;

for i=1:G.griddim
   bc = addBC([], bcl{i}.face, 'pressure', dp_scale);
   bc = addBC(bc, bcr{i}.face, 'pressure', 0);

   state  = psolver(state,G,fluid_pure,bc,rock);

   flux = sum(state.flux(bcr{i}.face));
   vmat(i) = (flux/sum(G.faces.areas(bcr{i}.face)));
   dp_mat(i) = dp_scale/L(i);
end

perm= vmat./dp_mat;
end
