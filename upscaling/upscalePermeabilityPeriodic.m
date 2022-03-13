function perm = upscalePermeabilityPeriodic(Gp, bcp, dp_scale, psolver, fluid, rock, L)
% function for doing permeability upscaling on periodic grids
% perm = upscalePermeability(Gp,G,bcp,dp_scale)
%
%  Gp       - periodic grid
%  bcp      - boundary conditions
%  dp_scale - scale of pressure drop
%  Transp   - transmissibilities on periodic grid
%  fluid - one phase fluid

% if size(fluid.properties(),2)
%     error('Fluid must be a single phase fluid!')
% end

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


d      = Gp.griddim;
dfaces = cell(d,1);
vmat   = nan(d);

for j=1:d
    dfaces{j}=bcp.face(bcp.tags==j);
end

state   = initResSol(Gp, 100*barsa, 1);
dp_mat  = dp_scale*eye(d);
bcp_new = bcp;

for i=1:d;
    for j=1:d;
        bcp_new.value(bcp.tags==j)=dp_mat(j,i);%/L(j);
    end
    state  = psolver(state,Gp,fluid,bcp_new,rock);

    for j=1:d;
        flux=sum(state.flux(dfaces{j}).*bcp.sign(bcp.tags==j));
        vmat(j,i)=flux/sum(Gp.faces.areas(dfaces{j}));
    end
end

for j=1:d;
    dp_mat(j,:)=dp_mat(j,:)/L(j);
end

perm=-vmat/dp_mat;
end

