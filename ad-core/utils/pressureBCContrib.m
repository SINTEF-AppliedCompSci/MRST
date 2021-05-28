function  [bXqXbc,bc_cell] = pressureBCContrib(G,s,pX,rhoX,mobX,bX,bc)
% LEGACY FUNCTION: Intentionally undocumented.
    assert(all(strcmp(bc.type,'pressure')),'only pressure bc allowed');
    Tbc=s.T_all(bc.face);
    assert(all(sum(G.faces.neighbors(bc.face,:)>0,2)==1),'bc on internal boundary');
    bc_cell=sum(G.faces.neighbors(bc.face,:),2);        
    dzbc=(G.cells.centroids(bc_cell,end)-G.faces.centroids(bc.face,end));
    g=norm(gravity);
    bXqXbc=cell(numel(pX),1);
    for i=1:numel(pX)
        pXbc=pX{i}(bc_cell);
        rhoXbc=rhoX{i}(bc_cell);
        dpbc=bc.value-pXbc+g*(rhoXbc.*dzbc);
        bXmobXbc = bX{i}(bc_cell).*mobX{i}(bc_cell);
        bXmobXbc(dpbc>0)=0;
        %NB some definition of in flow has to be done
        %if(any(dpbc_w>0))
        %    bWmobWbc(dpbc_w>0)=bW(bc_cell(dpbc_w>0)).*(mobW(bc_cell(dpbc_w>0))+mobO(bc_cell(dpbc_w>0)));     
        %end
        bXqXbc{i}  = -bXmobXbc.*Tbc.*(dpbc);
    end      
end

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
