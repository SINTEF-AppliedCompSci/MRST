function [uu, extra] = LSMIM_linElast(G, C, el_bc, load, varargin)
%Undocumented Utility Function

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

opt = struct('solvetype','direct');
opts =  merge_options(opt,varargin{:});
opt=struct('stabterm','inter','as_type','CN_type','invertBlocks', 'mex');
tic;
lsmim = localStressMimetic(G, C,'theta',0.3,'invertBlocks',opt.invertBlocks,'stabterm',opt.stabterm,'as_type',opt.as_type);
vec=[1 0];
%profile viewer

mm='direct';
switch mm
    case 'mixed'
        ngam=size(lsmim.mixed.gama1,1);
        rhs=[zeros(size(lsmim.mixed.B,1),1);...
            reshape(bsxfun(@times,load(G.cells.centroids),G.cells.volumes)',[],1);...
            zeros(ngam,1)];      
        mm1=[lsmim.mixed.CT',lsmim.mixed.gama1'];
        mm2=[lsmim.mixed.CT',lsmim.mixed.gama2'];
        A=[lsmim.mixed.B,mm1;mm2',sparse(size(mm1,2),size(mm1,2))];
        u=A\rhs;
        utmp=u(size(lsmim.mixed.B,1)+1:size(lsmim.mixed.B,1)+G.cells.num*G.griddim);
        uu=-reshape(utmp,2,[])';
        if(ngam==G.nodes.num)
            figure(),plotNodeData(G,u(end-ngam:end)),colorbar
        end
    case 'direct'
        u=lsmim.direct.A\reshape(bsxfun(@times,load(G.cells.centroids),G.cells.volumes)',[],1);
        uu=reshape(u,2,[])';
    otherwise 
        error()
end

extra=lsmim;
end
