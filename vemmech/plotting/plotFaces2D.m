function plotFaces2D(G, faces, varargin)
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

opt=struct('col','r','uu',[],'shift',[],'linewidth',4);
opt=merge_options(opt,varargin{:});
if(~isempty(opt.uu))
    G.nodes.coords=G.nodes.coords+opt.uu;
end
opt=merge_options(opt,varargin{:});
if(G.griddim==2)
      nodes=G.faces.nodes(mcolon(G.faces.nodePos(faces),G.faces.nodePos(faces+1)-1));
      nodes=reshape(nodes,2,[])';
      
      for j=1:size(nodes,1)
          %line(G.nodes.coords(nodes(1,:),1),G.nodes.coords(nodes(j,:),2),varargin{:})
          pos1=G.nodes.coords(nodes(j,:),1);
          pos2=G.nodes.coords(nodes(j,:),2);
          if(~isempty(opt.shift))
          pos1=pos1+bsxfun(@times,G.faces.normals(faces(j),1),opt.shift(j));
          pos2=pos2+bsxfun(@times,G.faces.normals(faces(j),2),opt.shift(j));
          end
          line(pos1,pos2,'Color',opt.col,'LineWidth',opt.linewidth);
      end
      return
else
    error('This function do only work for 2D')
end
