function s = systemStressMimetic(g, CC, varargin)
%Undocumented Utilty Function

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
      
       %s = computeMimeticIP(g,rock,varargin{:});
       dims=g.griddim;
       no2dof =@(no) reshape(bsxfun(@plus,g.griddim*(no-1),[1:g.griddim])',[],1);
       
       s = computeStressMimeticIP(g,CC,'type','hybrid');
       stmp= computeStressMimeticIP(g,CC,'type','mixed');
       s.B=stmp.B;
       cellNo  = rldecode(1:g.cells.num, diff(g.cells.facePos), 2) .';
       s.C     = sparse(no2dof(1:numel(cellNo)), no2dof(cellNo), 1);
       s.D     = sparse(no2dof(1:numel(cellNo)), no2dof(double(g.cells.faces(:,1))), 1, ...
           numel(cellNo)*dims, g.faces.num*dims);
       s.sizeB = repmat(size(g.cells.faces, 1), [1,2]);
       s.sizeC = size(s.C);
       s.sizeD = size(s.D);
end
