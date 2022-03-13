function [ep,ef] = L2_Norm(G,ua,fa,uh,fh)
%Compute L2 norm of pressure errors and flux errors. used for convergence
%study
%   G - Grid structure of MRST
%   ua - analytical pressure solution at cell centroids
%   fa - analytical flux solution at face centroids
%   uh - numerical pressure solution at cell centroids
%   fh - numerical flux solution at face centroids
%   ep - L2 norm of pressure error
%   ef - L2 nomr of flux error

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

A=G.cells.volumes;
ep=sqrt(sum(A.*(ua-uh).^2)/sum(A));

fa=abs(fa./G.faces.areas);
fh=abs(fh./G.faces.areas);
A=zeros(G.faces.num,1);
c1=G.faces.neighbors(:,1);
c2=G.faces.neighbors(:,2);
ind=all(G.faces.neighbors~=0,2);
A(ind)=0.5*(G.cells.volumes(c1(ind))+G.cells.volumes(c2(ind)));
A(~ind)=G.cells.volumes(max(G.faces.neighbors(~ind,:),[],2));
ef=sqrt(sum(A.*(fa-fh).^2)/sum(A));

% A=G.cells.volumes;
% ep=sqrt(sum(A.*(ua-uh).^2)/sum(A.*ua.^2));
% fa=abs(fa./G.faces.areas);
% fh=abs(fh./G.faces.areas);
% A=zeros(G.faces.num,1);
% c1=G.faces.neighbors(:,1);
% c2=G.faces.neighbors(:,2);
% ind=all(G.faces.neighbors~=0,2);
% A(ind)=0.5*(G.cells.volumes(c1(ind))+G.cells.volumes(c2(ind)));
% A(~ind)=G.cells.volumes(max(G.faces.neighbors(~ind,:),[],2));
% ef=sqrt(sum(A.*(fa-fh).^2)/sum(A.*fa.^2));
end
