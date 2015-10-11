%% Example 1: Cartesian grid
[nx,ny] = deal(15,10);
[x,y] = meshgrid(linspace(0,1,nx+1),linspace(0,1,ny+1));
p = [x(:) y(:)];
n = (nx+1)*(ny+1);
I = reshape(1:n,ny+1,nx+1);
T = [
    reshape(I(1:end-1,1:end-1),[],1)';
    reshape(I(1:end-1,2:end  ),[],1)';
    reshape(I(2:end,  2:end  ),[],1)';
    reshape(I(2:end,  1:end-1),[],1)'
    ]';

G = tessellationGrid(p, T);
clf, plotGrid(G);

%%
[dx, dy, dPhi] = deal(cos(pi/3),sin(pi/3), pi*40/180);
v = pi/180*[0 180 120 -60 240 60]';
dv = [cos(v-dPhi) sin(v-dPhi) cos(v+dPhi) sin(v+dPhi)]/3.5;
P = [...
    0           0           0           0;
    0+dv(1,1)   0+dv(1,2)   0+dv(1,3)   0+dv(1,4);
    1+dv(2,1)   0+dv(2,2)   1+dv(2,3)   0+dv(2,4);
    1           0           1           0;
    1+dv(3,1)   0+dv(3,2)   1+dv(3,3)   0+dv(3,4);
    dx+dv(4,1)  dy+dv(4,2)  dx+dv(4,3)  dy+dv(4,4);
    dx          dy          dx          dy;
    dx+dv(5,1)  dy+dv(5,2)  dx+dv(5,3)  dy+dv(5,4);
     0+dv(6,1)  0+dv(6,2)   0+dv(6,3)   0+dv(6,4);
    ];
m  = size(P,1);
P1 = P(:,1:2);
P2 = [P([1 m:-1:2],3) -P([1 m:-1:2],4)];
T  = reshape(1:4*m,m,4)';

[p,t,n] = deal([],[],0);
for j=0:2
    for i=0:4
        p = [p; bsxfun(@plus,P1,[i 2*j*dy])];                              %#ok<AGROW>
        p = [p; bsxfun(@plus,P2,[i-dx (2*j+1)*dy])];                       %#ok<AGROW>
        p = [p; bsxfun(@plus,P1,[i-dx (2*j+1)*dy])];                       %#ok<AGROW>
        p = [p; bsxfun(@plus,P2,[i 2*(j+1)*dy])];                          %#ok<AGROW>
        t = [t; T+n]; n=n+4*m;                                             %#ok<AGROW>
    end
end

[p,ia,ic] = unique(round(p*1e5)/1e5,'rows');
G = tessellationGrid(p, ic(t));
i=repmat((1:2)',G.cells.num/2,1);
clf
plotCellData(G,i(:));
plotFaces(G,find(any(G.faces.neighbors==0,2)),'EdgeColor','r','LineWidth',2);
axis tight off;

%%
[dx, dy] = deal(cos(pi/3),sin(pi/3));
P1 = [0 0; 1 0; dx dy];
P2 = [P1([1 3 2],1) -P1([1 3 2],2)];
[p,t,n] = deal([],[],0);
for j=0:2
    for i=0:4
        p = [p; bsxfun(@plus,P1,[i 2*j*dy])];                              %#ok<AGROW>
        p = [p; bsxfun(@plus,P2,[i-dx (2*j+1)*dy])];                       %#ok<AGROW>
        p = [p; bsxfun(@plus,P1,[i-dx (2*j+1)*dy])];                       %#ok<AGROW>
        p = [p; bsxfun(@plus,P2,[i 2*(j+1)*dy])];                          %#ok<AGROW>
        t = [t; [1:3; 4:6; 7:9; 10:12]+n]; n=n+12;                         %#ok<AGROW>
    end
end

dPhi   = pi*45/180;
e      = reshape(t(:,[1 2 2 3 3 1])',2,[])';
v      = p(e(:,2),:) - p(e(:,1),:);
phi    = atan2(v(:,2),v(:,1));
P      = p;

pn     = p(e(:,1),:) + 1/3*[cos(phi-dPhi)  sin(phi-dPhi)];
e      = e(:,[1 1 2]);
e(:,2) = size(P,1)+(1:size(pn,1)).';
P      = [P; pn];

pn     = p(e(:,1),:) + 1/3*[cos(phi-.2*dPhi)  sin(phi-.2*dPhi)];
e      = e(:,[1:3 3]);
e(:,3) = size(P,1)+(1:size(pn,1)).';
P      = [P; pn];

phi    = atan2(-v(:,2),-v(:,1));
pn     = p(e(:,4),:) + 1/3*[cos(phi-.2*dPhi)  sin(phi-.2*dPhi)];
e      = e(:,[1:4 4]);
e(:,4) = size(P,1)+(1:size(pn,1)).';
P      = [P; pn];

e      = e(:,[1:5 5]);
pn     = p(e(:,5),:) + 1/3*[cos(phi-dPhi) sin(phi-dPhi)];
e(:,5) = size(P,1)+(1:size(pn,1)).';
P      = [P; pn];

T         = reshape(e',18,[])';
T         = T(:,[1:5 7:11 13:17]);
[P,~,ic] = unique(round(P*1e5)/1e5, 'rows','stable');
G         = tessellationGrid(P,ic(T)); 

i=repmat((1:2)',G.cells.num/2,1);
clf
plotCellData(G,i(:));
plotFaces(G,find(any(G.faces.neighbors==0,2)),'EdgeColor','r','LineWidth',2);
axis tight off;

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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
