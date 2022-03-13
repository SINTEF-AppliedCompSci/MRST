function [C,N, CC] = maxTPFAGravityMatrix(Gtop,varargin)
%Returns a matrix defining the maximal gravity connection of the top surface
%
% SYNOPSIS:
%    C          = maxTPFAGravityMatrix(Gtop)
%    [C, N ]    = maxTPFAGravityMatrix(Gtop)
%    [C, N, CC] = maxTPFAGravityMatrix(Gtop)
%
% VARIABLES:
%    Gtop - grid structure for a top-surface grid
%
%    C    - matrix with maximal gravity connection
%
%    N    - 2*num_cells matrix defining the neibouring cells
%
%    CC   - matrix with maximal gravity connection including conections to
%           boundary. Boundary is defind as cell num_cells+1
%
%  'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%            supported options are:
%
%              use_multipoint -- locical to extend neighbours to nodebased
%                                stencile neighbours.
%
% DESCRIPTION:
%    The routine constructs a matrix that represents the maximal gravity
%    connection of a top-surface grid. To this end, we define a two-point
%    connection between cells (i.e., two cells are connected only if they
%    share a common face) and from this graph, we define a subgraph in
%    which there exists a directed edge between a cell and the neighbour
%    that has the largest height difference in cell centroids. In addition,
%    bi-directional edges are added to all neighbouring cells having a zero
%    height difference. The resulting graph defines an approximation to the
%    direction in which a light fluid will migrate at an infinitesimaly
%    slow rate.
%
% EXAMPLE:
%    To demonstrate how to use the maximal gravity connection, we will use
%    graph routines from the MatlabBGL to investigate the connections.
%
%    C = maxTPFAGravityMatrix(Gtop);
%
%    Find the accumulation or drainage area of each local maximum. To
%    this end, we replace the unidirectional connections by bi-directional
%    connections and find the strongly connected components of the
%    resulting undirected graph
%
%    [ci,ss] = components(C+C');
%
%    Find all cells that are in the accumulation/drainage area of a given
%    cell. To this end, we can either use a depth-first search or a
%    breadth-first search
%
%    cc=dfs(C',cell)
%    cc=bfs(C',cell)
%
%    Likewise, we can determine cells that are upstream of a given cell
%
%    cc=dfs(C,cell)
%    cc=bfs(C,cell)

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

opt = struct('use_multipoint',false);
opt = merge_options(opt, varargin{:});


if(opt.use_multipoint)
    N = neighboursByNodes(Gtop,'include_bc',false);
else
    N=Gtop.faces.neighbors;
end

intern=all(N>0,2);
ni1=N(intern,1);ni2=N(intern,2);
dist=sqrt(sum((Gtop.cells.centroids(ni2,:)-Gtop.cells.centroids(ni1,:)).^2,2));
z_diff   = Gtop.cells.z(ni2)-Gtop.cells.z(ni1);
z_diff   = z_diff./dist;

cnum   = Gtop.cells.num;
B      = sparse(double(ni2), double(ni1), z_diff, cnum, cnum);
B      = B + sparse(double(ni1), double(ni2), -z_diff, cnum, cnum);
[d,ii] = max(B, [],2);



indy   = ii(d>0);
indx   = 1:cnum;
indx   = indx(d>0);

C      = sparse(indx,indy,1,cnum,cnum);


C      = C + sparse(double(ni2),double(ni1),z_diff==0,cnum,cnum) + ...
    sparse(double(ni1),double(ni2),z_diff==0,cnum,cnum);
if(nargout>2)
    %
    CC      = sparse(indx,indy,1,cnum+1,cnum+1);
    CC      = CC + sparse(double(ni2),double(ni1),z_diff==0,cnum+1,cnum+1) + ...
        sparse(double(ni1),double(ni2),z_diff==0,cnum+1,cnum+1);
    
    b_cells=sum(N(~intern,:),2);
    ind=full(sum(CC(b_cells,:),2))==0;
    b_out=b_cells(ind);
    %
    CC      = CC + sparse(b_out,cnum+1,1,cnum+1,cnum+1);
end

end

