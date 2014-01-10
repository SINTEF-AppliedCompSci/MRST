function nodetype = nodeType(g,nodes,box,h)
% return node type
% 0: inner node 1: edge node 2: corner node
%
% SYNAPSES
% boundaryType = numBoundaries(g,nodes,box,h)
%
% PARAMETERS
%       g - mrst grid structure
%       nodes - a list of nodes to check
%
%       OPTIONAL
%       box  - bounding box of the domain
%       h - a distance tolerance to include close nodes
%           default: sqrt(eps)
%
% OUTPUT
%       nodeType - a vector giving the node types of the input nodes
%                0: inner node 1: edge node 2: corner node
%
% Copyright 2011-2012 University of Bergen, 2013 IRIS AS
%
% This file is licensed under the GNU General Public License v3.0.

% find the bounding box if it is not given
if nargin < 3
    maxb = max(g.nodes.coords);
    minb = min(g.nodes.coords);
    box = [minb(1),maxb(1),minb(2),maxb(2)];
end

% use a default tolerance if not given
if nargin < 4
    h = sqrt(eps);
end

% check whether the x coordinate lies on the boundary
c1 = isequal(g.nodes.coords(nodes,1),box(1),h);
c2 = isequal(g.nodes.coords(nodes,1),box(2),h);

% check whether the y coordinate lies on the boundary
c3 = isequal(g.nodes.coords(nodes,2),box(3),h);
c4 = isequal(g.nodes.coords(nodes,2),box(4),h);

% find the boundary type
nodetype =  sum([c1,c2,c3,c4],2);

end

% pick the close nodes
function is = isequal(a,b,h)
is = sum((a - b).^2,2)<h;
end