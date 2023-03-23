function [vertices,edges] = removeDuplicates(vertices,edges)
% Remove duplicated edges and vertices
% The first tag is kept
% TODO: keep track of the lost tags
%
% SYNAPSES
% [vertices,edges] = removeDuplicates(vertices,edges)
%
% PARAMETERS
%       vertices - original vertices list
%       edges - original edge list
%
% OUTPUT
%       vertices - vertices list without duplicates
%       edges - edge list without duplicates
%
% Copyright 2011-2012 University of Bergen, 2013 IRIS AS
%
% This file is licensed under the GNU General Public License v3.0.

% take care of the tags
tags = edges(:,3:end);

% remove tags
edges = edges(:,1:2);

% remove duplicated nodes
[vertices,~,J] = unique(vertices,'rows','first');
edges(:,1:2) = J(edges(:,1:2));

% remove duplicated edges
[~,I] = unique(sort(edges,2),'rows','first');
edges = edges(I,:);

% keep the first tag
tags = tags(I,:);

% remove edges going nowhere
keep = edges(:,1)~=edges(:,2);
edges = edges(keep,:);
tags = tags(keep,:);

% add back the tags
edges = [edges,tags];