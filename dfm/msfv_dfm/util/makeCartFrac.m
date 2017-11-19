function makeCartFrac(n,ftag,box,filename)
% saves a Cartesian grid in filename
%
% SYNOPSIS
%  makeCartFrac(n,ftag,box,filename)
%
% PARAMETERS
% n:    a matrix containing the number of internal points for each set of
%       Cartesian grid i.e. different tags can be used for the different
%       sets
% ftag: tag numbers must have the size of the rows of n.
% box:  bounding box
% filename: file name where the vertices and edges are stored
%
% Copyright 2013 IRIS AS
%
% This file is licensed under the GNU General Public License v3.0.



% number of internal points for each set
numf = sum(n,1);

% find the endpoints of the sets
x = [];
y = [];
tags = [];

for i = 1 : numel(ftag)
    x_tmp = linspace(box(1,1),box(2,1),n(i,1)+2)';
    x_tmp = x_tmp(2:n(i,1)+1);
    y_tmp = linspace(box(1,2),box(2,2),n(i,2)+2)';
    y_tmp = y_tmp(2:n(i,2)+1);
    tags_tmp = ftag(i)*ones(n(i),1);
    x = [x ; x_tmp];
    y = [y ; y_tmp];
    tags = [tags ; tags_tmp];
end

% compute vertices and edges
vertices = [rldecode(x,2),repmat([box(1,2);box(2,2)],numf(1),1);  ...
    repmat([box(1,1);box(2,1)],numf(2),1),rldecode(y,2)];
edges = reshape((1:sum(numf)*2)',2,[])';
tags = repmat(tags,2,1);
edges = [edges,tags];

% store them in "filename"
save(filename,'vertices','edges','box');