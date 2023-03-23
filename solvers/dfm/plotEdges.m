function plotEdges(coords,edges,varargin)
% plot lines
% SYNOPSIS
%   plotEdges(coords,edges,varargin)
%
% PARAMETERS
%       coords - coordinates of the the end points of the edges
%       edges - end points of the edges
%
%       OPTIONAL
%
%       varargin - argument directly applied in plot.m
%               see help plot for all the options
%
% Copyright 2011-2012 University of Bergen, 2013 IRIS AS
%
% This file is licensed under the GNU General Public License v3.0.
numEdges = size(edges,1);
edgeCoords = [coords(edges(:,1),:),coords(edges(:,2),:)];
xedge = reshape([edgeCoords(:,[1,3]),nan(numEdges,1)]',[],1);
yedge = reshape([edgeCoords(:,[2,4]),nan(numEdges,1)]',[],1);
plot(xedge,yedge,varargin{:});