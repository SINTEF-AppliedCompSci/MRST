function [] = startup()


d = fileparts(mfilename('fullpath'));

vor2D = strcat(d,'/voronoi2D');
dataset = strcat(d,'/voronoi2D/examples/datasets');
dist  = strcat(d,'/distmesh');
util  = strcat(d,'/util');
addpath(vor2D)
addpath(dist)
addpath(util)
addpath(dataset)
end