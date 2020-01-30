function [] = startup()


d = fileparts(mfilename('fullpath'));

vor2D = strcat(d,'/pebi2D');
vor3D = strcat(d,'/pebi3D');
dist  = strcat(d,'/distmesh');
util  = strcat(d,'/util');
addpath(vor2D)
addpath(vor3D)
addpath(dist)
addpath(util)
end