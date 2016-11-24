%% Grid-orientation and anisotropy effects
% This script contains two examples that originate from the first MRST
% paper: Lie et al, "Open source MATLAB implementation of consistent
% discretisations on complex grids". Comput. Geosci., 16(2):297-322, 2012.
% DOI: 10.1007/s10596-011-9244-4
addpath(fullfile(fileparts(mfilename('fullpath')), 'src'))
mrstModule add incomp mimetic mpfa

%% First example
% This example corresponds to Figure 7 in the paper, which illustrates
% grid-orientation effects for the TPFA scheme and reproduction of linear
% flow for the mimetic and the MPFA-O method on a perturbed grid for a
% homogeneous, diagonal permeability tensor with entries Kx=1 and Ky=1000.

% Grid and permeability
G = cartGrid([51, 51]);
G.nodes.coords = twister(G.nodes.coords, 0.03);

% Permeability
K = diag([1, 1000]);

% Seed for streamline tracing
seed = (ceil(G.cartDims(1)/2):4*G.cartDims(1):prod(G.cartDims))';

% Run example
showMonotonicityExample(G, K, seed, true);
colormap(.75*jet(32) + .25*ones(32,3));

%% Second example
% This example corresponds to Figure 8 in the paper,  which illustrates
% montonicity effects for the TPFA, mimetic, and MPFA schemes. Same setup
% as in the first example, but now with the anisotropy ration of 1:1000
% making 30 degree angle with the grid directions.

% Grid and permeability
G = cartGrid([21, 21]);
G.nodes.coords = twister(G.nodes.coords, 0.03);

% Permeability
t  =  30*pi/180;
U  = [ cos(t), sin(t); -sin(t), cos(t)];
Kd = diag([1,1000]);
K  =  U'*Kd*U;

% Start points for streamline tracing
seed = (ceil(G.cartDims(1)/2):G.cartDims(1):prod(G.cartDims))';

% Run example
showMonotonicityExample(G, K, seed, true);