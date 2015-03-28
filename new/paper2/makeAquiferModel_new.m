function aquifer = makeAquiferModel_new(varargin)
% Make a model of a 1D antiform aquifer with small-scale caprock undulations
%
% SYNOPSIS
%   aquifer = makeAquiferModel()
%   aquifer = makeAquiferModel('pn1','pv1')
%
% PARAMETERS:
%   'pn'/pv - List of property names/property values that match those given
%   in the description of the model below
%
% OUTPUT:
%   aquifer - struct with fields: G, Gt, rock, rock2D, W
%
% DESCRIPTION:
%   The function will make a model of a 1D antiform aquifer with
%   small-scale caprock undulations. The model is realized on a 1D grid
%   described by the following parameters (with default values)
%      Lx  = 30 km;
%      Ly  = 10 km;
%      H   = 50 m;
%      nx  = 1000;
%   The form of the caprock is given by the following equation
%      z = D - L1 * sin(x/L1)*tan(phi) + A sin(2*pi*x/L2)
%   Default values:
%      D   =  2300 m;
%      A   =  2 m;
%      L1  =  20 km;
%      L2  =  300 m;
%      phi =  0.3;
%   In addition, we set a porosity 0.2 and a 100 mD permeability.
%   CO2 is injected in the deepest cell and the pressure is held constant
%   at 300 bar at the shallowest cell (using an artificial producing well).
%
% SEE ALSO:
%   makeFluidModel

   opt = struct(...
       'Lx', 30e3, 'Ly', 10e3, 'H', 50, 'nx', 1000,...
       'D', 2300, 'A', 2, 'L1', 20e3, 'L2',.3e3, 'phi', 0.03); 

   opt = merge_options(opt, varargin{:});

   % Create the grid
   G = cartGrid([opt.nx, 1, 1], [opt.Lx, opt.Ly, opt.H]); 
   x = G.nodes.coords(:, 1); 
   G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + opt.D ...
       - opt.L1 * sin(x / opt.L1) * tan(opt.phi) ...
       + opt.A * sin(2 * pi * x / opt.L2); 

   G = computeGeometry(G); 
   Gt = topSurfaceGrid(G); 

   % Create petrophysical model
   rock = struct('perm', 100 * milli * darcy * ones(G.cells.num, 1),...
                 'poro', 0.2 * ones(G.cells.num, 1)); 
   rock2D = averageRock(rock, Gt);

   % Set well position
   W = createSampleWell_new([], G, rock, floor(0.1 * opt.nx),...
                            'Type', 'rate', 'Val', 1e6 / year,...
                            'Radius', 0.125, 'Name', 'I', 'Comp_i', [0 1]); 
   W = createSampleWell_new(W, G, rock, G.cartDims(1),...
                            'Type', 'bhp', 'Val', 300 * barsa,...
                            'Radius', 0.125, 'Name', 'P', 'Sign', -1, 'Comp_i', [1 0]); 

   aquifer.G = G; 
   aquifer.Gt = Gt; 
   aquifer.rock = rock; 
   aquifer.rock2D = rock2D; 
   aquifer.W = W; 
end

