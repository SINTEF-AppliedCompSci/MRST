function [model, W, state0] = simpleModelForHMExample(varargin)
opt = struct('cellDim',        [40, 20, 15], ...
             'physDim',        [40, 20, 15], ...
             'fluid',                    [], ...
             'seed',                      0);

opt = merge_options(opt, varargin{:});         
%% Make model
% We make a small model that consists of two different facies with
% contrasting petrophysical properties. Two injectors and two producers are
% placed in the corner at different depths.
rng(opt.seed);
G  = computeGeometry(cartGrid(opt.cellDim, opt.physDim));
K1 = gaussianField(G.cartDims, [200 2000]);
p1 = K1(:)*1e-4 + .2;
K1 = K1(:)*milli*darcy;
K2 = gaussianField(G.cartDims, [10 500]);
p2 = K2(:)*1e-4 + .2;
K2 = K2(:)*milli*darcy;

rad1 = G.cells.centroids(:,1).^2 + .5*G.cells.centroids(:,2).^2 ...
   + (G.cells.centroids(:,3)-2).^2;
rad2 = .5*(G.cells.centroids(:,1)-40).^2 + 2*G.cells.centroids(:,2).^2 ...
   + 2*(G.cells.centroids(:,3)-2).^2;

ind = ((rad1>600) & (rad1<1500)) | ((rad2>700) & (rad2<1400));
rock.perm = K2(:);
rock.perm(ind) = K1(ind);
rock.perm = bsxfun(@times, rock.perm, [1 1 1]);
rock.poro = p2(:);
rock.poro(ind) = p1(ind);

%% wells
W = verticalWell([], G, rock, 4, 17, 4:10, 'Comp_i', [1 0], ...
   'Type', 'rate', 'Val', 100*stb/day, 'Name', 'I1');
W = verticalWell(W, G, rock, 2, 3, 1:6, 'Comp_i', [1 0], ...
    'Type', 'rate', 'Val', 100*stb/day, 'Name', 'I2');
W = verticalWell(W,  G, rock, 35, 3, 1:8, 'Comp_i', [0 1], ...
   'Type', 'bhp', 'Val', 95*barsa, 'Name', 'P1');
W = verticalWell(W,  G, rock, 35, 17, 8:15, 'Comp_i', [0 1], ...
   'Type', 'bhp', 'Val', 95*barsa, 'Name', 'P2');

%% Fluid model
fluid = opt.fluid;
if isempty(fluid)
    fluid = initSimpleADIFluid('phases', 'WO',... %Fluid phase
                               'mu' , [1, 10]*centi*poise     , ...%Viscosity
                               'rho', [1014, 859]*kilogram/meter^3, ...%Surface density [kg/m^3]
                               'n', [2 2]);
end
%% initial state and model
sw_conn = fluid.krPts.w(1);
state0  = initState(G, W, 100*barsa, [sw_conn, 1-sw_conn]);                    
% Initial state
% s0 = [0.2,0.8];
% p0 = 100*barsa;
% state0  = initState(G, W,p0,s0);

% Create model
model = GenericBlackOilModel(G, rock, fluid, 'gas', false);