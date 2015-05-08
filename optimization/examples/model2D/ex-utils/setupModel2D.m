% setupModel:  construct model for example analyseModel2D

% set up grid and rock-properties for simple 2D-model
nxyz = [40, 40, 1];
Dxyz = [400, 400, 10];

G = cartGrid(nxyz, Dxyz);
G = computeGeometry(G);

rock = SPE10_rock(1:40, 101:140, 4);
rock.perm = convertFrom(rock.perm, milli*darcy);

% fluid
pRef = 200*barsa;

fluid = initSimpleADIFluid('mu',    [.3, 5, 0]*centi*poise, ...
                           'rho',   [1000, 700, 0]*kilogram/meter^3, ...
                           'n',     [2, 2, 0]);
c = 1e-5/barsa;
p_ref = 200*barsa;
fluid.bO = @(p) exp((p - p_ref)*c);

W = [];
% Injectors (lower-left and upper-right)
ci(1) = 1;
ci(2) = G.cells.num;
for k  = 1:2
    W = addWell(W, G, rock, ci(k), 'Type' , 'rate', ...
                                   'Val'  , 300*meter^3/day, ...
                                   'Name' , sprintf('I%d', k), ...
                                   'comp_i', [1 0], ...
                                   'Sign' , 1);
end
% Producers (upper-left and -right)
cp(1) = G.cartDims(1);
cp(2) = 1 + (G.cartDims(2)-1)*G.cartDims(1);
for k  = 1:2
    W = addWell(W, G, rock, cp(k), 'Type', 'bhp', ...
                                   'Val' , 150*barsa, ...
                                   'Name', sprintf('P%d', k), ...
                                   'comp_i', [0 1], ...
                                   'Sign', -1);
end


% Set up 4 control-steps each 150 days
ts = { [1 1 3 5 5 10 10 10 15 15 15 15 15 15 15]'*day, ...
                   repmat(150/10, 10, 1)*day, ...
                   repmat(150/6, 6, 1)*day, ...
                   repmat(150/6, 6, 1)*day};
       
numCnt = numel(ts);
[schedule.control(1:numCnt).W] = deal(W);
schedule.step.control = rldecode((1:4)', cellfun(@numel, ts));
schedule.step.val     = vertcat(ts{:});

gravity on


