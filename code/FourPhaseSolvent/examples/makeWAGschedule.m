function [schedule, W_G, W_W] = makeWAGschedule(model, injectors, producers, nCycles, varargin)

opt = struct('wag_end', 0.5, ...
             'gas_end'  , 0.5, ...
             'T', 1*year, ...
             'nStep', 100, ...
             'gRate', [], ...
             'wRate', [], ...
             'W',     []);
         
opt = merge_options(opt, varargin{:});

G = model.G;
rock = model.rock;
pv = poreVolume(G, rock);

gRate = opt.gRate;
wRate = opt.wRate;
T = opt.T;
if isempty(opt.gRate)
    gRate = sum(pv)/T;
end
if isempty(opt.gRate)
    wRate = sum(pv)/T;
end
wag_end = opt.wag_end;
T = opt.T;
Tc = T*wag_end/nCycles;
dT = T/opt.nStep;
gas_end = opt.gas_end;

dT_G = rampupTimesteps(gas_end*Tc, dT, 0);
dT_W  = rampupTimesteps((1-gas_end)*Tc, dT, 0);

dT_WAG = repmat([dT_G; dT_W], nCycles, 1);
dT_Water = rampupTimesteps(T - sum(dT_WAG), dT, 0);
dT = [dT_WAG; dT_Water];


[W_G, W_W] = deal([]);

W = opt.W;
if isempty(W)

    for inj = 1:numel(injectors)

        W_G = addWell(W_G, G, rock, injectors{inj}, ...
                      'comp_i', [0, 0, 0, 1], ...
                      'type'  , 'rate', ...
                      'val'   , gRate);

        W_W = addWell(W_W, G, rock, injectors{inj}, ...
                      'comp_i', [1, 0, 0, 0], ...
                      'type'  , 'rate', ...
                      'val'   , wRate);
    end

    for prod = 1:numel(producers)

        W_G = addWell(W_G, G, rock, producers{prod}, ...
                      'comp_i', [1, 0, 0, 0], ...
                      'type'  , 'bhp', ...
                      'val'   , 0);

        W_W = addWell(W_W, G, rock, producers{prod}, ...
                      'comp_i', [1, 0, 0, 0], ...
                      'type'  , 'bhp', ...
                      'val'   , 0);
    end

else
    
    [W_G, W_W] = deal(W);
    
    for i = 1:numel(W)
        W_G(i).compi = [0, 0, 0, 1];
        W_W(i).compi = [1, 0, 0, 0];
    end
    
end
    
control(1).W = W_G;
control(2).W = W_W;

step.control = [repmat([1*ones(numel(dT_G),1); 2*ones(numel(dT_W),1)],nCycles,1); ...
                2*ones(numel(dT_Water),1)];
step.val = dT;

schedule.control = control;
schedule.step = step;

end

