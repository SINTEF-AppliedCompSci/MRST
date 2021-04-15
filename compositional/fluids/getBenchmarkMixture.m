function [fluid, info] = getBenchmarkMixture(name, varargin)
% Grab bag of different fluids for use in examples etc
    descr = name;
    switch(lower(name))
        case 'lumped_1'
            % From SPE 79691 ex 5
            descr = 'Lumped mixture, SPE 79691, example 5, Mallison et al.';
            T_c = [189.515, 304.200, 387.607, 597.497, 698.515, 875.00];
            P_c = [45.2012, 72.90, 40.4196, 33.0150, 17.4525, 11.5372]*atm;
            mw = [16.1594, 44.01, 45.5725, 117.740, 248.827, 481.520]/1000;
            acc = [0.00854, 0.228, 0.16733, 0.38609, 0.80784, 1.23141];
            Z_c = [0.28981, 0.27055, 0.27588, 0.25668, 0.21967, 0.18250];
            
            V_c = Z_c.*8.314.*T_c./P_c;
            
            names = {'N2/CH4', 'CO2', 'C2-5', 'C6-13', 'C14-24', 'C25-80'};
            fluid = CompositionalMixture(names, T_c, P_c, V_c, acc, mw);
            
            bic = [0.11883, 0, 0, 0,0, 0;...
                   0.00070981, 0.15,0, 0, 0, 0; ...
                   0.00077754, 0.15,0, 0, 0, 0; ...
                   0.01, 0.15, 0, 0,0, 0; ...
                   0.011, 0.15, 0, 0,0, 0; ...
                   0.011, 0.15, 0, 0,0, 0];
               
            bic = bic + tril(bic, -1)';
            
            fluid = fluid.setBinaryInteraction(bic);
            info = makeInfo('injection', [0, 1, 0, 0, 0, 0], ...
                            'initial',   [0.463, 0.01640, 0.20520, 0.19108, 0.08113, 0.04319], ...
                            'pressure', 225*atm, ...
                            'T', 387.45);
        case 'spe5'
            descr = 'SPE5 benchmark';
            % Fifth SPE benchmark
            T_c = [343, 665.7, 913.4, 1111.8, 1270.0, 1380.0]*Rankine;
            P_c = [667.8, 616.3, 436.9, 304.0, 200.0, 162.0]*psia;
            mw = [16.040, 44.100, 86.180, 142.290, 206.000, 282.000]/1000;
            acc = [0.0130, 0.1524, 0.3007, 0.4885, 0.6500, 0.8500];
            Z_c = [0.290, 0.277, 0.264, 0.257, 0.245, 0.235];
            
            V_c = Z_c.*8.314.*T_c./P_c;
            
            names = {'C1', 'C3', 'C6', 'C10', 'C15', 'C20'};
            fluid = CompositionalMixture(names, T_c, P_c, V_c, acc, mw);
            
            ncomp = numel(names);
            bic = zeros(ncomp, ncomp);
            bic(5, 1) = 0.05;
            bic(6, 1) = 0.05;
            bic(5, 2) = 0.005;
            bic(6, 2) = 0.005;
            
            bic = bic + tril(bic, -1)';
            
            fluid = fluid.setBinaryInteraction(bic);
            info = makeInfo('injection', [0.77, 0.20, 0.03, 0, 0, 0], ...
                            'initial',   [0.5, 0.03, 0.07, 0.20, 0.15, 0.05], ...
                            'pressure', 4000*psia, ...
                            'T', 344.26);

        case 'fourmix'
            names = {'Methane', 'Nitrogen', 'n-Pentane', 'n-Decane'};
            
            info = makeInfo('injection', [0, 1, 0, 0], ...
                            'initial',   [0.1, 0.0, 0.45, 0.45], ...
                            'pressure', 75*barsa, ...
                            'T', 273.15 + 50);
                        
            fluid = TableCompositionalMixture(names);
        case '20components'
            names = {'Methane', 'Nitrogen', 'n-Pentane', 'n-Decane', 'CycloHexane', ...
                     'CarbonMonoxide', '1-Butene', 'Ammonia', 'HydrogenSulfide', 'IsoButane', ...
                     'n-Dodecane', 'n-Heptane', 'n-Hexane', 'n-Nonane', 'n-Octane', ...
                     'p-Xylene', 'n-Propane', 'n-Undecane', 'Oxygen', 'Water'};
            info = makeInfo('injection', (1:20)/20, ...
                            'initial',   ones(1, 20)/20, ...
                            'pressure', 75*barsa, ...
                            'T', 273.15 + 50);
                        
            fluid = TableCompositionalMixture(names);

        case {'liquid_initial', 'vapor_initial'}
            names = {'Oxygen', 'Water'};
            if strcmpi(name, 'liquid_initial')
                inj =  [1, 0];
                init = [0, 1];
            else
                inj =  [0, 1];
                init = [1, 0];
            end            
            info = makeInfo('injection', inj, ...
                            'initial',   init, ...
                            'pressure', 75*barsa, ...
                            'T', 273.15 + 30);
                        
            fluid = TableCompositionalMixture(names);
            fluid.names{2} = 'H2O';
            fluid.names{1} = 'O2';
        case 'watertracer'
            fluid = TableCompositionalMixture({'water', 'water'});
            info = makeInfo('injection', [1, 0], ...
                            'initial',   [0, 1], ...
                            'pressure', 75*barsa, ...
                            'T', 273.15 + 30);
            fluid.names = {'Tracer1', 'Tracer2'};
        case 'simple'
            names = {'Methane', 'CarbonDioxide', 'n-Decane'};
            
            info = makeInfo('injection', [0.1, 0.9, 0], ...
                            'initial',   [0.3, 0.1, 0.6], ...
                            'pressure', 75*barsa, ...
                            'T', 273.15 + 150);
                        
            order = 1:3;
            info.injection = info.injection(order);
            info.initial = info.initial(order);
            fluid = TableCompositionalMixture(names(order));
        case 'verysimple'
            names = {'CarbonDioxide', 'n-Decane'};
            fluid = TableCompositionalMixture(names);
            info = makeInfo('injection', [1, 0], ...
                            'initial',   [0.1, 0.9], ...
                            'pressure', 75*barsa, ...
                            'T', 273.15 + 150);
        case 'onlydecane'
            names = {'n-Decane'};
            fluid = TableCompositionalMixture(names);
            info = makeInfo('injection', [1], ...
                            'initial',   [1], ...
                            'pressure', 75*barsa, ...
                            'T', 273.15 + 150);
        case 'watermethane'
            names = {'Methane', 'Water'};
            fluid = TableCompositionalMixture(names);
            info = makeInfo('injection', [1, 0], ...
                            'initial',   [0.01, 0.99], ...
                            'pressure', 75*barsa, ...
                            'T', 273.15 + 150);
        otherwise
            error('Unknown case')
    end
    fluid.name = descr;
end

function info = makeInfo(varargin)
    info = struct('injection', [], 'initial', [],...
                  'pressure',  [], 'temp',    [], 'T', []);
    info = merge_options(info, varargin{:});
    info.temp = info.T;
end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
