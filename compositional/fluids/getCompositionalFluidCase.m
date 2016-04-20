function [fluid, info] = getCompositionalFluidCase(name, varargin)
% Grab bag of different fluids for use in examples etc
    switch(lower(name))
        case 'lumped_1'
            % From SPE 79691 ex 5
            T_c = [189.515, 304.200, 387.607, 597.497, 698.515, 875.00];
            P_c = [45.2012, 72.90, 40.4196, 33.0150, 17.4525, 11.5372]*atm;
            mw = [16.1594, 44.01, 45.5725, 117.740, 248.827, 481.520]/1000;
            acc = [0.00854, 0.228, 0.16733, 0.38609, 0.80784, 1.23141];
            Z_c = [0.28981, 0.27055, 0.27588, 0.25668, 0.21967, 0.18250];
            
            V_c = Z_c.*mw*8.314.*T_c./P_c;
            T_ref = 273.15 + 25;
            P_ref = 1*atm;
            
            names = {'N2/CH4', 'CO2', 'C2-5', 'C8-13', 'C14-24', 'C25-80'};
            fluid = CompositionalFluid(names, T_c, P_c, V_c, acc, mw, T_ref, P_ref);
            
            bic = [0.11883, 0, 0, 0,0, 0;...
                   0.00070981, 0.15,0, 0, 0, 0; ...
                   0.00077754, 0.15,0, 0, 0, 0; ...
                   0.01, 0.15, 0, 0,0, 0; ...
                   0.011, 0.15, 0, 0,0, 0; ...
                   0.011, 0.15, 0, 0,0, 0];
               
            bic = bic + tril(bic, -1)';
            
            fluid.setBinaryInteraction(bic);
            info = makeInfo('injection', [0, 1, 0, 0, 0, 0], ...
                            'initial',   [0.463, 0.01640, 0.20520, 0.19108, 0.08113, 0.04319], ...
                            'pressure', 225*atm, ...
                            'temp', 387.45);
        case 'simple'
            names = {'Methane', 'CarbonDioxide', 'n-Decane'};
            
            info = makeInfo('injection', [0.1, 0.9, 0], ...
                            'initial',   [0.3, 0.1, 0.6], ...
                            'pressure', 75*barsa, ...
                            'temp', 273.15 + 150);
                        
            order = 1:3;
            info.injection = info.injection(order);
            info.initial = info.initial(order);
            fluid = CoolPropsCompositionalFluid(names(order));
        case 'verysimple'
            names = {'CarbonDioxide', 'n-Decane'};
            fluid = CoolPropsCompositionalFluid(names);
            info = makeInfo('injection', [1, 0], ...
                            'initial',   [0.1, 0.9], ...
                            'pressure', 75*barsa, ...
                            'temp', 273.15 + 150);
        case 'onlydecane'
            names = {'n-Decane'};
            fluid = CoolPropsCompositionalFluid(names);
            info = makeInfo('injection', [1], ...
                            'initial',   [1], ...
                            'pressure', 75*barsa, ...
                            'temp', 273.15 + 150);

        otherwise
            error('Unknown case')
    end
end

function info = makeInfo(varargin)
    info = struct('injection', [], 'initial', [],...
                  'pressure',  [], 'temp',    []);
    info = merge_options(info, varargin{:});
end