function states3D = VEstates23D(statesVE, Gt, fluid, varargin)
% for the moment, only converts pressure and saturation

%{
Copyright 2009-2026 SINTEF Digital, Mathematics & Cybernetics.

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

    opt.poro3D = [];
    opt.tol = 1e-5;
    opt = merge_options(opt, varargin{:});
    
    G = Gt.parent;
    states3D = cell(numel(statesVE), 1);
    
    for i = 1:numel(statesVE)
        
        state = statesVE{i};
        
        s3D.pressure = pVE_to_p3D(Gt, state.pressure, fluid);
        
        if strcmpi(fluid.relperm_model, 'sharp_interface_simple')
            % vertically homogeneous sharp-interface model
            [h, hmax] = upscaledSat2height(state.s(:,2), state.sGmax, Gt, ...
                                           'resSat', [fluid.res_water, fluid.res_gas]);
            
        elseif strcmpi(fluid.relperm_model, 'sharp_interface_integrated')
            % vertically heterogeneous sharp-interface model
            if isempty(opt.poro3D)
                error(['To convert upscaled saturation to height for a vertically ' ...
                       'heterogeneous model, the fine-scale porosity field ' ...
                       'must be provided  in the `poro3D` optional ' ...
                       'argument']);
            end
            
            [h, hmax] = upscaledSat2height(state.s(:,2), state.sGmax, Gt, ...
                                           'resSat', [fluid.res_water, fluid.res_gas], ...
                                           'poro', opt.poro3D, ...
                                           'tol', opt.tol);
        else 
            % capillary fringe model
            [h, hmax] = upscaledSat2height(state.s(:,2), state.sGmax, Gt, ...
                                           'pcWG', fluid.pcWG, ...
                                           'rhoW', fluid.rhoW, ...
                                           'rhoG', fluid.rhoG, ...
                                           'p', state.pressure);
        end
        
        if is_sharp_interface(fluid)
            [s3D.s, s3D.sMax] = height2finescaleSat(h, hmax, Gt, fluid.res_water, fluid.res_gas);
        else
            % capillary fringe model
            [s3D.s, s3D.sMax] = height2finescaleSat(h, hmax, Gt, fluid.res_water, fluid.res_gas, ...
                                                    'invPc3D', fluid.invPc3D, ...
                                                    'hyst_model', fluid.hyst_model, ...
                                                    'rhoW', fluid.rhoW(state.pressure), ...
                                                    'rhoG', fluid.rhoG(state.pressure));
        end
        
        s3D.s = [1-s3D.s, s3D.s];

        if isfield(state, 'rs')
            % dissolution of CO2 in water is included.  Convert upscaled dissolved CO2 saturation
            % to fine-scale
            rs_max = fluid.dis_max;
            fraction = state.rs ./ rs_max; % fraction of maximum dissolved CO2 saturation
            poro3D = opt.poro3D;
            if isempty(poro3D)
                poro3D = ones(G.cells.num, 1); % asssumed vertically uniform (value will be
                                               % normalized away, so won't affect results)
            end

            % compute total amount of brine in each column
            bvol3D = poro3D .* s3D.s(:, 1);
            bvol2D = integrateVertically(bvol3D, Gt.cells.H, Gt);
            bvol2D_saturated = fraction .* bvol2D; % volume of brine that is  saturated with
                                                   % dissolved CO2 in each column

            % compute how far the dissolved CO2 has penetrated downinto the column.
            % Find the root of the equation using a simple Newton-Raphson method
            h_guess = fraction .* Gt.cells.H; % initial guess for penetration depth
            tol = 1e-10; % should be enough to determine interface depth to satisfactory accuracy
            error = inf;
            max_iter = 1000; % usually we should never need this
            cur_iter = 0;
            while cur_iter < max_iter
                cur_iter = cur_iter + 1;
                [f, df] = integrateVertically(bvol3D, h_guess, Gt);
                f = f - bvol2D_saturated; % value of the equation at current guess
                error = abs(f); % update error
                if error < tol
                    break; % solution found
                end
                h_guess = h_guess - f ./ df; % update guess using Newton-Raphson step
            end
            if cur_iter == max_iter
                warning('Maximum number of iterations reached while computing dissolved CO2 penetration depth. Results may be inaccurate.');
            end

            % use the computed penetration depth to determine the fine-scale
            % dissolved CO2 saturation.  We can hijack height2finescaleSat
            % function for this, by treating the dissolved CO2 saturation as a
            % sharp-interface model with zero residual saturation and zero
            % maximum saturation.
            s3D.rs = height2finescaleSat(h_guess, h_guess, Gt, 0, 0) .* rs_max; 
            
        end

        
        states3D{i} = s3D;
    end
end

% ----------------------------------------------------------------------------
function res = is_sharp_interface(fluid)
    res = strcmpi(fluid.relperm_model, 'sharp_interface_simple') || ...
          strcmpi(fluid.relperm_model, 'sharp_interface_integrated');
end
