function states3D = VEstates23D(statesVE, Gt, fluid, varargin)
% for the moment, only converts pressure and saturation

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
            [s3D.s, s3D.sMax] = height2Sat(h, hmax, Gt, fluid.res_water, fluid.res_gas);
        else
            % capillary fringe model
            [s3D.s, s3D.sMax] = height2Sat(h, hmax, Gt, fluid.res_water, fluid.res_gas, ...
                                           'invPc3D', fluid.invPc3D, ...
                                           'rhoW', fluid.rhoW(state.pressure), ...
                                           'rhoG', fluid.rhoG(state.pressure));
        end
        
        s3D.s = [1-s3D.s, s3D.s];
        
        states3D{i} = s3D;
    end
end

% ----------------------------------------------------------------------------
function res = is_sharp_interface(fluid)
    res = strcmpi(fluid.relperm_model, 'sharp_interface_simple') || ...
          strcmpi(fluid.relperm_model, 'sharp_interface_integrated');
end
