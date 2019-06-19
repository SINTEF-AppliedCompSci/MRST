classdef SurfactantRelativePermeability < BaseRelativePermeability

    properties
        zeroSurf
        fullSurf
    end
    
        
    methods
        function prop = SurfactantRelativePermeability(model, satreg, surfreg, varargin)
            prop@BaseRelativePermeability(model, varargin{:});
            prop = prop.dependsOn({'sat', 'c', 'capillaryNumber'}, 'state');
            prop.zeroSurf = StateFunction(model, satreg);
            prop.fullSurf = StateFunction(model, surfreg);
        end
        
        function kr = evaluateOnDomain(prop, model, state)

            fluid = model.fluid;
            [sat, cs, Nc] = model.getProps(state, 'sat', 'surfactant', 'capillaryNumber');

            isSft = (value(cs) > 0);
            m = 0*cs;
            if nnz(isSft) > 0
                logNc = log(Nc(isSft))/log(10);
                % We cap logNc (as done in Eclipse)
                logNc = min(max(-20, logNc), 20);
                m(isSft) = fluid.miscfact(logNc, 'cellInx', find(isSft));
            end
            
            [sW, sO] = model.getProps(state, 'sw', 'so');

            satreg  = model.rock.regions.saturation; 
            surfreg = model.rock.regions.surfactant;

            sWcon_noSft = fluid.krPts.w(satreg, 2);  % Residual water saturation   without surfactant
            sOres_noSft = fluid.krPts.w(satreg, 2);  % Residual oil saturation     without surfactant
            sWcon_Sft   = fluid.krPts.w(surfreg, 2); % Residual water saturation   with    surfactant
            sOres_Sft   = fluid.krPts.w(surfreg, 2); % Residual oil saturation     with    surfactant

            % Interpolated water/oil residual saturations
            sNcWcon = m.*sWcon_Sft + (1 - m).*sWcon_noSft;
            sNcOres = m.*sOres_Sft + (1 - m).*sOres_noSft;

            sNcWEff = (sW - sNcWcon)./(1 - sNcWcon - sNcOres);
            sNcOEff = (sO - sNcOres)./(1 - sNcWcon - sNcOres);

            % Rescaling of the saturation - without surfactant
            sW_noSft  = (1 - sWcon_noSft - sOres_noSft).*sNcWEff + sWcon_noSft;
            sO_noSft  = (1 - sWcon_noSft - sOres_noSft).*sNcOEff + sOres_noSft;
            % Compute rel perm - without surfactant
            krW_noSft = prop.zeroSurf.evaluateFunctionOnDomainWithArguments(fluid.krW, sW_noSft);
            if isfield(fluid, 'krO')
                krO_noSft = prop.zeroSurf.evaluateFunctionOnDomainWithArguments(fluid.krO, ...
                                                                 sO_noSft);
            else
                krO_noSft = prop.zeroSurf.evaluateFunctionOnDomainWithArguments(fluid.krOW, ...
                                                                 sO_noSft);
            end
            
            % Rescaling of the saturation - with surfactant
            sW_Sft  = (1 - sWcon_Sft - sOres_Sft).*sNcWEff + sWcon_Sft;
            sO_Sft  = (1 - sWcon_Sft - sOres_Sft).*sNcOEff + sOres_Sft;
            % Compute rel perm - with surfactant
            krW_Sft = prop.zeroSurf.evaluateFunctionOnDomainWithArguments(fluid.krW, sW_Sft);
            if isfield(fluid, 'krO')
                krO_Sft = prop.fullSurf.evaluateFunctionOnDomainWithArguments(fluid.krO, ...
                                                               sO_Sft);
            else
                krO_Sft = prop.fullSurf.evaluateFunctionOnDomainWithArguments(fluid.krOW, ...
                                                               sO_Sft);
            end

            % Interpolate relperm, with and without surfactant
            
            krW  = m.*krW_Sft + (1 - m).*krW_noSft;
            krO  = m.*krO_Sft + (1 - m).*krO_noSft;

            kr = {krW, krO};
        end
    end
end