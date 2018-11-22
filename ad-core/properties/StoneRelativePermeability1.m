classdef StoneRelativePermeability1 < GridProperty
    properties
    end
    
    methods
        function kr = evaluateOnGrid(prop, model, state)
            if model.water && model.gas && model.oil
                kr = prop.relPermWOG(model, state);
            elseif model.water && model.oil
                
            elseif model.water && model.gas
                
            elseif model.oil && model.gas
                
            end
            
        end

    % --------------------------------------------------------------------%
    function kr = relPermWOG(prop, model, state)
        [sw, so, sg] = model.getProps(state, 'sw', 'so', 'sg');
        f = model.fluid;
        swcon = 0;
        if isfield(f, 'sWcon')
            swcon = f.sWcon(prop.regions);
        end
        swcon = min(swcon, double(sw)-1e-5);

        d  = (sg+sw-swcon);
        ww = (sw-swcon)./d;
        
        krW = prop.evaluateFunctionOnGrid(f.krW, sw);
        krG = prop.evaluateFunctionOnGrid(f.krG, sg);

        wg = 1-ww;
        
        if isfield(f, 'krO')
            krO = prop.evaluateFunctionOnGrid(f.krO, so);
        else
            krow = prop.evaluateFunctionOnGrid(f.krOW, so);
            krog = prop.evaluateFunctionOnGrid(f.krOG, so);
            krO  = wg.*krog + ww.*krow;
        end
        kr = {krW, krO, krG};
    end

    % --------------------------------------------------------------------%
    function [krW, krO] = relPermWO(sw, so, f, varargin)
        % Two-phase water-oil relative permeability function
        %
        % SYNOPSIS:
        %   [krW, krO] = model.relPermWO(sw, so, f);
        %
        %
        % PARAMETERS:
        %   sw  - Water saturation
        %   so  - Oil saturation
        %   f   - Struct representing the field. Fields that are used:
        %
        %            - `krW`: Water relperm function of water saturation.
        %            - `krO`: Oil relperm function of oil saturation.
        %            - `krOW`: Oil-water relperm function of oil saturation.
        %              Only used if `krO` is not found.
        % RETURNS:
        %   krW - Water relative permeability.
        %   krO - Oil relative permeability.
        %
        % NOTE:
        %   This function should typically not be called directly as its
        %   interface is subject to change. Instead, use `evaluateRelPerm`.
        krW = f.krW(sw, varargin{:});
        if isfield(f, 'krO')
            krO = f.krO(so, varargin{:});
        else
            krO = f.krOW(so, varargin{:});
        end
    end

    % --------------------------------------------------------------------%
    function [krO, krG] = relPermOG(so, sg, f, varargin)
        % Two-phase oil-gas relative permeability function
        %
        % SYNOPSIS:
        %   [krO, krG] = model.relPermOG(so, sg, f);
        %
        %
        % PARAMETERS:
        %   sw  - Water saturation
        %   sg  - Gas saturation
        %   f   - Struct representing the field. Fields that are used:
        %
        %           - `krO`: Oil relperm function of oil saturation.
        %           - `krOG`: Oil-gas relperm function of gas saturation.
        %             This function is only used if `krO` is not found.
        %           - `krG`: Gas relperm function of gas saturation.
        %
        % RETURNS:
        %   krO - Oil relative permeability.
        %   krG - Gas relative permeability
        %
        % NOTE:
        %   This function should typically not be called directly as its
        %   interface is subject to change. Instead, use `evaluateRelPerm`.
        krG = f.krG(sg, varargin{:});
        if isfield(f, 'krO')
            krO = f.krO(so, varargin{:});
        else
            krO = f.krOG(so, varargin{:});
        end
    end

    % --------------------------------------------------------------------%
    function [krW, krG] = relPermWG(sw, sg, f, varargin)
        % Two-phase water-gas relative permeability function
        %
        % SYNOPSIS:
        %   [krW, krG] = model.relPermWG(sw, sg, f);
        %
        %
        % PARAMETERS:
        %   sw  - Water saturation
        %   sg  - Gas saturation
        %   f   - Struct representing the field. Fields that are used:
        %
        %           - `krW`: Water relperm function of water saturation.
        %           - `krG`: Gas relperm function of gas saturation.
        %
        % RETURNS:
        %   krW - Water relative permeability.
        %   krG - Gas relative permeability
        %
        % NOTE:
        %   This function should typically not be called directly as its
        %   interface is subject to change. Instead, use `evaluateRelPerm`.
        krG = f.krG(sg, varargin{:});
        krW = f.krW(sw, varargin{:});
    end
end
end