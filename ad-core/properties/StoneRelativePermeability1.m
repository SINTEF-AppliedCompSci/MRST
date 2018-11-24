classdef StoneRelativePermeability1 < GridProperty
    properties
        relpermScaling = false;
        relpermPoints = 2;
        regions_imbibition = [];
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
        if isfield(f, 'krPts')
            swcon = reshape(f.krPts.w(prop.regions, 1), [], 1);
        end
        swcon = min(swcon, double(sw)-1e-5);

        d  = (sg+sw-swcon);
        ww = (sw-swcon)./d;
        
        krW = evaluatePhaseRelativePermeability(prop, model, 'w', sw);
        krG = evaluatePhaseRelativePermeability(prop, model, 'g', sg);

        wg = 1-ww;
        if isfield(f, 'krO')
            krO = prop.evaluateFunctionOnGrid(f.krO, so);
        else
            krow = evaluatePhaseRelativePermeability(prop, model, 'ow', so);
            krog = evaluatePhaseRelativePermeability(prop, model, 'og', so);
            krO  = wg.*krog + ww.*krow;
        end
        kr = {krW, krO, krG};
    end
    
    function kr = evaluatePhaseRelativePermeability(prop, model, phase, s)
        fn = model.fluid.(['kr', upper(phase)]);
        f = model.fluid;
        if prop.relpermScaling
            [ss, kr_max] = prop.scaleSaturation(model.rock.krscale.drainage, phase, prop.regions, f, s);
            kr = kr_max.*prop.evaluateFunctionOnGrid(fn, ss);
        else
            kr = prop.evaluateFunctionOnGrid(fn, s);
        end
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
    
    function [s_scale, k_max] = scaleSaturation(prop, pts, phase, reg, f, s)
        k_max = pts.(phase)(:, 4);
        if prop.relpermPoints == 2 % 2-point
            [m, c, p] = getTwoPointScalers(pts, phase, reg, f);
            ix1 = s < p{1};
            ix2 = s >= p{2};
            ix  = ~or(ix1,ix2);
            s_scale = (ix.*m).*s + (ix.*c + ix2);
        elseif prop.relpermPoints == 3
            [m, c, p] = getThreePointScalers(pts, phase, reg, f);
            ix1 = and(s >= p{1}, s < p{2});
            ix2 = and(s >= p{2}, s < p{3});
            ix3 = s >= p{3};
            a = ix1.*m{1} + ix2.*m{2};
            s_scale = a.*s + ix1.*c{1} + ix2.*c{2} + ix3;
        else
            error('Unknown number of scaling points');
        end
    end
end
end


function [m, c, p] = getTwoPointScalers(scal, tab, ph, f)
    [get, CR, U] = getSatPointPicker(f, pts, reg);
    switch ph
        case {'w', 'g'}
            [su, SU] = get(ph, U);
        case {'ow', 'og'}
            [swl, SWL] = get('w', U);
            [sgl, SGL] = get('g', U);
            SU   = 1 - SWL - SGL;
            su = 1 - swl - sgl;
    end
    [scr, SCR] = get(ph, CR);
    p = cell(1, 2);
    p{1} = SCR;
    m    = (su-scr)./(SU-SCR);
    c    = scr - SCR.*m;
    p{2} = SU;
end

function [m, c, p] = getThreePointScalers(pts, ph, reg, f)
    [get, CR, U, L] = getSatPointPicker(f, pts, reg);
    switch ph
        case 'w'
            [sowcr, SOWCR] = get('ow', CR);
            [sgl, SGL] = get('g', L);
            
            SR = 1 - SOWCR - SGL;
            sr = 1 - sowcr - sgl;
            
            [su, SU] = get('w', U);  
        case 'g'
            [sogcr, SOGCR] = get('og', CR);
            [swl, SWL] = get('w', L);
            SR = 1 - SOGCR - SWL;
            sr = 1 - sogcr - swl;
            
            [su, SU] = get('g', U);  
        case 'ow'
            [swcr, SWCR] = get('w', CR);
            [sgl, SGL] = get('g', L);
            
            SR = 1 - SWCR - SGL;
            sr = 1 - swcr - sgl;

            [swl, SWL] = get('w', L);
            [sgl, SGL] = get('g', L);
            
            SU = 1 - SWL - SGL;
            su = 1 - swl - sgl;
        case 'og'
            [sgcr, SGCR] = get('g', CR);
            [swl, SWL] = get('w', L);
            [sgl, SGL] = get('g', L);
            SR   = 1 - SGCR - SWL;
            sr   = 1 - sgcr - swl;

            SU   = 1 - SGL - SWL;
            su = 1 - sgl - swl;
    end
    [scr, SCR] = get(ph, CR);

    m = cell(1, 2);
    p = cell(1, 3);
    c = cell(1, 2);

    p{1} = SCR;
    m{1} = (sr-scr)./(SR-SCR);
    c{1} = scr - SCR.*m{1};
    p{2} = SR;

    p{3} = SU;

    m{2} = (su-sr)./(SU-SR);
    c{2} = sr - SR.*m{2};
    ix = SU <= SR;
    if nnz(ix) > 0
        m{2}(ix) = 0;
        c{2}(ix) = 0;
    end
end

function [get, CR, U, L] = getSatPointPicker(f, pts, reg)
    L = 1;
    CR = 2;
    U = 3;
    
    tbl = @(phase, index) f.krPts.(phase)(reg, index);
    scal = @(phase, index) pts.(phase)(:, index);
    get = @(phase, index) getPair(phase, index, tbl, scal);
    
end

function [v1, v2] = getPair(phase, index, fn1, fn2)
    v1 = fn1(phase, index);
    v2 = fn2(phase, index);
end