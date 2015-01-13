function [pvMult, transMult, mobMult, pvMult0, transMult0, mobMult0] = getMultipliers(fluid, p, p0)

    [pvMult, transMult, mobMult, pvMult0, transMult0, mobMult0] = deal(1);
    % Pressure dependent pore volume multiplier
    if isfield(fluid, 'pvMultR')
        pvMult =  fluid.pvMultR(p);
        if nargout > 3
            pvMult0 = fluid.pvMultR(p0);
        end
    end
    % Pressure dependent mobility multiplier 
    if isfield(fluid, 'tranMultR')
        mobMult = fluid.tranMultR(p);
        if nargout > 4
            mobMult0 = fluid.tranMultR(p0);
        end
    end
    % Pressure dependent transmissibility multiplier
    if isfield(fluid, 'transMult')
       transMult = fluid.transMult(p); 
       if nargout > 5
           transMult0 = fluid.transMult(p0); 
       end
    end
end