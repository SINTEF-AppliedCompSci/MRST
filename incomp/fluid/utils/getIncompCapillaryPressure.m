function varargout = getIncompCapillaryPressure(state, fluid, varargin)
    varargout = cell(1, nargout);
    if isfield(fluid, 'properties')
        [varargout{:}] = getPropsLegacy(state, fluid);
    else
        [varargout{:}] = getPropsAD(state, fluid);
    end
end

function [pc, dpc] = getPropsLegacy(state, fluid)
    if isfield(fluid, 'pc')
        if nargout > 1
            [pc, dpc] = fluid.pc(state);
        else
            pc = fluid.pc(state);
        end
    else
        pc = [];
        dpc = [];
    end
end

function [pc, dpc] = getPropsAD(state, fluid)
    if isfield(fluid, 'pcOW')
        getDer = nargout > 1;
        s = state.s(:, 1);
        if getDer
            s = initVariablesAD_diagonal(s);
        end
        pc = fluid.pcOW(s);
        if getDer
            dpc = pc.jac{1}.diagonal; 
        end
    else
        pc = [];
        dpc = [];
    end
end