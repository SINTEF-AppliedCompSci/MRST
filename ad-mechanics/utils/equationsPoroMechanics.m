function [eqs, names, types] = equationsPoroMechanics(x, model, fluidp)
%
%
% SYNOPSIS:
%   function [eqs, names, types] = equationsPoroMechanics(x, model, fluidp)
%
% DESCRIPTION: Assemble the residual equations for the mechanical system. The
% function takes fluid input, given by the pore pressure.
%
% PARAMETERS:
%   x      - Displacement
%   model  - Model class instance that is used
%   fluidp - Fluid pressure
%
% RETURNS:
%   eqs   - The residual values as ADI variables (that is with the Jacobian)
%           if the inputs were also ADI.
%   names - The name of each equations
%   types - The type of each equations
%
% EXAMPLE:
%
% SEE ALSO:
%

    
    G = model.G;
    s = model.operators;
    alpha = model.rock.alpha;

    eqs{1} = s.A * x - s.gradP * (alpha .* fluidp) - s.rhs;

    % normalization constant
    fac =  1 / (1e6 * mean(G.cells.volumes));
    eqs{1} = eqs{1} * fac;
    names = {'disp'};
    types = {'disp_dofs'};

end