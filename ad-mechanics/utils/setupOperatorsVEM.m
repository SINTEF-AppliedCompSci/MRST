function operators = setupOperatorsVEM(G, C, el_bc, load, alpha_scaling, S)
%
%
% SYNOPSIS:
%   function operators = setupOperatorsVEM(G, C, el_bc, load, alpha_scaling, S)
%
% DESCRIPTION: Assemble the discretization operators for the mechanical
% system, using VEM.
%
% PARAMETERS:
%   G             - Grid structure
%   C             - Stiffness tensor
%   el_bc         - Structure describing the boundary conditions (see VEM_linElast)
%   load          - Structure giving the volumetric forces (see VEM_linElast)
%   alpha_scaling - Coefficient of the stabilisation term (default 1, see VEM_linElast)
%   S             - Stabilization matrix to use (used only in very special
%                   cases, experimental feature, see VEM_linElast)
%
% RETURNS:
%   operators - Structure with the VEM discretization operators.
%
% EXAMPLE:
%
% SEE ALSO: VEM_linElast
%


    [~, extra] = VEM_linElast(G                    , ...
                              C                   , ...
                              el_bc               , ...
                              load                , ...
                              'alpha_scaling', alpha_scaling , ...
                              'S', S                         , ...
                              'linsolve', @(A, rhs) 0 * rhs);


    operators               = extra.disc;
    operators.global_strain = extra.WC' * extra.assemb';
    operators.strain        = operators.global_strain(:, ~operators.isdirdofs);
    operators.global_stress = extra.D * operators.global_strain;
    operators.stress        = extra.D * operators.strain;

    griddim = G.griddim;

    if (griddim == 3)
        vtrace = [1; 1; 1; 0; 0; 0];
    else
        vtrace = [1; 1;  0];
    end

    operators.trace = @(tensor) (tensor*[1; 1; 1; 0; 0; 0]);

end