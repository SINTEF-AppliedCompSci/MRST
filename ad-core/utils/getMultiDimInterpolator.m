function [fn, F, epsilons] = getMultiDimInterpolator(x, Fx, extrap)
% Get a multidimensional interpolator (with support for ADI varibles)
%
% SYNOPSIS:
%       fn = getMultiDimInterpolator(x, Y);
%
% PARAMETERS:
%   x       - Cell array containing the vectors for the points of the
%             function. Must have equal length to the number of dimensions
%             in Y.
%
%   Y       - A matrix of dimension equal to x, representing a structured
%             dataset for multdimensional interpolation
%
% RETURNS:
%   fn      - A function for interpolation with variable arguments equal to
%             the length of x.
%
%   F       - griddedInterpolant class instace
%
%   epsilons- Epsilon values for each input argument representing one half
%             of the minimum distance between elements. Useful for
%             computing numerical derivatives of the interpolant since
%             Matlab does not expose the slopes.
%
% SEE ALSO:
%   `interpTable`

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

    if nargin == 2
        extrap = 'linear';
    end
    
    sz = size(Fx);
    assert(iscell(x));
    if numel(x) > 1
        assert(numel(sz) == numel(x));
        assert(all(cellfun(@numel, x) == sz));
    end
    
    nvar = numel(x);
    dFdx = cell(nvar, 1);
    
    for i = 1:nvar
        if diff(x{i}(1:2)) < 0
            x{i} = x{i}(end:-1:1);
            Fx = flip(Fx, i);
        end
    end
    
    F = griddedInterpolant(x, Fx, 'linear', extrap);
    epsilons = cellfun(@(x) min(abs(diff(x))), x)./2;
    assert(min(epsilons) > 1e-18, 'Difference between points must be larger than 2e-18');
    
    fn = @(varargin) interpTableND(F, epsilons, varargin{:});
end

function Ye = interpTableND(Y, epsilon, varargin)
    Yq = varargin;
    nvar = numel(epsilon);
    assert(numel(Yq) == nvar);
    
    isad = cellfun(@(x) isa(x, 'ADI'), Yq);
    isnewad = cellfun(@(x) isa(x, 'GenericAD'), Yq);
    if any(isad)
        % Get a sample variable
        ad = Yq(isad);
        ad = ad{1};
        Y.ExtrapolationMethod = 'linear';
        
        % Create double as input for the function evaluation
        Yqd = cellfun(@value, Yq, 'UniformOutput', false);
        Ye = Y(Yqd{:});
        Yev = Ye;
        % Cast to ADI
        if any(isnewad)
            Ye = double2GenericAD(Ye, ad);
        else
            Ye = double2ADI(Ye, ad);
        end
        for i = 1:nvar
            if isad(i)
                % If it is AD, compute derivative using chain rule
                permuted = Yqd;
                e = epsilon(i);
                permuted{i} = permuted{i} + e;
                Yd = Y(permuted{:});
                dydx = (Yd - Yev)./e;
                
                ix = (1:numel(dydx))';
                d = sparse(ix, ix, dydx);
                for j = 1:numel(Ye.jac)
                    Ye.jac{j} = Ye.jac{j} + d*Yq{i}.jac{j};
                end
            end
        end
    else
        % Just use interpolator directly
        Ye = Y(Yq{:});
    end
end
