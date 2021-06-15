function v = discreteDivergence(acc, flux, options)
% Discrete divergence for the GenericAD library

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

    hasAcc = ~isempty(acc);
    accAD = isa(acc, 'GenericAD');
    fluxAD = isa(flux, 'GenericAD');
    if accAD || fluxAD
        % We have either as AD. Pick one of them to be our AD instance to
        % avoid re-instantiating a class.
        if accAD
            v = acc;
        else
            v = flux;
        end
        v.val = accumulate(value(acc), value(flux), options);
        if hasAcc && accAD
            % If the accumulation term is the only AD variable, we do not
            % need to do anything. If the fluxes are AD, we need to perform
            % acc + div.
            if fluxAD
                % Both present, both are AD
                for i = 1:numel(v.jac)
                    v.jac{i} = discreteDivergenceDiagonalJac(acc.jac{i}, flux.jac{i}, options);
                end
            end
        else
            % Just the flux - accumulation terms are doubles
            for i = 1:numel(v.jac)
                v.jac{i} = discreteDivergenceDiagonalJac([], flux.jac{i}, options);
            end
        end
    else
        assert(isnumeric(flux), 'Expected numeric vector, but got ''%s''\n', class(flux))
        v = accumulate(acc, flux, options);
    end
end
function v = accumulate(acc, v, opt)
    if opt.useMex
        v = mexDiscreteDivergenceVal(acc, v, opt.N, opt.nc, opt.mex.facePos, opt.mex.faces);
    else
        v = accumarray(opt.N(:, 1), v, [opt.nc, 1]) - accumarray(opt.N(:, 2), v, [opt.nc, 1]);
        if ~isempty(acc)
            v = v + acc;
        end
    end
end

