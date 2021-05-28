function varargout = initVariablesAD_diagonalRowMajor(varargin)
% Diagonal AD initializer

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

    assert(nargout == nargin || nargout == nargin - 1);
    
    n_in = nargin;
    if n_in > nargout
        opts = varargin{nargout+1};
        varargin = varargin(1:end-1);
        n_in = n_in - 1;
        types = opts.types;
        useMex = opts.useMex;
    else
        types = [];
        useMex = false;
    end
    numVals = cellfun(@numel, varargin)';

    if isempty(types)
        types = ones(n_in, 1);
        for i = 2:n_in
            if numVals(i) ~= numVals(i-1)
                types(i) = types(i-1) + 1;
            else
                types(i) = types(i-1);
            end
        end
    end
    offset = zeros(max(types)+1, 1);
    ctr = 2;
    offset(1) = 1;
    for i = 2:numel(types)
        if types(i) ~= types(i-1)
            offset(ctr, 1) = i;
            ctr = ctr+1;
        end
    end
    offset(end) = numel(types) + 1;
    
    ntypes = max(types);

    n = nargout;
    varargout = cell([1, n]);
    for i = 1:n
        zerojac = cell(1, ntypes);
        for j = 1:ntypes
            type = j;
            sub = types == type;
            nv = numVals(sub);
            nval = nv(1);
            dim = [nval, nnz(sub)];
            d = zeros(0, numVals(i));
            zerojac{j} = DiagonalJacobian(d, dim, zeros(numVals(i), 1), useMex, true);
        end
        varargout{i} = GenericAD(varargin{i}, zerojac, numVals, offset, useMex);
    end
    for type = 1:ntypes
        sub = find(types == type);
        nv = numVals(sub);
        nval = nv(1);
        nsub = numel(sub);
        
        assert(all(nv == nval));
        djac = DiagonalJacobian(zeros(nsub, nval), [nval, nsub], [], useMex, true);
        for i = 1:nsub
            s = sub(i);
            J = djac;
            J.diagonal(i, :) = 1;
            varargout{s}.jac{type} = J;
        end
    end
end
