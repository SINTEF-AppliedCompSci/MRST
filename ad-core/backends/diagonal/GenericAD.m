classdef GenericAD < ADI
    % GenericAD is the testbed for future updates to the ADI class. All
    % features herein are subject to rapid change.
    properties
        numVars
        offsets
        useMex = false;
    end
    methods
        function ad = GenericAD(val, jac, numVars, offset, useMex)
            ad@ADI();
            if nargin
                ad.val = val;
                if iscell(jac)
                    ad.jac = jac;
                else
                    ad.jac = {jac};
                end
                if nargin > 2
                    ad.numVars = numVars;
                    if nargin > 3
                        ad.offsets = offset;
                        if nargin > 4
                            ad.useMex = useMex;
                        end
                    end
                end
            end
        end
        
        function u = convertDouble(x, v)
            assert(isa(v, 'double'));
            nval  = numel(v);
            nj = numel(x.jac);
            jac = cell(1, nj);
            for i = 1:nj
                jx = x.jac{i};
                if issparse(jx)
                    jac{i} = sparse([], [], [], nval, size(jx, 2));
                else
                    jac{i} = jx.toZero(nval);
                end
            end
            u = x;
            u.val = v;
            u.jac = jac;
        end
        
        function x = incrementSubset(x, subs, v)
            if isa(x, 'GenericAD')
                x.val(subs) = x.val(subs) + value(v);
                if isa(v, 'GenericAD')
                    for i = 1:numel(x.jac)
                        x.jac{i} = incrementSubset(x.jac{i}, subs, v.jac{i});
                    end
                end
            else
                % V is AD
                x(subs) = x(subs) + v;
            end
        end
        
        function h = vertcat(varargin)
            isD = cellfun(@isnumeric, varargin);
            if any(isD)
                sampleAD = varargin(~isD);
                sampleAD = sampleAD{1};
                for i = 1:numel(isD)
                    if isD(i)
                        varargin{i} = double2GenericAD(varargin{i}, sampleAD);
                    end
                end
            end
            nv    = numel(varargin);
            nj    = numel(varargin{1}.jac);
            vals  = cell(1,nv);
            jacs  = cell(1,nj);
            sjacs = cell(1,nv);
            for k = 1:nv
                vals{k} = varargin{k}.val;
            end
            for k = 1:nj
                for k1 = 1:nv
                    sjacs{k1} = varargin{k1}.jac{k};
                end
                jacs{k} = GenericAD.vertcatJac(sjacs{:});
            end
            h = varargin{1};
            h.val = vertcat(vals{:});
            h.jac = jacs;
        end

        function h = combineEquations(varargin)
            isD = cellfun(@isnumeric, varargin);
            if any(isD)
                sampleAD = varargin(~isD);
                sampleAD = sampleAD{1};
                for i = 1:numel(isD)
                    if isD(i)
                        varargin{i} = double2GenericAD(varargin{i}, sampleAD);
                    end
                end
            end
            
            if false && nargin > 1
                nj = numel(varargin{1}.jac);
                nv = nargin;
                
                I = cell(nv, nj);
                J = cell(nv, nj);
                V = cell(nv, nj);
                N = 0;
                for i = 1:nv
                    M = 0;
                    for j = 1:nj
                        [I{i, j}, J{i, j}, V{i, j}, n, m] = getSparseArguments(varargin{i}.jac{j}, N, M);
                        M = M + m;
                    end
                    N = N + n;
                end
                h = varargin{1};
                v = cellfun(@(x) x.val, varargin, 'UniformOutput', false);
                h.val = vertcat(v{:});
                h.jac = {sparse(vertcat(I{:}), vertcat(J{:}), vertcat(V{:}), N, M)};
            else
                for i = 1:numel(varargin)
                    varargin{i} = varargin{i}.castJacToSparse();
                end
                h = cat(varargin{:});
            end
        end
        
        function AD = castJacToSparse(AD)
            for i = 1:numel(AD.jac)
                if ~issparse(AD.jac{i})
                    AD.jac{i} = AD.jac{i}.sparse();
                end
            end
        end
        
        function h = cat(varargin)
            if 1
                val = cellfun(@(x) x.val, varargin, 'Unif', false);
                jac = cellfun(@(x) horzcat(x.jac{:}), varargin, 'Unif', false);
                h = varargin{1};
                h.val = vertcat(val{:});
                h.jac = {vertcat(jac{:})};
            else
                n = nargin;
                nj = numel(varargin{1}.jac);
                njac = cellfun(@(x) matrixDims(x, 2), varargin{1}.jac);
                nv = cellfun(@(x) numel(x.val), varargin);
                nzeros = zeros(n, nj);
                
                vals = zeros(sum(nv), 1);
                offset = 0;
                for i = 1:n
                    nzeros(i, :) = cellfun(@nnz, varargin{i}.jac);
                    nloc = numel(varargin{i}.val);
                    vals(offset + (1:nloc)) = varargin{i}.val;
                    offset = offset + nloc;
                end
                [I, J, V] = deal(zeros(sum(sum(nzeros)), 1));
                shift = 0;
                
                csj = cumsum([0, njac]);
                csn = cumsum([0, nv]);
                for jacNo = 1:nj
                    for eqNo = 1:n
                        nz = nzeros(eqNo, jacNo);
                        offset = (shift+1):(shift + nz);
                        if 0
                            [ii, jj, vv, nl, ml] = getSparseArguments(varargin{eqNo}.jac{jacNo});
                            I(offset) = ii + sum(nv(1:eqNo-1));
                            J(offset) = jj + sum(njac(1:jacNo-1));
                            V(offset) = vv;
                        else
                            [I(offset), J(offset), V(offset)] = getSparseArguments(varargin{eqNo}.jac{jacNo}, ...
                                csn(eqNo), csj(jacNo));
                        end
                        shift = shift + nz;
                    end
                end
                h = varargin{1};
                h.jac = {sparse(I, J, V, sum(nv), sum(njac))};
                h.val = vals;
            end
        end

        function numVars = getNumVars(ad)
            numVars = ad.numVars;
        end
        function u = reduceToDouble(u)
            isZ = cellfun(@nnz, u.jac) == 0;
            if all(isZ)
                u = u.val;
            end
        end
        
        function h = rdivide(u,v)
            % Right element-wise division: `h = u./v`
            if ~isa(v, 'ADI')
                % Divide by double
                h = u;
                v_inv = 1./v;
                h.val = h.val.*v_inv;
                if numel(v_inv) == 1
                    for i = 1:numel(h.jac)
                        h.jac{i} = h.jac{i}*v_inv;
                    end
                else
                    h.jac = h.lMultDiag(v_inv, h.jac);
                end
            elseif isa(u, 'ADI')
                % Both are AD
                if numel(u.val) == numel(v.val)
                    uv = u.val;
                    vv = v.val;
                    h = u;
                    v_inv = 1./vv;
                    h.val = uv.*v_inv;
                    h.jac = h.timesJac(-h.val.*v_inv, v_inv, u.jac, v.jac);
                else
                    h = rdivide@ADI(u,v);
                end
            else
                % v is AD
                vv = v.val;
                h = v;
                h.val = u./vv;
                h.jac = h.lMultDiag(-u./(vv.^2), h.jac);
            end
        end
        
        function h = minus(u,v)
            if ~isa(u,'ADI') %u is a vector/scalar and v is ADI
                if numel(u) <= numel(v.val)
                    h = uminus(v);
                    h.val = h.val + u;
                elseif numel(v.val) == 1
                    h = minus(u, repmat(v,[numel(u), 1]));
                else
                    error('Vectors have different lengths')
                end
            elseif ~isa(v,'ADI')   %v is a vector/scalar and u is ADI
                if numel(v) <= numel(u.val)
                    h = u;
                    h.val = h.val - v;
                elseif numel(u.val) == 1
                    h = minus(repmat(u,[numel(v), 1]), v);
                else
                    error('Vectors have different lengths')
                end
            else
                % Both variables are ADI
                if numel(u.val) == numel(v.val)
                    h = u;
                    h.val = u.val - v.val;
                    h.jac = u.minusJac(h.jac, v.jac);
                elseif numel(u.val) == 1
                    h = minus(repmat(u, [numel(v.val), 1]), v);
                elseif numel(v.val) == 1
                    h = minus(u, repmat(v, [numel(u.val), 1]));
                else
                    error('Vectors have different lengths')
                end
            end
        end
    end
    
    methods (Static)

        %--------------------------------------------------------------------------
        
        function J = plusJac(J1, J2)
            nv1 = matrixDims(J1{1},1);
            nv2 = matrixDims(J2{1},1);
            if  nv1 == nv2
                J = J1;
                for i = 1:numel(J1)
                    J{i} = J{i} + J2{i};
                end
            else     % only other legal option is that nv1 = 1 or nv2 =1
                if nv1 == 1
                    J = cell(1, numel(J1));
                    for k = 1:numel(J)
                        J{k} = repmat(J1{k}, [nv2, 1]) + J2{k};
                    end
                else % nv2 = 1
                    assert(nv2 == 1)
                    J = GenericAD.plusJac(J2, J1);
                end
            end
        end
        
        %--------------------------------------------------------------------------
        
        function J = minusJac(J1, J2)
            nv1 = matrixDims(J1{1},1);
            nv2 = matrixDims(J2{1},1);
            if  nv1 == nv2
                J = J1;
                for i = 1:numel(J)
                    J{i} = J{i} - J2{i};
                end
            else     % only other legal option is that nv1 = 1 or nv2 =1
                if nv1 == 1
                    J = cell(1, numel(J1));
                    for k = 1:numel(J)
                        J{k} = repmat(J1{k}, [nv2, 1]) - J2{k};
                    end
                else % nv2 = 1
                    assert(nv2 == 1)
                    J = GenericAD.minusJac(J2, J1);
                end
            end
        end
        
        %--------------------------------------------------------------------------
        
        function J = mtimesJac(M, J1)
            J = cell(1, numel(J1));
            for k = 1:numel(J)
                J{k} = M*J1{k};
            end
        end
        
        %--------------------------------------------------------------------------
        
        function J = mtimesScalarJac(J1, J2)
            nv1 = size(J1{1},1);
            nv2 = size(J2{1},1);
            if nv1 == 1
                J = cell(1, numel(J1));
                for k = 1:numel(J)
                    J{k} = repmat(J1{k}, [nv2, 1])*J2{k};
                end
            elseif nv2 == 1
                J = GenericAD.mtimesScalarJac(J2, J1);
            else
                error('Not supported')
            end
        end
        
        %--------------------------------------------------------------------------
        
        function J = lMultDiag(d, J)
            D = [];
            for k = 1:numel(J)
                [J{k}, D] = diagMult(d, J{k}, D);
            end
        end
        
        %--------------------------------------------------------------------------
        
        function J = timesJacUnit(v_unit, v2, J_unit, J2)
            nj = numel(J_unit);
            J = cell(1, nj);
            for k = 1:nj
                J{k} = v_unit*J2{k} + sparse(v2)*J_unit{k};
            end
        end
        
        function J = timesJac(v1, v2, J1, J2)
            D1 = [];
            D2 = [];
            nj = numel(J1);
            anyPair1 = any(v1);
            anyPair2 = any(v2);
            if anyPair1 && anyPair2
                for k = 1:nj
                    [J1{k}, D1, D2] = diagProductMult(v1, v2, J1{k}, J2{k}, D1, D2);
                end
                J = J1;
            elseif anyPair1
                for k = 1:nj
                    [J1{k}, D1] = diagMult(v1, J2{k}, D1);
                end
                J = J1;
            else
                for k = 1:nj
                    [J2{k}, D2] = diagMult(v2, J1{k}, D2);
                end
                J = J2;
            end
        end
    end
end

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
