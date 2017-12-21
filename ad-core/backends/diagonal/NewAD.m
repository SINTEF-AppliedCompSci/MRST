classdef NewAD < ADI
    % NewAD is the testbed for future updates to the ADI class. All
    % features herein are subject to rapid change.
    properties
        numVars
        offsets
    end
    methods
        function ad = NewAD(varargin)
            ad = ad@ADI(varargin{:});
        end
        
      %--------------------------------------------------------------------

      function h = plus(u,v)
         if ~isa(u,'ADI')       %u is a vector/scalar
             if numel(u) <= numel(v.val)
                 h = v;
                 h.val = h.val + u;
             elseif numel(v.val) == 1
                 h = plus(u, repmat(v,[numel(u), 1]));
             else
                 error('Vectors have different lengths')
             end
         elseif ~isa(v,'ADI')   %v is a vector/scalar
             if numel(v) <= numel(u.val)
                 h = u;
                 h.val = h.val + v;
             elseif numel(u.val) == 1
                 h = plus(repmat(u,[numel(v), 1]), v);
             else
                 error('Vectors have different lengths')
             end
         else
             if numel(u.val) == numel(v.val)
                 h = u;
                 h.val = u.val + v.val;
                 h.jac = NewAD.plusJac(h.jac, v.jac);
                 if isempty(h.jac)
                     h = h.val;
                 end
             elseif numel(u.val) == 1
                 h = plus(repmat(u, [numel(v.val), 1]), v);
             elseif numel(v.val) == 1
                 h = plus(u, repmat(v, [numel(u.val), 1]));
             else
                 error('Vectors have different lengths')
             end
         end
      end
      function h = mtimes(u,v) % '*'
          if ~isa(u,'ADI') %u is a scalar/matrix
              h = v;
              h.val = u*h.val;
              h.jac = NewAD.mtimesJac(u, h.jac);
              if isempty(h.jac)
                  h = h.val;
              end
          elseif ~isa(v,'ADI') %v is a scalar
              h = mtimes(v,u);
          else % special case where either u or v has single value
              if numel(u.val) == 1
                  h = u;
                  h.val = times(u.val, v.val);
                  h.jac = NewAD.timesJacUnit(u.val, v.val, u.jac, v.jac);
              elseif numel(v.val) == 1
                  h = u;
                  h.val = times(u.val, v.val);
                  h.jac = NewAD.timesJacUnit(v.val, u.val, v.jac, u.jac);
              else
                  error('Operation not supported');
              end
          end
      end
      function h = times(u,v) % '.*'
         if ~isa(u,'ADI') %u is a scalar/vector
             if numel(u)==numel(v.val)
                 h = v;
                 h.val = u.*h.val;
                 h.jac = NewAD.lMultDiag(u, h.jac);
             else
                 h = mtimes(u,v);
             end
         elseif ~isa(v,'ADI') %v is a scalar/vector
             h = times(v,u);
         else
             if numel(u.val)==numel(v.val)
                 h = u;
                 h.jac = NewAD.timesJac(h.val, v.val, h.jac, v.jac);
                 h.val = h.val.*v.val;
             elseif numel(v.val)==1 || numel(u.val)==1
                 h = mtimes(u,v);
             else
                 error('Operation not supported');
             end
         end
      end
      function h = power(u,v) % '.^'
         if ~isa(v,'ADI') % v is a scalar
             h = u;
             h.val = h.val.^v;
             h.jac = NewAD.lMultDiag(v.*u.val.^(v-1), u.jac);
         elseif ~isa(u,'ADI') % u is a scalar
             h = v;
             h.val = u.^v.val;
             h.jac = NewAD.lMultDiag((u.^v.val).*log(u), v.jac);
         else % u and v are both ADI
             h = u;
             h.val = u.val.^v.val;
             h.jac = NewAD.plusJac( ...
               NewAD.lMultDiag((u.val.^v.val).*(v.val./u.val), u.jac), ...
               NewAD.lMultDiag((u.val.^v.val).*log(u.val),     v.jac) );
         end
      end
      
      
      function u = subsasgn(u,s,v)
          if strcmp(s(1).type, '.')
              u = builtin('subsasgn',u,s,v);
          else
              switch s(1).type
                  case '()'
                      subs  = s.subs{:};
                      if isa(v, 'ADI') % v is a ADI
                          u.jac = u.subsasgnJac(u.jac, subs, v.jac);
                          u.val(subs) = v.val;
                      else
                          u.jac = u.subsasgnJac(u.jac, subs); % set rows to zero
                          u.val(subs) = v;
                      end
                  case '{}'
                      error('Operation not supported');
              end
          end
      end
      
      function x = incrementSubset(x, subs, v)
          if isa(x, 'NewAD')
              x.val(subs) = x.val(subs) + double(v);
              if isa(v, 'NewAD')
                  for i = 1:numel(x.jac)
                      x.jac{i} = incrementSubset(x.jac{i}, subs, v.jac{i});
                  end
              end
          else
              % V is AD
              x(subs) = x(subs) + v;
          end
      end
      
      function h = exp(u)
          eu = exp(u.val);
          h = u;
          h.val = eu;
          h.jac = NewAD.lMultDiag(eu, u.jac);
      end
      %
      function h = log(u)
          logu = log(u.val);
          h = u;
          h.val = logu;
          h.jac = NewAD.lMultDiag(1./u.val, u.jac);
      end
      function h = max(u,v) % this function should be expanded
          if(nargin==1)
              assert(isa(u,'ADI'));
              [value,i] = max(u.val);
              jacs      = u.subsrefJac(u.jac, i);
              h = u;
              h.val = value;
              h.jac = jacs;
              return;
          end
          assert(nargin==2, 'Max function implemented for up to 2 variables only.');
          if ~isa(u, 'ADI'), % u is a DOUBLE
              value =  bsxfun(@max, u, v.val);
              inx   = ~bsxfun(@gt,  u, v.val) + 1; % Pick 'v' if u <= v
              h  = v;
              h.val = value;
              h.jac = NewAD.lMultDiag(inx==2, v.jac);
             if isempty(h.jac)
                 h = h.val;
             end
          elseif ~isa(v,'ADI') %v is a vector
              h = max(v,u);
          else % both ADI, should have same number of values
              value = max(u.val, v.val);
              inx   = u.val > v.val;
              h = u;
              h.val = value;
              h.jac = NewAD.plusJac(NewAD.lMultDiag(inx, u.jac),NewAD.lMultDiag(~inx, v.jac));
              if isempty(h.jac)
                  h = h.val;
              end
          end
      end
      function u = abs(u)
          u.jac = NewAD.lMultDiag(sign(u.val), u.jac);
          u.val = abs(u.val);
          if isempty(u.jac)
               u = u.val;
          end
      end

      function h = vertcat(varargin)
          isD = cellfun(@isnumeric, varargin);
          if any(isD)
              sampleAD = varargin(~isD);
              sampleAD = sampleAD{1};
              for i = 1:numel(isD)
                  if isD(i)
                      varargin{i} = double2NewAD(varargin{i}, sampleAD);
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
              jacs{k} = NewAD.vertcatJac(sjacs{:});
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
                      varargin{i} = double2NewAD(varargin{i}, sampleAD);
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
              h = vertcat(varargin{:});
              h.jac = {h.horzcatJac(h.jac{:})};
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

      function u = interpReg(T, u, reginx)
          [y, dydu] = interpReg(T, u.val, reginx);
          u.val = y;
          u.jac = NewAD.lMultDiag(dydu, u.jac);
         if isempty(u.jac)
             u = u.val;
         end
      end
      function h = interpRegPVT(T, x, v, flag, reginx)

          if ~isa(x,'ADI') %u is a scalar/matrix
              h = v;
              [h.val, dydx] = interpRegPVT(T, x, v.val, flag, reginx);
              h.jac = NewAD.lMultDiag(dydx, v.jac);
          elseif ~isa(v,'ADI') %v is a scalar
              h = x;
              [h.val, dydx] = interpRegPVT(T, x.val, v, flag, reginx);
              h.jac = NewAD.lMultDiag(dydx, x.jac);
          else
              h = x;
              [h.val, dydx, dydv] = interpRegPVT(T, x.val, v.val, flag, reginx);
              h.jac = NewAD.timesJac(dydx, dydv, v.jac, x.jac); % note order of input
          end
         if isempty(h.jac)
             h = h.val;
         end
      end

      function h = interpTable(X, Y, x, varargin)
         h = x;
         h.val = interpTable(X, Y, x.val, varargin{:});
         h.jac = NewAD.lMultDiag(dinterpTable(X,Y, x.val, varargin{:}), x.jac);
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

    end
    
    methods (Access=protected, Static)

        function J = uminusJac(J1)
        J = cellfun(@uminus, J1, 'UniformOutput', false);
        end

        %--------------------------------------------------------------------------

        function J = plusJac(J1, J2)
        nv1 = matrixDims(J1{1},1);
        nv2 = matrixDims(J2{1},1);
        if  nv1 == nv2
            J = cellfun(@plus, J1, J2, 'UniformOutput', false);
        else     % only other legal option is that nv1 = 1 or nv2 =1
            if nv1 == 1
                J = cell(1, numel(J1));
                for k = 1:numel(J)
                    J{k} = repmat(J1{k}, [nv2, 1]) + J2{k};
                end
            else % nv2 = 1
                assert(nv2 == 1)
                J = NewAD.plusJac(J2, J1);
            end
        end
        end


        %--------------------------------------------------------------------------

        function J = mtimesJac(M, J1)
        if nnz(M) == 0
            J = {};
            return
        end
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
            J = NewAD.mtimesScalarJac(J2, J1);
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
        J = cell(1, nj);
        anyPair1 = any(v1);
        anyPair2 = any(v2);

        if anyPair1 && anyPair2
            for k = 1:nj
                [J{k}, D1] = diagMult(v1, J2{k}, D1);
                [tmp2, D2] = diagMult(v2, J1{k}, D2);
                J{k} = J{k} + tmp2;
            end
        elseif anyPair1
            for k = 1:nj
                [J{k}, D1] = diagMult(v1, J2{k}, D1);
            end
        elseif anyPair2
            for k = 1:nj
                [J{k}, D1] = diagMult(v2, J1{k}, D1);
            end
        else
            for k = 1:nj
                if issparse(J1{k})
                    J{k} = 0*J1{k};
                else
                    J{k} = J1{k}.toZero();
                end
            end
        end
        end
    end
end