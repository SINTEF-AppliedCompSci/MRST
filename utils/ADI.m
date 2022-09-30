classdef ADI
    % ADI class: simple implementation of automatic differentiation for easy construction of jacobian matrices.
    %
    % SYNOPSIS:
    %   x = ADI(value, jacobian)
    %
    % PARAMETERS:
    %   value    - The numerical value of the object
    %
    %   jacobian - The Jacobian of the object.
    %
    % RETURNS:
    %   u - ADI object.
    %
    % NOTE:
    %  This class is typically instansiated for a set of different variables
    %  using `initVariablesADI`. The file contains a worked example
    %  demonstrating the usage for several variables.
    %
    % SEE ALSO:
    %   `initVariablesADI`, `PhysicalModel`


   properties
      val  % function value as a column vector of doubles
      jac  % cell array of sparse jacobian matrices
   end

   methods
      function obj = ADI(a,b)
         % ADI class constructor
         if nargin == 0 % empty constructor
            obj.val     = [];
            obj.jac     = {};
         elseif nargin == 1
             % Only allowed for values that are already ADI.
             if isa(a, 'ADI')
                 obj = a;
             else
                 error('Contructor requires 2 inputs')
             end
         elseif nargin == 2 % values + jacobians are supplied
             obj.val = a; % value
             if ~iscell(b)
                 b = {b};
             end
             obj.jac = b; % jacobian or list of jacobians
         else
             error('Input to constructor not valid')
         end
      end
      %--------------------------------------------------------------------
      function u = convertDouble(x, v)
          % Convert numeric type v into AD with zero derivatives with
          % respect to the same variables as x (which is already AD)
          assert(isa(v, 'double'));
          nval  = numel(v);
          nj = numel(x.jac);
          newjac = cell(1, nj);
          for i = 1:nj
              newjac{i} = sparse([], [], [], nval, size(x.jac{i}, 2));
          end
          u = x;
          u.val = v;
          u.jac = newjac;
      end

      %--------------------------------------------------------------------
      function h = numelValue(u)
          % Get number of values. Equivalent of `numel` for doubles. We do
          % not overload numel directly as it is not recommended by
          % Mathworks.
          h = numel(u.val);
      end

      %--------------------------------------------------------------------
      function h = double(u)
          % Cast to double and thereby remove derivatives:
          iname = inputname(1);
          if isempty(iname), iname = 'expression'; end

          warning('Future:Deprecation', ...
                 ['Method ADI/double may become deprecated. Use ', ...
                  '''value(%s)'' instead.'], iname);

          h = value(u);
      end

      %--------------------------------------------------------------------
      function h = value(u, varargin)
          % Cast to double and thereby remove derivatives:
          h = u.val;
      end
      %--------------------------------------------------------------------

      function h = ge(u, v)
          % Greater than or equal: `u>=v`
          h = ge(value(u), value(v));
      end

      %--------------------------------------------------------------------

      function h = gt(u, v)
          % Greater than: `u>v`
          h = gt(value(u), value(v));
      end

      %--------------------------------------------------------------------

      function h = le(u, v)
          % Less than or equal: `u<=v`
          h = le(value(u), value(v));
      end

      %--------------------------------------------------------------------

      function h = lt(u, v)
          % Less than: `u < v`
          h = lt(value(u), value(v));
      end
      %--------------------------------------------------------------------

      function h = uplus(u)
          % Unitary plus: `h = +u`
          h = u;
      end

      %--------------------------------------------------------------------

      function u = uminus(u)
          % Unitary minus: `h = -u`
         u.val = -u.val;
         u.jac =  u.uminusJac(u.jac);
      end

      %--------------------------------------------------------------------

      function h = plus(u,v)
          % Addition ot two values, where either value can be ADI of
          % appropriate dimensions: `h = u + v`
         if ~isa(u,'ADI') %u is a vector/scalar and v is ADI
             if numel(u) <= numel(v.val)
                 h = v;
                 h.val = h.val + u;
             elseif numel(v.val) == 1
                 h = plus(u, repmat(v,[numel(u), 1]));
             else
                 error('Vectors have different lengths')
             end
         elseif ~isa(v,'ADI')   %v is a vector/scalar and u is ADI
             if numel(v) <= numel(u.val)
                 h = u;
                 h.val = h.val + v;
             elseif numel(u.val) == 1
                 h = plus(repmat(u,[numel(v), 1]), v);
             else
                 error('Vectors have different lengths')
             end
         else
             % Both variables are ADI
             if numel(u.val) == numel(v.val)
                 h = u;
                 h.val = u.val + v.val;
                 h.jac = u.plusJac(h.jac, v.jac);
             elseif numel(u.val) == 1
                 h = plus(repmat(u, [numel(v.val), 1]), v);
             elseif numel(v.val) == 1
                 h = plus(u, repmat(v, [numel(u.val), 1]));
             else
                 error('Vectors have different lengths')
             end
         end
      end

      %--------------------------------------------------------------------

      function h = minus(u,v)
          % Subtraction with two elements: `h = u - v`
         h = plus(u, uminus(v));
      end

      %--------------------------------------------------------------------

      function h = mtimes(u,v)
          % Multiplication with matrix or scalar: `h = u*v`
          if ~isa(u,'ADI') %u is a scalar/matrix
              h = v;
              h.val = u*h.val;
              h.jac = h.mtimesJac(u, h.jac);
          elseif ~isa(v,'ADI') && isscalar(v)
              h = mtimes(v,u);
          else % special case where either u or v has single value
              if numel(u.val) == 1
                  h = u;
                  h.val = times(u.val, v.val);
                  h.jac = u.timesJacUnit(u.val, v.val, u.jac, v.jac);
              elseif numel(v.val) == 1
                  h = u;
                  h.val = times(u.val, v.val);
                  h.jac = u.timesJacUnit(v.val, u.val, v.jac, u.jac);
              else
                  error('Operation not supported');
              end
          end
      end

      %--------------------------------------------------------------------

      function h = times(u,v)
         % Element-wise multiplication: h = u.*v
         if ~isa(u,'ADI') %u is a scalar/vector and v is ADI.
             if numel(u)==numel(v.val)
                 h = v;
                 h.val = u.*h.val;
                 h.jac = h.lMultDiag(u, h.jac);
             else
                 h = mtimes(u,v);
             end
         elseif ~isa(v,'ADI') %v is a scalar/vector and u is ADI.
             h = times(v,u);
         else
             % Both u and v are ADI
             if numel(u.val)==numel(v.val)
                 h = u;
                 h.jac = h.timesJac(h.val, v.val, h.jac, v.jac);
                 h.val = h.val.*v.val;
             elseif numel(v.val)==1 || numel(u.val)==1
                 h = mtimes(u,v);
             else
                 error('Operation not supported');
             end
         end
      end

      %--------------------------------------------------------------------

      function h = mrdivide(u,v)
          % Right matrix divide: `h=u/v`
         if ~isa(v,'ADI') && isscalar(v) 
            h = mtimes(u, 1/v);
         else
            error('Operation not supported');
         end
      end

      %--------------------------------------------------------------------

      function h = mldivide(u,v)
          % Left matrix divide: `h=u\v`
          if ~isa(u,'ADI') %u is a scalar/matrix
              h.val = u\v.val;
              h.jac = h.mldivideJac(u, h.jac);
          else
              error('Operation not supported');
          end
      end

      %--------------------------------------------------------------------

      function h = power(u,v)
      % Element-wise power. `h=u.^v`.
          nu = numel(value(u));
          nv = numel(value(v));
          if nu == 1 && nv > 1
              u = repmat(u, [nv, 1]);
          end

          if nv == 1 && nu > 1
              v = repmat(v, [nu, 1]);
          end

          if ~isa(v,'ADI') % v is a scalar and u is ADI
              h = u;
              h.val = h.val.^v;
              h.jac = u.lMultDiag(v.*u.val.^(v-1), u.jac);
          elseif ~isa(u,'ADI') % u is a scalar and v is ADI
              h = v;
              h.val = u.^v.val;
              h.jac = v.lMultDiag((u.^v.val).*log(u), v.jac);
          else % u and v are both ADI
              h = u;
              h.val = u.val.^v.val;
              h.jac = h.plusJac( ...
                  h.lMultDiag((u.val.^v.val).*(v.val./u.val), u.jac), ...
                  h.lMultDiag((u.val.^v.val).*log(u.val),     v.jac) );
          end
      end

      %--------------------------------------------------------------------
      function h = polyval(p,v)
         if ~isa(p,'ADI') % p should be intigers and v adi
             h = v;
             h.val = polyval(p,v.val);
             h.jac = v.lMultDiag(polyval(polyder(p),v.val),v.jac);
             assert(isa(v,'ADI'));
         else
             error('polyval with adi coefficents not valid')
         end
      end
      
      %--------------------------------------------------------------------

      function h = rdivide(u,v)
          % Right element-wise division: `h = u./v`
          h = times(u, power(v, -1));
      end

      %--------------------------------------------------------------------

      function h = ldivide(u,v)
          % Left element-wise division: `h = u.\v`
          h = rdivide(v,u);
      end

      %--------------------------------------------------------------------
      function numVars = getNumVars(ad)
          % Get number of derivatives in each Jacobian block.
          numVars = cellfun(@(x) size(x, 2), ad.jac)';
      end

      function h = subsref(u,s)
          % Subscripted reference. Called for `h = u(v)`.
          if strcmp(s(1).type, '.')
              h = builtin('subsref',u,s);
          else
              switch s(1).type
                  case '()'
                      ns = numel(s(1).subs);
                      if ns > 1
                          if ns == 2
                              assert(s(1).subs{2} == 1 || ischar(s(1).subs{2}),...
                                  'Invalid indexing for AD-variable. Object can only be indexed as a column-vector');
                          else
                              error('More than 2 inputs recieved to AD variable.');
                          end
                      end
                      subs  = s(1).subs{1};
                      if ischar(s) && strcmp(subs, ':')
                          h = u;
                      else
                          if islogical(subs)
                              subs = find(subs);
                          end
                          h = u;
                          h.val = h.val(subs);
                          h.jac = h.subsrefJac(h.jac, subs);
                      end
                      if numel(s) > 1
                          % Recursively handle next operation
                          h = subsref(h, s(2:end));
                      end
                  case '{}'
                      error('Operation not supported');
              end
          end
      end

      %--------------------------------------------------------------------

      function u = subsasgn(u,s,v)
          % Subscripted reference. Called for `u(s) = v`
          if strcmp(s(1).type, '.')
              u = builtin('subsasgn',u,s,v);
          else
              switch s(1).type
                  case '()'
                      subs  = s.subs{:};
                      if ~isa(u, 'ADI') % u is a vector
                          warning('This place in the code is not reachable!!!')
                          u = double2AD(u, v.jac);
                      end
                      
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
      
      %--------------------------------------------------------------------

      function ix = end(u, k, n)
          assert(k == 1, 'ADI objects only support vector indexing.');
          ix = numel(u.val);
      end
      
      %--------------------------------------------------------------------

      function h = exp(u)
          % Element-wise exponential: `h=exp(u)`.
          eu = exp(u.val);
          h = u;
          h.val = eu;
          h.jac = h.lMultDiag(eu, u.jac);
      end
      %--------------------------------------------------------------------
      function h = log(u)
          % Element-wise natural logarithm: `h=log(u)`
          logu = log(u.val);
          h = u;
          h.val = logu;
          h.jac = h.lMultDiag(1./u.val, u.jac);
      end

      %--------------------------------------------------------------------
      function h = tanh(u)
          % Element-wise hyperbolic tangent: `h=tanh(u)`
          tanhu = tanh(u.val);
          h = u;
          h.val = tanhu;
          h.jac = h.lMultDiag(1./(cosh(u.val).^2), u.jac);
      end

      %--------------------------------------------------------------------
      function h = sin(u)
          % Element-wise sine: `h=sin(u)`
          sinu = sin(u.val);
          h = u;
          h.val = sinu;
          h.jac = h.lMultDiag(cos(u.val), u.jac);
      end
          
      %--------------------------------------------------------------------
      function h = cos(u)
          % Element-wise cosine: `h=cos(u)`
          cosu = cos(u.val);
          h = u;
          h.val = cosu;
          h.jac = h.lMultDiag(-sin(u.val), u.jac);
      end
         
      %-------------------------------------------------------------------- 
      function h = asin(u)
         asinu = asin(u.val);
         h = u;
         h.val = asinu;
         h.jac = h.lMultDiag(1./sqrt(1-u.val.^2), u.jac);
      end
      
      %-------------------------------------------------------------------- 
      function h = acos(u)
         acosu = acos(u.val);
         h = u;
         h.val = acosu;
         h.jac = h.lMultDiag(-1./sqrt(1-u.val.^2), u.jac);
      end
      
      %--------------------------------------------------------------------
      function h = atan(u)
         atanu = atan(u.val);
         h = u;
         h.val = atanu;
         h.jac = h.lMultDiag(1./(1 + u.val.^2), u.jac);
      end

      %--------------------------------------------------------------------

      function h = max(u,v)
          % Take the element-wise maximum value of two objects
          if(nargin==1)
              assert(isa(u,'ADI'));
              [value,i] = max(u.val);
              h = u;
              jacs  = h.subsrefJac(u.jac, i);
              h.val = value;
              h.jac = jacs;
              return;
          end
          assert(nargin==2, 'Max function implemented for up to 2 variables only.');
          if ~isa(u, 'ADI') % u is a DOUBLE
              value =  bsxfun(@max, u, v.val);
              inx   = ~bsxfun(@gt,  u, v.val) + 1; % Pick 'v' if u <= v
              h  = v;
              h.val = value;
              h.jac = h.lMultDiag(inx==2, v.jac);
          elseif ~isa(v,'ADI') %v is a vector
              h = max(v,u);
          else % both ADI, should have same number of values
              value = max(u.val, v.val);
              inx   = u.val > v.val;
              h = u;
              h.val = value;
              h.jac = h.plusJac(h.lMultDiag(inx, u.jac), h.lMultDiag(~inx, v.jac));
          end
      end

      %--------------------------------------------------------------------

      function h = min(u, v)
          % Element-wise minimum value of two objects.
          
          if(nargin==1)
              % Use def. of maximum to handle this
              h = -max(-u);
              return;
          end
          h = -max(-u, -v);
      end

      %--------------------------------------------------------------------

      function u = sum(u)
          % Sum of vector
          u.val = sum(u.val);
          u.jac = u.sumJac(u.jac);
      end

      %--------------------------------------------------------------------
      function u = cumsum(u)
          % Cumulative sum of vector
          u.val = cumsum(u.val);
          u.jac = u.cumsumJac(u.jac);
      end

      %--------------------------------------------------------------------

      function h = sign(u)
         % Element-wise sign of vector
         h = sign(u.val);
      end

      %--------------------------------------------------------------------

      function u = abs(u)
          % Absolute value
          u.jac = u.lMultDiag(sign(u.val), u.jac);
          u.val = abs(u.val);
      end

      %--------------------------------------------------------------------

      function u = repmat(u, varargin)
          % Replicate and tile array of values.
          % NOTE:
          %   Only allowed in the first (column) dimension for ADI objects.
          u.val = repmat(u.val, varargin{:});
          u.jac = u.repmatJac(u.jac, varargin{:});
      end

      %--------------------------------------------------------------------
      function h = vertcat(varargin)
          % Vertical concatentation of both ADI and double types.
          isD = cellfun(@isnumeric, varargin);
          if any(isD)
              sampleAD = varargin(~isD);
              sampleAD = sampleAD{1};
              for i = 1:numel(isD)
                  if isD(i)
                      varargin{i} = double2ADI(varargin{i}, sampleAD);
                  end
              end
          end
          nv    = numel(varargin);
          nj    = numel(varargin{1}.jac);
          vals  = cell(1,nv);
          jacs  = cell(1,nj);
          sjacs = cell(1, nj);
          for k = 1:nv
              vals{k} = varargin{k}.val;
          end
          for k = 1:nj
              for k1 = 1:nv
                  sjacs{k1} = varargin{k1}.jac{k};
              end
              jacs{k} = ADI.vertcatJac(sjacs{:});
          end
          h = varargin{1};
          h.val = vertcat(vals{:});
          h.jac = jacs;
      end

       function h = combineEquations(varargin)
          % Combine a set of equations of either ADI or double type to a
          % single equation with ADI type. The resulting equation will have
          % a single assembled sparse Jacobian containing all derivatives.
          isD = cellfun(@isnumeric, varargin);
          if any(isD)
              sampleAD = varargin(~isD);
              sampleAD = sampleAD{1};
              for i = 1:numel(isD)
                  if isD(i)
                      varargin{i} = double2ADI(varargin{i}, sampleAD);
                  end
              end
          end
          h = cat(varargin{:});
       end
      %--------------------------------------------------------------------

      function h = cat(varargin)
          h = vertcat(varargin{:});
          h.jac = {ADI.horzcatJac(h.jac{:})};
      end

      %--------------------------------------------------------------------

      function h = horzcat(varargin)
          if nargin > 1
              error('Operation horzcat is not supported for class ADI.');
          else
              h = varargin{1};
          end
      end

      %--------------------------------------------------------------------

      function u = interpReg(T, u, reginx)
          % Interpolate property with region support
          [y, dydu] = interpReg(T, u.val, reginx);
          u.val = y;
          u.jac = u.lMultDiag(dydu, u.jac);
      end

      function h = interpTable(X, Y, x, varargin)
          % Interpolate in a table
          h = x;
          [h.val, dyidx] = interpTable(X, Y, x.val, varargin{:});
          h.jac = h.lMultDiag(dyidx, x.jac);
      end
      
      function u = reduceToDouble(u)
          % Switch to double representation if no derivatives are present.
          if sum(cellfun(@nnz, u.jac)) == 0
              u = u.val;
          end
      end
      
      function u = subsetPlus(u, v, subs)
          if isa(u, 'ADI')
              u.val(subs) = u.val(subs) + value(v);
              if isa(v, 'ADI')
                  % Both are ADI. We need to update Jacobians
                  for i = 1:numel(u.jac)
                      u.jac{i} = subsetPlus(u.jac{i}, v.jac{i}, subs);
                  end
              end
          else
              % Adding ADI into double subset
              u = v.convertDouble(u);
              % Recursively handle this case
              u = subsetPlus(u, v, subs);
          end
      end
      
      function u = subsetMinus(u, v, subs)
          u = subsetPlus(u, -v, subs);
      end
   end
   
   methods (Static)
        %**************************************************************************
        %-------- Helper functions involving Jacobians  ---------------------------
        %**************************************************************************
        function J = uminusJac(J)
            for i = 1:numel(J)
                J{i} = -J{i};
            end
        end

        %--------------------------------------------------------------------------

        function J = plusJac(J1, J2)
        if isempty(J1)
            J = J2;
            return
        end
        if isempty(J2)
            J = J1;
            return
        end

        nv1 = size(J1{1},1);
        nv2 = size(J2{1},1);
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
                J = ADI.plusJac(J2, J1);
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
            J = ADI.mtimesScalarJac(J2, J1);
        else
            error('Not supported')
        end
        end

        %--------------------------------------------------------------------------

        function J = lMultDiag(d, J1)
        n = numel(d);
        if any(d)
            ix = (1:n)';
            D = sparse(ix, ix, d, n, n);
        else
            D = 0;
        end
        J = cell(1, numel(J1));
        for k = 1:numel(J)
            J{k} = D*J1{k};
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
        n = numel(v1);
        ix = (1:n)';
        D1 = sparse(ix, ix, v1, n, n);
        D2 = sparse(ix, ix, v2, n, n);

        nj = numel(J1);
        J = cell(1, nj);
        for k = 1:nj
            J{k} = D1*J2{k} + D2*J1{k};
        end
        end

        %--------------------------------------------------------------------------

        function J = mldivideJac(M, J1)
        J = cell(1, numel(J1));
        for k = 1:numel(J)
            J{k} = M\J1{k};
        end
        end

        %--------------------------------------------------------------------------

        function J = subsrefJac(J1, subs)
        J = cell(1, numel(J1));
        for k = 1:numel(J)
            J{k} = J1{k}(subs,:);
        end
        end

        %--------------------------------------------------------------------------
        function J = sumJac(J)
            J = cellfun(@(j1) sum(j1, 1), J, 'UniformOutput', false);
        end

        %--------------------------------------------------------------------------
        function J = cumsumJac(J1)
            J = cellfun(@(j1) cumsum(j1, 1), J1, 'UniformOutput', false);
        end

        %--------------------------------------------------------------------------

        function J = repmatJac(J1, varargin)
        J   = cell(1, numel(J1));
        if (varargin{end}(end)==1)
            for k = 1:numel(J)
                J{k} = repmat(J1{k}, varargin{:});
            end
        else
            error('Only vertical concatenation allowed for class ADI objects');
        end
        end

        %--------------------------------------------------------------------------

        function u = double2AD(u, J1)
        % u is vector, J reference jacobian
        nr = numel(u);
        J  = cell(1, numel(J1));
        for k = 1:numel(J)
            nc   = size(J1{k}, 2);
            J{k} = sparse(nr, nc);
        end
        u = ADI(u, J);
        end

        %--------------------------------------------------------------------------

        function J = subsasgnJac(J, subs, J1)
        if nargin == 3
            for k = 1:numel(J)
                J{k}(subs,:) = J1{k};
            end
        else
            for k = 1:numel(J)
                J{k}(subs,:) = 0;
            end
        end
        end

        %--------------------------------------------------------------------------

        function J = vertcatJac(varargin)
        J = vertcat(varargin{:});
        end

        %--------------------------------------------------------------------------

        function J = horzcatJac(varargin)
        J = horzcat(varargin{:});
        end
   end
end


%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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



%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
