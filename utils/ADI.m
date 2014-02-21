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
%   ADI object.
%
% COMMENTS:
%  This class is typically instansiated for a set of different variables
%  using initVariablesADI. The file contains a worked example demonstrating
%  the usage for several variables.
%
% SEE ALSO:
%   initVariablesADI

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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


   properties
      val  %function value
      jac  %list of sparse jacobian matrices
   end

   methods
      function obj = ADI(a,b)
         %ADI class constructor
         if nargin == 0 % empty constructor
            obj.val     = [];
            obj.jac     = {};
         elseif nargin == 1 %
             if isa(a, 'ADI')
                 obj = a;
             else
                 error('Contructor requires 2 inputs')
             end
         elseif nargin == 2 % values + jacobians
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
      function h = numval(u)
          h = numel(u.val);
      end

      %--------------------------------------------------------------------
      function h = double(u)
          h = u.val;
      end

      %--------------------------------------------------------------------

      function h = ge(u, v)
          h = ge(double(u), double(v));
      end

      %--------------------------------------------------------------------

      function h = gt(u, v)
          h = gt(double(u), double(v));
      end

      %--------------------------------------------------------------------

      function h = le(u, v)
          h = le(double(u), double(v));
      end

      %--------------------------------------------------------------------

      function h = lt(u, v)
          h = lt(double(u), double(v));
      end
      %--------------------------------------------------------------------

      function h = uplus(u)
          h = u;
      end

      %--------------------------------------------------------------------

      function h = uminus(u)
         h = ADI(-u.val, uminusJac(u.jac));
      end

      %--------------------------------------------------------------------

      function h = plus(u,v)
         if ~isa(u,'ADI')       %u is a vector/scalar
             if numel(u) <= numel(v.val)
                 h = ADI(u+v.val, v.jac);
             elseif numel(v.val) == 1
                 h = plus(u, repmat(v,[numel(u), 1]));
             else
                 error('Vectors have different lengths')
             end
         elseif ~isa(v,'ADI')   %v is a vector/scalar
             if numel(v) <= numel(u.val)
                 h = ADI(u.val + v, u.jac);
             elseif numel(u.val) == 1
                 h = plus(repmat(u,[numel(v), 1]), v);
             else
                 error('Vectors have different lengths')
             end
         else
             if numel(u.val) == numel(v.val)
                 h = ADI(u.val+v.val, plusJac(u.jac, v.jac) );
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
         h = plus(u, uminus(v));
      end

      %--------------------------------------------------------------------

      function h = mtimes(u,v)% '*'
          if ~isa(u,'ADI') %u is a scalar/matrix
              h = ADI(u*v.val, mtimesJac(u, v.jac));
          elseif ~isa(v,'ADI') %v is a scalar
              h = mtimes(v,u);
          else % special case where either u or v has single value
              if numel(u.val) == 1
                  h = times(repmat(u, [numel(v.val), 1]), v);
              elseif numel(v.val) == 1
                  h = times(u, repmat(v, [numel(u.val), 1]));
              else
                  error('Operation not supported');
              end
          end
      end

      %--------------------------------------------------------------------

      function h = times(u,v)% '.*'
         if ~isa(u,'ADI') %u is a scalar/vector
             if numel(u)==numel(v.val)
                 h = ADI(u.*v.val, lMultDiag(u, v.jac));
             else
                 h = mtimes(u,v);
             end
         elseif ~isa(v,'ADI') %v is a scalar/vector
             h = times(v,u);
         else
             if numel(u.val)==numel(v.val)
                 h = ADI(u.val.*v.val, timesJac(u.val, v.val, u.jac, v.jac));
             elseif numel(v.val)==1||numel(u.val)==1
                 h = mtimes(u,v);
             else
                 error('Operation not supported');
             end
         end
      end

      %--------------------------------------------------------------------

      function h = mrdivide(u,v)% '/'
         if ~isa(v,'ADI') %v is a scalar
            h = mtimes(u, 1/v);
         else
            error('Operation not supported');
         end
      end

      %--------------------------------------------------------------------

      function h = mldivide(u,v)% '\'
          if ~isa(u,'ADI') %u is a scalar/matrix
              h = ADI(u\v.val, mldivideJac(u, v.jac));
          else
              error('Operation not supported');
          end
      end

      %--------------------------------------------------------------------

      function h = power(u,v)% '.^'
         if ~isa(v,'ADI') % v is a scalar
            h = ADI(u.val.^v, lMultDiag(v.*u.val.^(v-1), u.jac));
         elseif ~isa(u,'ADI') % u is a scalar
            h = ADI(u.^v.val, lMultDiag((u.^v.val).*log(u), v.jac) );
         else % u and v are both ADI
            h = ADI(u.val.^v.val, plusJac( ...
               lMultDiag((u.val.^v.val).*(v.val./u.val), u.jac), ...
               lMultDiag((u.val.^v.val).*log(u.val),     v.jac) ) );
         end
      end

      %--------------------------------------------------------------------

      function h = rdivide(u,v)% './'
          h = times(u, power(v, -1));
      end

      %--------------------------------------------------------------------

      function h = ldivide(u,v)% '.\'
          h = rdivide(v,u);
      end

      %--------------------------------------------------------------------

      function h = subsref(u,s)
          switch s(1).type
              case '.'
                  h = builtin('subsref',u,s);
              case '()'
                  assert(numel(s(1).subs) == 1, ...
                      'Expected single index, got %d', numel(s(1).subs))
                  subs  = s(1).subs{1};
                  if ischar(s) && strcmp(subs, ':'),
                      h = u;
                  else
                      if islogical(subs), subs = find(subs); end
                      h = ADI(u.val(subs), subsrefJac(u.jac, subs));
                  end
                  if numel(s) > 1
                      % Recursively handle next operation
                      h = subsref(h, s(2:end));
                  end
              case '{}'
                  error('Operation not supported');
          end
      end

      %--------------------------------------------------------------------

      function u = subsasgn(u,s,v)
          switch s(1).type
              case '.'
                  u = builtin('subsasgn',u,s,v);
              case '()'
                  subs  = s.subs{:};
                  if ~isa(u, 'ADI') % u is a vector
                      warning('This place in the code is not reachable!!!')
                      u = double2AD(u, v.jac);
                  end
                  u.val(subs) = v;
                  if ~isa(v, 'ADI') % v is a constant vector
                      u.jac = subsasgnJac(u.jac, subs); % set rows to zero
                  else
                      u.jac = subsasgnJac(u.jac, subs, v.jac);
                  end
              case '{}'
                  error('Operation not supported');
          end
      end

      %--------------------------------------------------------------------

      function h = exp(u)
          eu = exp(u.val);
          h  = ADI(eu, lMultDiag(eu, u.jac));
      end
      %
      function h = log(u)
          logu = log(u.val);
          h  = ADI(logu, lMultDiag(1./u.val, u.jac));
      end
      %--------------------------------------------------------------------

      function h = max(u,v) % this function should be expanded
          if ~isa(u,'ADI') %u is a vector
              [value, inx] = max([u v.val], [], 2);
              h  = ADI(value, lMultDiag(inx==2, v.jac));
          elseif ~isa(v,'ADI') %v is a vector
              h = max(v,u);
          else
              error('Not yet implemented ...');
          end
      end

      %--------------------------------------------------------------------

      function h = sum(u)
         h = ADI(sum(u.val), sumJac(u.jac));
      end

      %--------------------------------------------------------------------
      function h = cumsum(u)
         h = ADI(cumsum(u.val), cumsumJac(u.jac));
      end

      %--------------------------------------------------------------------

      function h = sign(u)
         h = sign(u.val);
      end

      %--------------------------------------------------------------------

      function h = abs(u)
         h = ADI(abs(u.val), lMultDiag(sign(u.val), u.jac));
      end

      %--------------------------------------------------------------------

      function h = repmat(u, varargin)  % only makes sense if second dim =1
          h = ADI(repmat(u.val, varargin{:}), repmatJac(u.jac, varargin{:}));
      end

      %--------------------------------------------------------------------
      function h = vertcat(varargin)
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
              jacs{k} = vertcatJac(sjacs{:});
          end
          h = ADI(vertcat(vals{:}), jacs);
      end


      %--------------------------------------------------------------------

      function h = cat(varargin)
          h = vertcat(varargin{:});
          h = ADI(h.val, horzcatJac(h.jac{:}));
      end

      %--------------------------------------------------------------------

      function horzcat(varargin)
          error('horzcat doesn not make sense for class ADI')
      end

      %--------------------------------------------------------------------

      function h = interpReg(T, u, reginx)
          [y, dydu] = interpReg(T, u.val, reginx);
          h = ADI(y, lMultDiag(dydu, u.jac));
      end

      %--------------------------------------------------------------------

      function h = interpRegPVT(T, x, v, flag, reginx)

          if ~isa(x,'ADI') %u is a scalar/matrix
              [y, dydx, dydv] = interpRegPVT(T, x, v.val, flag, reginx);
              h = ADI(y, lMultDiag(dydx, v.jac));
          elseif ~isa(v,'ADI') %v is a scalar
              [y, dydx, dydv] = interpRegPVT(T, x.val, v, flag, reginx);
              h = ADI(y, lMultDiag(dydx, x.jac));
          else
              [y, dydx, dydv] = interpRegPVT(T, x.val, v.val, flag, reginx);
              h = ADI(y, timesJac(dydx, dydv, v.jac, x.jac)); %note order of input
          end

      end

      function h = interpTable(X, Y, x, varargin)
         y = interpTable(X, Y, x.val, varargin{:});
         dydx  = dinterpTable(X,Y, x.val, varargin{:});
         h = ADI(y,lMultDiag(dydx, x.jac));
      end

      %--------------------------------------------------------------------

%       function u = addToVals(u, inx, v)
%           % adds v to u(inx)
%           assert(numel(inx)==numel(v.val));
%           u.val(inx) = u.val(inx) + v.val;
%           for k = 1:numel(u.jac)
%               u.jac{k} = addToRows(u.jac{k}, inx, v.jac{k});
%           end
%       end
      %--------------------------------------------------------------------
   end
end

%**************************************************************************
%-------- Helper functions involving Jacobians  ---------------------------
%**************************************************************************

function J = uminusJac(J1)
J = cellfun(@uminus, J1, 'UniformOutput', false);
end

%--------------------------------------------------------------------------

function J = plusJac(J1, J2)
nv1 = size(J1{1},1);
nv2 = size(J2{1},1);
if  nv1 == nv2
    J = cellfun(@plus, J1, J2, 'UniformOutput', false);
else     % only other legal option is that nv1 = 1 or nv2 =1
    if nv1 == 1
        J = cell(1, numel(J1));
        for k = 1:numel(J)
            J{k} = repmat(J1{k}, [nv2, 1]) + J2{k};
        end
    else % nv2 = 1
        J = plusJac(J2, J1);
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
    J = mtimesScalarJac(J2, J1);
else
    error('Not supported')
end
end

%--------------------------------------------------------------------------

function J = lMultDiag(d, J1)
n = numel(d);
D = sparse((1:n)', (1:n)', d, n, n);
J = cell(1, numel(J1));
for k = 1:numel(J)
    J{k} = D*J1{k};
end
end

%--------------------------------------------------------------------------

function J = timesJac(v1, v2, J1, J2)
n = numel(v1);
D1 = sparse((1:n)', (1:n)', v1, n, n);
D2 = sparse((1:n)', (1:n)', v2, n, n);
J = cell(1, numel(J1));
for k = 1:numel(J)
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
function J = sumJac(J1)
J = cellfun(@sum, J1, 'UniformOutput', false);
end

%--------------------------------------------------------------------------
function J = cumsumJac(J1)
J = cellfun(@cumsum, J1, 'UniformOutput', false);
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

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
