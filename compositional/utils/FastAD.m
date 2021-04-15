classdef FastAD
    % A very limited AD class for quick EOS assembly
    properties
        val
        jac
    end
    
    methods
        function v = FastAD(val, jac)
            v.val = val;
            v.jac = jac;
        end

        function h = rdivide(u,v)% './'
            h = times(u, power(v, -1));
        end
        
        function u = sum(u)
            u.val = sum(u.val);
            u.jac = sum(u.jac, 1);
        end
        
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
        
        function u = plus(u, v)
            if isa(u, 'FastAD') && isa(v, 'FastAD')
                u.val = u.val + v.val;
                u.jac = bsxfun(@plus, u.jac, v.jac);
            elseif isa(u, 'FastAD')
                u.val = u.val + v;
            else
                v.val = u + v.val;
                u = v;
            end
        end
        
        function u = minus(u, v)
            u = plus(u, -v);
        end
        
        function u = uminus(u)
            u.val = -u.val;
            u.jac = -u.jac;
        end

        function u = max(u, v)
            if nargin == 1
                [u.val, index] = max(u.val);
                u.jac = u.jac(index, :);
            else
                uLarger = double(u) > double(v);
                u = uLarger.*u + ~uLarger.*v;
            end
        end

        function h = log(u)
            logu = log(u.val);
            h  = FastAD(logu, bsxfun(@times, 1./u.val, u.jac));
        end
        
      function h = exp(u)
          eu = exp(u.val);
          h  = FastAD(eu, bsxfun(@times, eu, u.jac));
      end
        
        function h = times(u, v)
            if isa(u, 'FastAD') && isa(v, 'FastAD')
                vv = u.val.*v.val;
                jj = bsxfun(@plus, ...
                    bsxfun(@times, u.val, v.jac),  bsxfun(@times, v.val, u.jac));
            elseif isa(u, 'FastAD')
                vv = u.val.*v;
                jj = bsxfun(@times, u.jac, v);
            else
                vv = u.*v.val;
                jj = bsxfun(@times, v.jac, u);
            end
            h = FastAD(vv, jj);
        end
        
        function v = double(v)
            v = v.val;
        end
        
        function v = value(v)
            v = v.val;
        end
        
        function h = power(u,v)% '.^'
            if ~isa(v,'FastAD') % v is a scalar
                vv = u.val.^v;
                jj = bsxfun(@times, v.*u.val.^(v-1), u.jac);
            elseif ~isa(u,'FastAD') % u is a scalar
                vv = u.^v.val;
                jj = bsxfun(@times, vv.*log(u), v.jac);
            else % u and v are both ADI
                vv = u.val.^v.val;
                jj = bsxfun(@mtimes, vv.*(v.val./u.val), u.jac) +...
                     bsxfun(@mtimes, vv.*log(u.val), v.jac);
            end
            h = FastAD(vv, jj);
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
