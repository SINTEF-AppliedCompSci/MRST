classdef FastAD
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
        
        function u = plus(u, v)
            if isa(u, 'FastAD') && isa(v, 'FastAD')
                u.val = u.val + v.val;
                u.jac = u.jac + v.jac;
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
        
        function h = power(u,v)% '.^'
            if ~isa(v,'FastAD') % v is a scalar
                vv = u.val.^v;
                jj = bsxfun(@times, v.*u.val.^(v-1), u.jac);
            elseif ~isa(u,'FastAD') % u is a scalar
                vv = u.^v.val;
                jj = bsxfun(@times, vv.*log(u), v.jac);
            else % u and v are both ADI
                vv = u.val.^v.val;
                jj = bsxfun(@plus, bsxfun(@times, vv.*(v.val./u.val), u.jac), ...
                                   bsxfun(@times, vv.*log(u.val), u.jac));
            end
            h = FastAD(vv, jj);
        end
    end
end