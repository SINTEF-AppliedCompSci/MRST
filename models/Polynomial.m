classdef Polynomial
    
    properties
        k
        w
        dim
    end
    
    methods
        
        function p = Polynomial(k, w)
            if nargin == 0
                p.k = [0,0,0];
                p.w = 1;
            elseif nargin == 1
                p.k = k;
                p.w = cell(size(k,1),1);
                p.w = ones(size(k,1),1);
            elseif nargin == 2
                p.k = k;
                p.w = w;
            else
                error('Contructor requires 1 or 2 inputs');
            end
            p.dim = size(k,2);
        end
        
        function varargout = subsref(p, x)
            
            switch x(1).type
                case '()'
                    x = x.subs{1};
                    
                    val = zeros(size(x,1), size(p,2));
                    for dNo = 1:size(p,2)
                        pk = reshape(repmat(p(dNo).k', size(x,1), 1), p(dNo).dim, [])';
                        pw = reshape(repmat(p(dNo).w', size(x,1), 1), 1, [])';
                        tmp = pw.*prod(repmat(x, numel(p(dNo).w), 1).^pk,2);
                        val(:, dNo) = accumarray(repmat((1:size(x,1))', numel(p(dNo).w), 1), tmp);
                    end
                    [varargout{1:nargout}] = val;
                    return
                otherwise
                    [varargout{1:nargout}] = builtin('subsref', p, x);
                    return
            end

        end
        
        function r = dx(p, i)
            r = p;
            r.w = r.w.*r.k(:,i);
            r.k(:,i) = max(r.k(:,i) - 1, 0);
        end
        function v = grad(p)
            for dNo = 1:p.dim
                v(dNo) = dx(p,dNo);
            end
        end
        
        function r = uplus(p)
            r = p;
        end
        
        function p = uminus(p)
            p.w = -p.w;
%             p.w = cellfun(@(w) -w, p.w, 'unif', false);
        end
        
        function r = plus(p,q)
            
            if ~isa(p, 'Polynomial')
                p = Polynomial(zeros(1, q.dim), {p});
            end
            if ~isa(q, 'Polynomial')
                q = Polynomial(zeros(1, p.dim), {q});
            end
            
            r = Polynomial([p.k; q.k], [p.w; q.w]);
            
        end
        
        function r = minus(p,q)
            r = plus(p,-q);
        end
        
        function r = mtimes(p,q)
            
            if ~isa(p, 'Polynomial')
                r = q;
                r.w = r.w*p;
            elseif ~isa(q, 'Polynomial')
                r = p;
                r.w = r.w*q;
            else
                pk = repmat(p.k, numel(q.w), 1);
                qk = rldecode(q.k, numel(p.w)*ones(numel(q.w),1), 1);
                rk = qk + pk;
                
                pw = repmat(p.w, numel(q.w), 1);
                qw = rldecode(q.w, numel(p.w)*ones(numel(q.w),1), 1);
                rw = qw.*pw;
                
                r = Polynomial(rk, rw);
            end
        end
        
        function r = tensorProduct(varargin)
            pdim = nargin;
            r = Polynomial(zeros(1,pdim),1);
            
            for pNo = 1:nargin
                p = varargin{pNo};
                pk = zeros(numel(p.w),pdim);
                pk(:, pNo) = p.k;
                q = Polynomial(pk, p.w);
                r = r*q;
            end
        end
                
                
        
        function r = rdivide(p,q)
            r = mtimes(p, power(q,-1));
        end
        
        
%         function r = times(p,q)
%             
%         end
        
        function r = power(p,q)
            if ~isa(q, 'Polynomial')
                r = p;
                r.k = r.k*q;
                r.w = cellfun(@(w) w.^q, r.w, 'unif', false);
            else
                error();
            end
        end
        
    end
    
end