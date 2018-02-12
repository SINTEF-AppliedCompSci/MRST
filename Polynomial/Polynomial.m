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
%                 p.w = {0};
            elseif nargin == 1
                p.k = k;
                p.w = cell(size(k,1),1);
                p.w = ones(size(k,1),1);
%                 [p.w{:}] = deal(1);
            elseif nargin == 2
%                 assert(isa(w, 'cell'));
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
                        tmp = pw.*prod(repmat(x, numel(p(dNo)), 1).^pk,2);
                        val(:, dNo) = accumarray(repmat((1:size(x,1))', numel(p(dNo)), 1), tmp);
                    end
%                     for wNo = 1:numel(p)
%                         val = val + prod(x.^p.k(wNo,:),2).*p.w{wNo};
%                     end
                    [varargout{1:nargout}] = val;
                otherwise
                    [varargout{1:nargout}] = builtin('subsref', p, x);
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
        
        function n = numel(p)
            n = numel(p.w);
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
                r.w = cellfun(@(w) w*p, r.w, 'unif', false);
            elseif ~isa(q, 'Polynomial')
                r = p;
                r.w = cellfun(@(w) w*q, r.w, 'unif', false);
            else
                pk = repmat(p.k, numel(q), 1);
                qk = rldecode(q.k, numel(p)*ones(numel(q),1), 1);
                rk = qk + pk;
                
                pw = repmat(p.w, numel(q), 1);
                qw = rldecode(q.w, numel(p)*ones(numel(q),1), 1);
                rw = qw.*pw;
                
                
%                 rw = cell(numel(p)*numel(q),1);
%                 for i = 1:numel(p)
%                     for j = 1:numel(q)
%                         rw{(i-1)*numel(q) + j} = p.w{i}*q.w{j};
%                     end
%                 end
%                 
                r = Polynomial(rk, rw);
            end
        end
        
        function r = times(p,q)
            
        end
        
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