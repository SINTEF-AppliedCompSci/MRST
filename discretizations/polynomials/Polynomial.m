classdef Polynomial
    % Polynomial class for basis functions in MRST
    properties
        k
        w
        dim
    end
    
    methods
        %-----------------------------------------------------------------%
        function poly = Polynomial(k, w)
            if nargin == 0
                poly.k = [0,0,0];
                poly.w = 1;
            elseif nargin == 1
                poly.k = k;
                poly.w = ones(size(k,1),1);
            elseif nargin == 2
                poly.k = k;
                poly.w = w;
            else
                error('Contructor requires 1 or 2 inputs');
            end 
            poly.dim = size(k,2);
        end
        
        %-----------------------------------------------------------------%
        function varargout = subsref(poly, x)
            switch x(1).type
                case '()'
                    x = x.subs{1};
                    npoly  = numel(poly);
                    npts   = size(x,1);
                    val    = zeros(npts, npoly);
                    for i = 1:numel(poly)
                        p = poly(i);
                        nterms = size(p.k,1);
                        v = 0;
                        for j = 1:nterms
                            v = v + p.w(j).*prod(x.^p.k(j,:),2);
                        end
                        val(:,i) = v;
                    end
                    [varargout{1:nargout}] = val;
                    return
                otherwise
                    [varargout{1:nargout}] = builtin('subsref', poly, x);
                    return
            end
        end
        
        %-----------------------------------------------------------------%
        function r = dx(p, i)
            r = p;
            r.w = r.w.*r.k(:,i);
            r.k(:,i) = max(r.k(:,i) - 1, 0);
        end
        
        %-----------------------------------------------------------------%
        function v = grad(p)
            for dNo = 1:p.dim
                v(dNo) = dx(p,dNo);
            end
        end
        
        %-----------------------------------------------------------------%
        function r = uplus(p)
            r = p;
        end
        
        %-----------------------------------------------------------------%
        function p = uminus(p)
            p.w = -p.w;
        end
        
        %-----------------------------------------------------------------%
        function r = plus(p,q)
            if ~isa(p, 'Polynomial')
                p = Polynomial(zeros(1, q.dim), {p});
            end
            if ~isa(q, 'Polynomial')
                q = Polynomial(zeros(1, p.dim), {q});
            end
            r = Polynomial([p.k; q.k], [p.w; q.w]);
        end
        
        %-----------------------------------------------------------------%
        function r = minus(p,q)
            r = plus(p,-q);
        end
        
        %-----------------------------------------------------------------%
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
        
        %-----------------------------------------------------------------% 
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
           
        %-----------------------------------------------------------------%
        function r = rdivide(p,q)
            r = mtimes(p, power(q,-1));
        end
        
        %-----------------------------------------------------------------%
        function r = power(p,q)
            if ~isa(q, 'Polynomial')
                r = p;
                r.k = r.k*q;
                r.w = cellfun(@(w) w.^q, r.w, 'unif', false);
            else
                error('Not supported');
            end
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
