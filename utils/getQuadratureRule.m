function [x, w, nq] = getQuadratureRule(degree, dim)

    switch dim
        case 1
            
            degree = degree + (degree == 0);
            l  = legendrePolynomials(degree);
            dl = cellfun(@(l) dx(l,1), l, 'unif', false);
            dl = dl{end};
            
            if degree <= 1
                x = 0;
                
            elseif degree == 2
                x = [-sqrt(3); sqrt(3)];
                
            elseif degree == 3
                x = [-sqrt(3/5); 0; sqrt(3/5)];
                
            elseif degree == 4
                x = [-sqrt(3/7 + 2/7*sqrt(6/5)), -sqrt(3/7 - 2/7*sqrt(6/5)), ...
                      sqrt(3/7 - 2/7*sqrt(6/5)), -sqrt(3/7 + 2/7*sqrt(6/5))];
                
            end
            
            w = 2./((1 - x.^2).*dl(x).^2);
            
%             if degree <= 3
%                 x = [-1; 0; 1];
%                 degree = 3;
% 
%             elseif degree <= 4
%                 x = [-1; -sqrt(5)/5; sqrt(5)/5; 1];
%                 degree = 4;    
% 
%             elseif degree <= 5
%                 x = [-1; -sqrt(21)/7; 0; sqrt(21)/7; 1];
% 
%             elseif degree <= 6
%                 x = [-1; sqrt(1/21*(7+2*sqrt(7))); sqrt(1/21*(7-2*sqrt(7))); ...
%                          sqrt(1/21*(7-2*sqrt(7))); sqrt(1/21*(7+2*sqrt(7))); 1];
%                 degree = 6;
% 
%             end

%             w = 2./(degree*(degree-1)*legendre(x, degree-1).^2);
            
        case 2
            
            if degree <= 1
                w = 1;
                x = [0,0];
            
            elseif degree <= 2
                w = [1;1;1]/3;
                x = [0,1,1; 1,0,1; 1,1,0]/2;
            end
    end
    
    nq = numel(w);
    
end

% function p = legendre(x, l)
% 
%     p = 0;
%     
%     for k = 0:floor(l/2)
%         p = p + (-1).^k.*nchoosek(l, k)*nchoosek(2*l-2*k, l)*x.^(l-2*k);
%     end
%     p = p./(2^l);
%     
% end