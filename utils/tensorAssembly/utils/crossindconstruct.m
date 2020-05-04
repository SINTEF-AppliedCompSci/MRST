function [crossind1, crossind2] = crossindconstruct(nb1, nb2)
%
%
% SYNOPSIS:
%   function [crossind1, crossind2] = crossindconstruct(nb1, nb2)
%
% DESCRIPTION: Construct table of cross-product indices. It is required that
% nb1 and nb2 have same dimension. We explain the construction by two
% examples
% 
% 1) example of cross product when dim(nb1) = dim(nb2) = 1
% (it corresponds to a rldecode for crossind1 and repmat for crossind2)
%    
%    nb1 = [2], nb2 = [3];
%    ind1 = [1; 2]    % the definition of ind1 is  (1 : sum(nb1))'
%    ind2 = [1; 2; 3] % the definition of ind2 is  (1 : sum(nb2))'
%    
%    crossind1 = [1; 1; 1; 2; 2; 2];
%    crossind2 = [1; 2; 1; 2; 1; 2];
%    
% 2) dim(nb1) = dim(nb2) > 1: We repeat the previous construction for each pair
% in nb1 and nb2 and concacenate the result.
%
%    nb1 = [2; 3], nb2 = [2; 2];
%    ind1 = [1; 2; 3; 4; 5] % the definition of ind1 is  (1 : sum(nb1))'
%    ind2 = [1; 2; 3; 4]    % the definition of ind2 is  (1 : sum(nb2))'
%    
%    crossind1 = [1; 1; 1; 2; 2; 2; 3; 3; 4; 4; 5; 5];
%    crossind2 = [1; 2; 1; 2; 1; 2; 3; 4; 3; 4; 3; 4];
%
%
% PARAMETERS:
%   nb1 - integer vector gives number of repetition for first table
%   nb2 - integer vector gives number of repetition for second table
%
% RETURNS:
%   crossind1 - first index vector
%   crossind2 - second index vector
%
% EXAMPLE: `computeVagTrans`
%
% SEE ALSO:
%

    
    n = sum(nb1.*nb2);
    
    ci = rldecode(nb2, nb1);
    apos = [1; cumsum(ci) + 1];
    apos = apos(1 : (end - 1));
    crossind1 = zeros(n, 1);
    crossind1(apos) = 1;
    crossind1 = cumsum(crossind1);
    
    a = ones(n, 1);
    cpos = cumsum(ci) + 1;
    cpos = cpos(1 : (end - 1));
    cval = -ci;
    cval = cval(1 : (end - 1));
    c = zeros(n, 1);
    c(cpos) = cval;

    d = zeros(n, 1);
    dpos = cumsum(nb1.*nb2) + 1;
    dpos = dpos(1 : (end - 1));
    dval = nb2;
    dval = dval(1 : (end - 1));
    d(dpos) = dval;

    crossind2 = cumsum(a + c + d);
    
    
end
