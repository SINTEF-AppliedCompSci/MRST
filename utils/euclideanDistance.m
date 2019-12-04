function D = euclideanDistance(X, Y, varargin)
% Calculate euclidean distance from one set to another
% Equivalent to the matlab function pdist2(X, Y, 'euclidean')
% See https://stackoverflow.com/questions/7696734/pdist2-equivalent-in-matlab-version-7
    D = bsxfun(@plus,sum(X.^2,2),sum(Y.^2,2)') - 2*(X*Y') ;
    D = abs(D);
    D = sqrt(D);
end