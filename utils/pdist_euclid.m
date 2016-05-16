function dist = pdist_euclid(x)
% Pairwise euclidian distance between pairs of objects.

assert(~isempty(x) && numel(size(x))<=2, 'x must be a 1 or 2 dimensional non-empty matrix');
ind = nchoosek(1:size(x,1),2);
Xi = ind(:,1);
Yi = ind(:,2);
X = x';
diff = X(:,Xi) - X(:,Yi);
if numel(diff) == 2
    dist = norm(diff);
else
    dist = sqrt(sumsq(diff,1));
end
return

function s = sumsq(x,dim) % 2D matrices only
if dim == 1
    s = zeros(1,size(x,2));
    for i = 1:size(x,2)
        s(i) = norm(x(:,i))^2;
    end
elseif dim == 2
    s = zeros(size(x,1),1);
    for i = 1:size(x,1)
        s(i) = norm(x(i,:))^2;
    end
else
    error('Input must be a 2D matrix and dim must be 1 or 2');
end
return