function d = pdist2(x,y)
  % Calculate distance from one set to another
  % Arguments:
  %    x        Set of n points
  %    y        Set of m points
  % Returns:
  %    d        Matrix where d(i,j) is the distance from
  %             x(i,:) to y(j,:)
  % Uses the fact that 
  %    ||x-y||^2 = ||x||^2 + ||y||^2 - 2*x.y
  % and takes the absolute value to avoid problems with roundoff.
   d = sqrt(abs(bsxfun(@plus, sum(x.^2,2),sum(y.^2,2)') - 2*(x*y')));
end

 
