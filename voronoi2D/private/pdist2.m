function d = pdist2(x,y)
  % Calculate distance from one set to another
  % Arguments:
  %    x        Set of n points
  %    y        Set of m points
  % Returns:
  %    d        Matrix where d(i,j) is the distance from
  %             x(i,:) to y(j,:)
  % Written by Halvor Møll Nilsen
      d=nan(size(x,1),size(y,1));
    for i=1:size(y,1)
        dl=sqrt(sum(bsxfun(@minus,x,y(i,:)).^2,2));
        d(:,i)=dl;
    end
end

 
