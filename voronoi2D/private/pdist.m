function d = pdist(x)
  % Equivalent to the matlab function pdist
  % Arguments:
  %    x         A set of points [x,y]
  % Returns:
  %    d         Distance, all to all
  % Halvor Møll Nielsen
    nc=size(x,1);
    d=zeros(1,nc*(nc-1)/2);
    nf=0;

    for i=1:nc-1
        dl=sqrt(sum(bsxfun(@minus,x(i,:),x(i+1:end,:)).^2,2));
        ndl=numel(dl);
        d(nf+1:nf+ndl)=dl;
        nf=nf+ndl;
    end
end

 
