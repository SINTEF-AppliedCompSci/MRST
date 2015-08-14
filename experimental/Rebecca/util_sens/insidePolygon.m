function a=insidePolygon(coord, x)
if(size(x,1)>1)
    a=nan(size(x,1),1);
    for i=1:size(x,1)
        a(i)=insidePolygon(coord,x(i,:));
    end
else
    a=bsxfun(@minus,coord,x);
    b=[a(2:end,:);a(1,:)];
    n=sqrt(sum(b.^2,2).*sum(a.^2,2));
    %angl=acos(bsxfun(@rdivide,sum(a.*b,2),n));
    angl=asin(bsxfun(@rdivide,cross2D(a,b),n));
    v=abs(sum(angl));% use abs so bouth orientation works
    if(abs(v-2*pi)<1e-7)
        a=true;
    else
        %assert(abs(v)<1e-3);
        a=false;
    end
end

end
function c=cross2D(a,b)
c=a(:,1).*b(:,2)-a(:,2).*b(:,1);
end