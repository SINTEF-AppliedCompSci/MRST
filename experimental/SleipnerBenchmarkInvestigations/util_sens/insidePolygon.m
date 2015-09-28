function [a,v]=insidePolygon(coord, x)
if(size(x,1)>1)
    a=nan(size(x,1),1);
    v=nan(size(x,1),1);
    for i=1:size(x,1)
        [aa,vv]=insidePolygon(coord,x(i,:));
        a(i)=aa;
        v(i)=vv;
    end
else
    a=bsxfun(@minus,coord,x);
    a=[a(1:end,:);a(1,:)];
    a=bsxfun(@rdivide,a,sqrt(sum(a.^2,2)));
    b=[a(2:end,:);a(1,:)];
    %n=sqrt(sum(b.^2,2).*sum(a.^2,2));
    %angl=acos(bsxfun(@rdivide,sum(a.*b,2),n));
    %angl=asin(bsxfun(@rdivide,cross2D(a,b),n));
    angl=asin(cross2D(a,b));
    inprod=sum(a.*b,2);
    ind=inprod<0;% 
    angl(ind)=sign(angl(ind))*pi-angl(ind);
    %{
    ag
    %}

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