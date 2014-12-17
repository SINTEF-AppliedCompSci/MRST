function  makePvtTable(ffname,ff,p_range,T_range)
    [p,T]=meshgrid(p_range,T_range);
    val=ff(p(:),T(:));
    val=reshape(val,size(p));
    table=struct('x',reshape(p_range,[],1),'y',reshape(T_range',[],1),'data',val);
    save(ffname,'table');
end
%
