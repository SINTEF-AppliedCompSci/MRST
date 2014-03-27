function obj= pvt2Table(fluid,p_range,T_range)
    properties={'density','enthalpy','viscosity'};
    [p,T]=meshgrid(p_range,T_range);
    for j=1:fluid
        for i=1:numel(properties)
            property=properties{i};
            val=fluid{j}.(property)(p(:),T(:));
            val=reshape(val,size(p));
            table=struct('x',p_range','y',T_range','data',val);
            obj{j}.(property) =@(pa,Ta) interp2DTable(table, pa, Ta);
        end
    end
end
%
