function obj= pvt2Table(fluid,p_range,T_range)
    properties={'density','enthalpy','viscosity'};
    [p,T]=meshgrid(p_range,T_range);
    for j=1:numel(fluid)
        for i=1:numel(properties)
            property=properties{i};
            val=fluid{j}.(property)(p(:),T(:));
            val=reshape(val,size(p))';
            %{
            
            table=struct('x',p_range','y',T_range','data',val);
            obj{j}.(property) =@(pa,Ta) interp2DTable(table, pa, Ta);
            %}
            % the below seems to be faster if temprature vary less than
            % pressure
            table=struct('x',p_range','y',T_range','data',val);
            dT=diff(T_range);
            dP=diff(p_range);            
            if((max(dP)-min(dP))/max(dP) <1e-6 && (max(dT)-min(dT))/max(dT) <1e-6)% && false)
                table.dy=dT(1);
                table.dx=dP(1);
                %table2=struct('x',reshape(T_range,[],1),'y',reshape(p_range,[],1),'dx',dT(1),'dy',dP(1),'data',val');                
                obj{j}.(property) =@(pa,Ta) interpRegular2DTable(table, pa, Ta);
            else
                table2=struct('x',reshape(T_range,[],1),'y',reshape(p_range,[],1),'data',val');
                obj{j}.(property) =@(pa,Ta) interp2DTable(table, pa, Ta);
            end
        end
    end
end
%
