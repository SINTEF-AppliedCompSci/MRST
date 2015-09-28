function states = addHeightData(states,Gt,fluid)
    for i=1:numel(states)
       states{i}.h=Gt.cells.H.*states{i}.s(:,2)./(1-fluid.res_water); 
    end
end