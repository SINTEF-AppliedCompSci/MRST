function states = addCO2HeightData(states, Gt, fluid)
% for each column under the caprock, this function computes the height of
% the free (mobile) CO2 plume, and the maximum height that was reached in
% the past by the migrating CO2 plume.

% Notes: The height of the free plume corresponds to a fluid column that is
% comprised of free CO2 and residually trapped water. The maximum height
% reached in the past corresponds to a fluid column that is comprised of
% the free plume portion in addition to a portion where CO2 is residually
% trapped. The residual trapping of CO2 occured when the CO2 plume migrated
% away, leaving behind residual CO2.

% See also:
%   addHeightData



    for i = 1:numel(states)
        
        % get the saturation of the CO2 in its free (mobile) state:
        states{i}.sCO2_free = free_sg(states{i}.s(:,2), states{i}.sGmax(:), fluid);
        
        % use the free CO2 saturation to determine the free plume height
        states{i}.h_free = Gt.cells.H.*states{i}.sCO2_free./(1-fluid.res_water);
        
        % use the maximum CO2 saturation which was reached in the past to
        % determine the maximum depth the CO2 has reached
        states{i}.h_max   = Gt.cells.H.*states{i}.sGmax(:)./(1-fluid.res_water);
        
    end
    
    
    
end