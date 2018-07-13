function suitable = aquiferConditionSuitability(p,t, varargin)
% Determine whether the aquifer conditions yeild suitable conditions for storage
%
% NB: t is in kelvin, p is in pascals
%
% "Dense phase" CO2 is liquid, not gas (when temp is below critical point)
% 
% Returns: 1 for dense phase (liquid) or above supercritical, 0 for
% non-dense phase (gas) or below supercritical

opt.plot = false;
opt = merge_options(opt, varargin{:});

[pc, tc] = CO2CriticalPoint();

% can work when p and t are cell arrays


    % get corresponding vapour-pressure for the input temperature (i.e.,
    % the pressure on that liquid-vapour boundary for that temperature).
    % NB: if t > tc, NaN is returned.
    pv = CO2VaporPressure(t);

    % is the input pressure greater than or less than that boundary pressure?
    %   - if greater than, co2 is Liquid (suitable)
    %   - if equal to or less than, co2 is Gas (non-suitable)
    suitable = p > pv;

    % if nan entries exist, we check whether or not input pressure is above
    % the critical pressure (73.8 bars).
    %   - if yes, conditions are supercritical (suitable)
    %   - if no, conditions are not supercritical (non-suitable)
    cells2check = find(isnan(pv));
    suitable(cells2check) = p(cells2check) > pc;
    
% optional plot:
if opt.plot
    
    % show input press and temp on diagram that also shows critical point
    % and liquid-vapour boundary
    figure, hold on
    
    % critical point
    plot(tc-273.15, convertTo(pc,barsa), 'o')
    xlabel('Temperature (C)')
    ylabel('Pressure (bar)')
    
    % liquid-vapour boundary
    tv = [273.15:1:tc];
    pv = CO2VaporPressure(tv);
    plot(tv-273.15, convertTo(pv,barsa), '+')
    
    % input values
    plot(t-273.15, convertTo(p,barsa), 'x')
    
end


end