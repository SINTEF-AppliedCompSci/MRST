function state = updateAdsorption(state0, state, model)
%
%
% SYNOPSIS:
%   function state = updateAdsorption(state0, state, model)
%
% DESCRIPTION: Update the adsorption value in the state variable. Used by the
% surfactant models.
%
% PARAMETERS:
%   state0 - State at previous time step
%   state  - State at current time state
%   model  - Instance of the model
%
% RETURNS:
%   state - State at current time state with updated adsorption values.
%
% EXAMPLE:
%
% SEE ALSO: ExplicitConcentrationModel, FullyImplicitOilWaterSurfactantModel, ImplicitExplicitOilWaterSurfactantModel
%


    c = model.getProps(state, 'surfactant');
    adsmax0 = model.getProps(state0, 'adsmax');

    ads = computeEffAds(c, adsmax0, model.fluid);
    adsmax = max(ads, adsmax0);
    state = model.setProp(state, 'ads', ads);
    state = model.setProp(state, 'adsmax', adsmax);

end


