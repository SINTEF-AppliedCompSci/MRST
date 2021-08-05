function  [description, options, state0, model, schedule, plotOptions] = norne_simple_wo(varargin)

    % One-line description
    description ...
            = ['The Norne field model used in the paper:  ', ...
               'Flow Diagnostics for Model Ensembles , https://doi.org/10.3997/2214-4609.202035133'];
    % Module dependencies
    require ad-core ad-props ad-blackoil fd-publications
    options = struct('realization', 1);
    options = merge_options(options, varargin{:});
    if nargout <= 2, return; end
    % Get ensemble
    out      = setupNorneFn(options.realization); 
    state0   = out.state0;
    model    = out.model;
    fluid = initSimpleADIFluid('phases', 'WO'                        , ...
                               'mu'    , [   1,   5]*centi*poise     , ...
                               'rho'   , [1014, 859]*kilogram/meter^3, ...
                               'n'     , [   2,   2]                 );
    model = GenericBlackOilModel(model.G, model.rock, fluid, ...
                                 'gas'      , false          , ...
                                 'inputdata', model.inputdata);
    schedule = out.schedule;
%     schedule.step.val = rampupTimesteps(12*year, 100*day);
%     schedule.step.control = ones(numel(schedule.step.val), 1);
%     W = schedule.control(1).W;
%     isinj = [W.sign] > 0;
%     [W(isinj).lims] = deal(struct('bhp', 1500*barsa));
%     schedule.control(1).W = W;
    % Plotting
    plotOptions = {'View'              , [85, 30]  , ...
                   'PlotBoxAspectRatio', [1,1,0.4] , ...
                   'Size'              , [700, 600]};
end