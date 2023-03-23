function setup = olympus_field_wo(varargin)
% Setup function for waterflooding simulation on the OLYMPUS reservoir
%
% SYNOPSIS:
%   setup = olympus_field_wo('pn1', pv1, ...)
%   setup = olympus_field_wo(fullSetup, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Setup of a waterflooding simulation on the OLYMPUS reservoir model from
%   the ISAPP challenge (https://www.isapp2.com/optimization-challenge/).
%
%   Configurable parameter:
%      'realization' - realization number (default: 1, range: 1 to 50)
%
% RETURNS:
%   setup - test case with the following fields: name, description,
%      options, state0, model, schedule, and plotOptions.
%      If the optional input fullSetup (see synopsis) is false, the
%      returned setup only contains name, description, and options.
%
% SEE ALSO:
%   TestCase, testcase_template, testSuiteTutorial.

    % Description
    description = ['OLYMPUS field model from the ISAPP challenge; ', ...
        'see https://www.isapp2.com/optimization-challenge'];

    % Optional input arguments
    options = struct('realization', 1);
    
    [options, fullSetup, setup] = processTestCaseInput(mfilename, ...
        options, description, varargin{:});
    if ~fullSetup, return; end
    
    % Define module dependencies
    require ad-core ad-props ad-blackoil

    % Setup initial state, model, and schedule
    gravity reset on
    fldr     = fullfile(mrstDataDirectory, 'olympus', 'ECLIPSE');
    caseName = ['OLYMPUS_', num2str(options.realization)];
    deck     = fullfile(fldr, caseName, [caseName, '.DATA']);    
    [state0, model, schedule] = initEclipseProblemAD(deck, 'useMex', true);
    
    % Fix schedule to use more reasonable time steps:   
    % The input conversion only understands the DATES keyword and hence
    % ends up with one-year time steps. We split each of these time steps
    % in twelve to get closer to the maximum time step of 30 days specified
    % by the TUNING keyword.
    ix  = repmat(1:numel(schedule.step.val), [12 1]);
    schedule.step.val = schedule.step.val(ix(:))/12;
    schedule.step.control = schedule.step.control(ix(:));

    % Plotting
    plotOptions = {'View'              , [-10, 40]  , ...
                   'Size'              , [900, 450], ...
                   'PlotBoxAspectRatio', [33.7444   16.8243    5.0000], ...
                   'Box'               , true};
    
    % Pack setup
    setup = packTestCaseSetup(mfilename,                  ...
                              'description', description, ...
                              'options'    , options    , ...
                              'state0'     , state0     , ...
                              'model'      , model      , ...
                              'schedule'   , schedule   , ...
                              'plotOptions', plotOptions);
end