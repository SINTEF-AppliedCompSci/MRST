function [t, v] = collectEnsembleGUIValues(g, name, varargin)
% Utility function to extract diagnostics values from EnsembleGUI-object
%
% SYNOPSIS:
%    [t,v]  = getSaturationFront(g, valueName, 'pn1', pv1, ...)
% 
% REQUIRED PARAMETERS:
%  g    - ensembleGUI-object.
%  name - name (string) of diagnostics value (as given in drop-down menues for 
%        cross-plot or line-plot)
%
%OPTIONAL PARAMETERS
%  'injectors' - vector of injector-indicies defining interaction region 
%              (defqult all)
%  'producers' - vector of producers-indicies defining interaction region 
%              (defqult all)
%  'members'   - vector of member (realization)-indicies (default all) 
%  'time'      - vector of times (at given unit) for which values will be
%               output:
%                * value of scalar type (e.g., 'Lorenz coefficient', 'Sweep', ...):
%                   - time cannot be defaulted
%                * value of dynamic type (e.g., 'Sweep vs time', 'Oil rates', ...)
%                   - if defaulted, output all computed values, otherwise interpolated
%  'timeUnit' - 'years' or 'pvi' (deafult: 'years')
% 
% 
% RETURNS:
%  t   - time-matrix. Either single column or one column per member
%  v   - value-matrix. One column per member.

opt = struct('injectors', (1:numel(g.injectorNames)), ...
             'producers', (1:numel(g.producerNames)), ...
             'members',   g.validMembers, ...
             'time',      [], ...
             'timeUnit',  []);
         
opt = merge_options(opt, varargin{:});

% check if requested output is scalar type 
isScalar = any(strcmp(name, g.crossPlotSelections));
if ~isScalar
    assert(any(strcmp(name, g.linePlotSelections)), ...
        'Requested values of unknown type: %s', name);
    assert(~strcmp(name, 'F-Phi'), ...
        'This function does not currently output F-Phi-values')
end

% setup GUI according to input
tmUnit = opt.timeUnit;
if isempty(tmUnit)
    warning('Time-unit is defaulted to ''years''');
    tmUnit = 'years';
end
assert(any(strcmp(tmUnit, {'years', 'pvi'})), ...
       'Unknown time-unit: %s', tmUnit);
tmItem = g.Menu.items{5};
tmItem.switchUnit(tmUnit);

iregItem = g.Menu.items{2};
iregItem.communicationLimit = 0; %don't threshold communication
iregItem.autoDetect = false;     %we set region explicitly
iregItem.injectorIx = opt.injectors;
iregItem.producerIx = opt.producers;

if isScalar
    t = opt.time(:);
    assert(~isempty(t), 'Cannot output values for empty time');
    v = nan(numel(t), numel(opt.members));
    for k = 1:numel(t)
        tmItem.Value = t(k);
        % recompute (maybe not needed but be sure)
        g.emptyCurrentDiagnostics();
        tmp = g.getDiagnosticsValues(name);
        v(k,:) = [tmp{opt.members}];
    end
else % dynamic/vector-output
    % recompute (maybe not needed but be sure)
    g.emptyCurrentDiagnostics();
    %(1) time-years (2) time pvi, (3) value
    tmp = g.getDiagnosticsValues(name);
    if strcmp(tmUnit, 'years')
        tix = 1;
    else
        tix = 2;
    end
    % for empty input time, give full vectors
    if isempty(opt.time)
        t = {tmp{tix}{opt.members}};
        v = {tmp{3}{opt.members}};
    else % interpolate to given times
        t = opt.time(:);
        v = nan(numel(t), numel(opt.members));
        for k = 1:numel(opt.members)
            mix = opt.members(k);
            v(:,k) = interp1(tmp{tix}{mix}, tmp{3}{mix}, t, 'linear', nan);
        end
    end
end
end
    
        

        