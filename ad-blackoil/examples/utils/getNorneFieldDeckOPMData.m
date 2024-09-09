function [deck, output] = getNorneFieldDeckOPMData()
%Undocumented Utility Function

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    pth = getDatasetPath('norne_field_model');
    opm = fullfile(pth,'opm-data');
    opm_tests = fullfile(pth,'opm-tests');
    name = 'norne';
    caseName = 'NORNE_ATW2013';
  
    fn = fullfile(opm, name, [caseName, '.DATA']);
    deck = readEclipseDeck(fn);
    
    units = {'METRIC', 'FIELD', 'SI', 'LAB', 'PVT_M'};
    
    unit = 'SI';
    for i = 1:numel(units)
        if isfield(deck.RUNSPEC, units{i})
            unit = units{i};
            break
        end
    end
    
    deck = convertDeckUnits(deck);
    fldr = fullfile(opm_tests, name);
    output = struct('opm', [], 'ecl', []);
    try
        output.opm = getOutputOPM(deck, fldr, caseName, unit);
    catch e
        fprintf('Error in parsing of OPM output. Reported error:\n %s\n',e.message);
    end
    try
        output.eclipse = getOutputEclipse(deck, fldr, caseName, unit);
    catch e
        fprintf('Error in parsing of Eclipse output. Reported error:\n %s\n',e.message);
    end
    
end


function output = getOutputOPM(deck, fldr, caseName, unit)
    G = initEclipseGrid(deck, 'SplitDisconnected', false);
    smry_prefix = fullfile(fldr, 'opm-simulation-reference', 'flow', caseName);
    if ~exist([smry_prefix, '.UNSMRY'], 'file')
        smry_prefix = fullfile(fldr, 'opm-simulation-reference', 'flow_legacy', caseName);
    end

    [ws, wstime] = convertSummaryToWellSols(smry_prefix, unit);
    states = convertRestartToStates(smry_prefix, G, 'use_opm', true, 'consistentWellSols', false);

    output = struct('wellSols', {ws}, 'states', {states}, 'G', G, 'location', smry_prefix, 'time', wstime);
end

function output = getOutputEclipse(deck, fldr, caseName, unit)
    G = initEclipseGrid(deck);
    smry_prefix = fullfile(fldr, 'eclipse-simulation', caseName);

    if ~exist([smry_prefix, '.SMSPEC'], 'file')
        % Norne hack...
        smry_prefix = fullfile(fldr, 'ECL.2014.2', caseName);
    end

    ws = convertSummaryToWellSols(smry_prefix, unit);
    try
        states = convertRestartToStates(smry_prefix, G, 'use_opm', false, 'consistentWellSols', false);
    catch e
        fprintf('Restartfiles not pressent. Reported error:\n %s\n',e.message);
        states=[];
    end

    output = struct('wellSols', {ws}, 'states', {states}, 'G', G);
end
