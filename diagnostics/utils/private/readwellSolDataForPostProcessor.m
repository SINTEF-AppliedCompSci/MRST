function [wsdata] = readwellSolDataForPostProcessor(problem, varargin)
% Utility function to read wellSol data from simulated MRST packed problem 
% and return a structure in the correct format to be read into
% PostProcessDiagnostics.m

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

opt = struct('wellSolFields', [],...
             'startdate',    [0 0 0]);
opt = merge_options(opt, varargin{:});


% Get data from problem
wellSolHandler = problem.OutputHandlers.wellSols;
ws = wellSolHandler(:);


schedule = problem.SimulatorSetup.schedule;
for i = 1:numel(schedule.step.val)
   W{i} = schedule.control(schedule.step.control(i)).W;
end
W = W';


% Read names
% Assumes that all wellSols contain the same fields at the wellSol at
% timestep 1.
namelist = arrayfun(@(x) x.name, W{1}, 'UniformOutput', false);


% Read keywords
if isempty(opt.wellSolFields)
    keywordlist = fieldnames(ws{1});
else
    keywordlist = opt.wellSolFields;
end

% Generate data matrix

propIdx = 1;
for i = 1:numel(keywordlist)
    
    kw = keywordlist{i};
    
    try 
        dat = [ws{1}.(kw)];
    catch
            continue
    end
    
    if ~isa([ws{1}.(kw)],'double') || isempty([ws{1}.(kw)]) || numel([ws{1}.(kw)])~= numel(namelist)
        continue
    end
    
    wsdata.(kw) = zeros(numel(ws),numel(ws{1}));
    [nm{propIdx},unit] = getUnit(kw);
    
    for j = 1:numel(ws)
         data =  convertTo([ws{j}.(kw)],unit).*[ws{j}.status];
         if strcmp(nm{propIdx},'m^3/day')
             data = abs(data);
         end
         wsdata.(kw)(j,:) = data;
    end
    
    props{propIdx} = kw;
    propIdx = propIdx + 1;
    
end

t = cumsum(problem.SimulatorSetup.schedule.step.val(:))./day';
timesteps = 1:1:numel(t);


wsdata.props = props;
wsdata.wellNames = namelist;
wsdata.times = t;
wsdata.timesteps = timesteps;
wsdata.units = nm;

end

function [nm,u] = getUnit(kw)
switch kw 
    case 'bhp'
        nm  = 'barsa';
        u = barsa();
    case 'qWs'
        nm = 'm^3/day';
        u = 1./day();
    case 'qOs'
        nm = 'm^3/day';
        u = 1./day();
    case 'qGs'
        nm = 'm^3/day';
        u = 1./day();
   case 'qWr'
        nm = 'm^3/day';
        u = 1./day();
    case 'qOr'
        nm = 'm^3/day';
        u = 1./day();
    case 'qGr'
        nm = 'm^3/day';
        u = 1./day();  
    case 'qTr'
        nm = 'm^3/day';
        u = 1./day();      
    case 'qTs'
        nm = 'm^3/day';
        u = 1./day();       
    case 'time'
        nm = 'day';
        u = day();
    otherwise 
        nm = ''
        u = '-';
end
end




