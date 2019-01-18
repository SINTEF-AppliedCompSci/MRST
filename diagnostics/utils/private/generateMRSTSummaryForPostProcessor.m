function [wsdata] = readwellSolDataForPostProcessor(problem, timedatinfo)
%Undocumented utility function

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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
kwlist{1} = 'bhp';
kwlist{2} = 'qWs';
kwlist{3} = 'qOs';
kwlist{4} = 'qGs';
kwlist{5} = 'TIME';
kwlist{6} = 'TIMESTEPS';



% Generate data matrix


for i = 1:numel(states)
    bhp(:,i) = [states{i}.wellSol.bhp]';
    qWs(:,i) = [states{i}.wellSol.qWs]';
    qOs(:,i) = [states{i}.wellSol.qOs]';
    qGs(:,i) = [states{i}.wellSol.qGs]';
end

t = info.time';
timesteps = 1:1:numel(t);


sdata = [bhp;qWs;qOs;qGs;t;timesteps];

smry.KEYWORDS = kwlist;
smry.WGNAMES = namelist;
smry.STARTDAT = info.date;
smry.data = sdata;
smry.MRST = 1;
smry = addSmryFuncs(smry);

end
%--------------------------------------------------------------------------

function smry = addSmryFuncs(smry)
smry.get    = @(nm,kw,ms)getData(smry, nm, kw, ms);
smry.getKws = @(nm)getKeywords(smry, nm);
smry.getUnit= @(nm,kw)getUnit(smry, nm, kw);
end


function kws = getKeywords(smry,~)
kws = smry.KEYWORDS;
end

function s = getData(smry, nm, kw, ms)

nInx = find(strcmp(smry.WGNAMES, nm));  

switch kw 
    case 'bhp'
        kInx = 1;  
        s = smry.data(kInx*nInx,:);        
    case 'qWs'
        kInx = 2; 
        s = smry.data(kInx*nInx,:);
    case 'qOs'
        kInx = 3;
        s = smry.data(kInx*nInx,:);
    case 'qGs'
        kInx = 4;
        nInx = find(strcmp(smry.WGNAMES, nm));    
        s = smry.data(kInx*nInx,:);
    case 'TIME'
        s = smry.data(4*numel(smry.WGNAMES)+1,:);
    case 'TIMESTEPS'
        s = smry.data(4*numel(smry.WGNAMES)+2,:);       
end

end


function u = getUnit(smry, nm, kw)
switch kw 
    case 'bhp'
        u = 'Pascal';
    case 'qWs'
        u = 'm^3/s';
    case 'qOs'
        u = 'm^3/s';
    case 'qGs'
        u = 'm^3/s';
    case 'TIME'
        u = 'day';
    case 'TIMESTEPS'
        u = '-';
end
end




