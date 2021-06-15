function wellSol = initWellSolLocal(W, state0, wellSolInit)
% model = someFunction(state)
% initialization should depend on model, for now just dstinguish between 2
% and three phases

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

   wellSolGiven =  (nargin == 3);

   if size(state0.s, 2) == 2
      model = 'OW';
   else
      model = '3P';
   end

   if wellSolGiven
      wellSol = wellSolInit;
   elseif isfield(state0, 'wellSol')
      wellSol = extractMatchingWells(state0, W, model);
   else
      wellSol = defaultWellSol(state0, W, model);
   end
   wellSol = assignFromSchedule(W, wellSol);
end

function ws = defaultWellSol(state, W, model)

   nw = numel(W);
   ws = [];
   
   for k = 1 : nw
      ws = [ws; singleDefaultWell(state, W(k), model)];
   end

   if nw == 0
      % Create fake well. Could be made neater...
      W = [];
      W.cells = 1;
      W.name = 'FAKE_WELL';
      W.sign = 1;
      W.compi = [1; 0; 0];
      ws = singleDefaultWell(state, W, model);
      % Make wellSol structure empty (but the structure fields are kept, they are
      % needed when equations are assembled.)
      ws = repmat(ws, 0, 1);
   end

   % additional fields depending on model
   if isfield(state, 'c') % polymer model
      ws(1).poly = [];
   end
   
end

function ws = singleDefaultWell(state, W, model)

   if strcmp(model, 'OW')
      actPh = [1,2];
   else
      actPh = [1,2,3];
   end

   ws = struct(...
      'name',   [],...
      'status', [],...
      'type',   [],...
      'val',    [],...
      'sign',   [],...
      'bhp',    [],...
      'qTs',    [],...
      'qWs',    [],...
      'qOs',    [],...
      'qGs',    [],...
      'mixs',   [],...
      'cstatus',[],...
      'cdp',    [],...
      'cqs',    []);

   nConn = numel(W.cells);
   nPh   = numel(actPh);
   ws.name = W.name;
   % To avoid switching off wells, we need to start with a bhp that makes
   % a producer produce and an injector inject. Hence, we intitialize the
   % bhp such that the top connection pressure is 5bar above/below the
   % corresponding well-cell pressure. If W.dZ is ~= 0, however, we
   % don't know wht a decent pressure is ...
   % The increment should depend on the problem and the 5bar could be a
   % pit-fall... (also used in initializeBHP in updateConnDP)
   %if W.dZ(1) == 0
   ws.bhp = state.pressure(W.cells(1)) + 5*W.sign*barsa;
   %else
   %    ws.bhp = -inf;
   %end
   irate = eps;
   ws.qTs  = 0;
   ws.qWs  = W.sign*irate;
   ws.qOs  = W.sign*irate;
   ws.qGs  = W.sign*irate;
   ws.mixs = W.compi(actPh);
   ws.qs   = W.sign*ones(1, nPh)*irate;
   ws.cdp  = zeros(nConn,1);
   ws.cqs  = zeros(nConn,nPh);
end

function wellSol = extractMatchingWells(state, W, model)
   ws = state.wellSol;
   wellSol = [];
   if numel(W) == 0 
      wellSol = defaultWellSol(state, W, model);
   else
      
      [lia, locb] = ismember({W.name}, {ws.name});
      for iw = 1 : numel(W)
         % We compare only number of cell connections to identify the wells since the
         % position of the connections is not given in the state0.wellSol structure.
         if lia(iw) && (numel(ws(locb(iw)).cdp) ==  numel(W(iw).cells))
            wellSol = [wellSol; ws(locb(iw))];
         else
            wellSol = [wellSol; singleDefaultWell(state, W(iw), model)];
         end
      end
   end
end

function ws = assignFromSchedule(W, ws)
% set fields that should be updated if control has changed
   for k = 1:numel(W)
      ws(k).status  = W(k).status;
      ws(k).type    = W(k).type;
      ws(k).val     = W(k).val;
      ws(k).sign    = W(k).sign;
      ws(k).cstatus = W(k).cstatus;

      tp = W(k).type;
      if ws(k).status
         v  = W(k).val;
      else
         v = 0;
         ws(k).bhp = 0;
         ws(k).val = 0;
      end
      switch tp
        case 'bhp'
          ws(k).bhp = v;
        case 'rate'
          ws(k).qWs = v*W(k).compi(1);
          ws(k).qOs = v*W(k).compi(2);
          ws(k).qGs = v*W(k).compi(3);
        case 'orat'
          ws(k).qOs = v;
        case 'wrat'
          ws(k).qWs = v;
        case 'grat'
          ws(k).qGs = v;
      end % No good guess for qOs, etc...
   end
end


