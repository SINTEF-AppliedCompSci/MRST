function W = processWellsDP(G, rock_matrix,rock_fracture, control, varargin)
%Construct MRST well structure from ECLIPSE input deck control.
%
% SYNOPSIS:
%   W = processWells(G, rock, control)
%   W = processWells(G, rock, control, 'pn1', pv1, ...)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%   rock    - Rock data structure.  Must contain valid field 'perm'.
%
%   control - Well control mode.  Assumed to be defined from function
%             'readSCHEDULE'.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%           - InnerProduct -- The inner product with which to define
%                             the mass matrix.
%                             String.  Default value = 'ip_tpf'.
%             Supported values are:
%               - 'ip_simple'
%               - 'ip_tpf'
%               - 'ip_quasitpf'
%               - 'ip_rt'
% 
%           - Verbose      -- Whether or not to emit informational messages
%                            if any requested well control mode is
%                            unsupported.
%                            Logical.  Default value = false.
%           - DepthReorder -- Boolean indicating if we should attempt to
%                             reorder the perforations by depth. Should
%                             only be used for vertical wells. Default off.
%           - strictParsing -- Throw errors if unsupported well controls
%                              are encountered. If disabled, wells with
%                              unsupported controls will recieve NaN values
%                              for control and phase composition.
%           - createDefaultWell -- If enabled, this function will not parse
%                              any wells, but simply return a default well
%                              structure suitable for further modification.

%           - DualPoro          true or false for processing DP wells
%
% RETURNS:
%   W  - Updated (or freshly created) well structure, each element of which
%        has the following fields:
%        cells   -- Grid cells perforated by this well.
%        type    -- Well control type.
%        val     -- Target control value.
%        r       -- Well bore radius.
%        dir     -- Well direction.
%        WI      -- Well index.
%        dZ      -- Vertical displacement of each well perforation measured
%                   from 'highest' horizontal contact (i.e., the 'TOP'
%                   contact with the minimum 'Z' value counted amongst all
%                   cells perforated by this well).
%       name     -- Well name.
%       compi    -- Fluid composition--only used for injectors.
%       sign     -- injection (+) or production (-) flag.
%       status   -- Boolean indicating if the well is open or shut.
%       cstatus  -- One entry per cell, indicating if the completion is
%                   open.
%       lims     -- Limits for the well. Contains subfields for the types
%                   of limits applicable to the well (bhp, rate, orat, ...)
%                   Injectors generally have upper limits, while producers
%                   have lower limits.
%       c        -- Component concentrations. Unit convention may vary.
%
% SEE ALSO:
%   readSCHEDULE, readWellKW, addWell.

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


   opt = struct('InnerProduct',     'ip_tpf', ...
                'Verbose',           mrstVerbose(), ...
                'DepthReorder',      false, ...
                'strictParsing',     true, ...
                'createDefaultWell', false, ...
                'DualPoro',             false, ...
                'DualPoroSplit',        false );
   opt = merge_options(opt, varargin{:});

   if opt.createDefaultWell
      W = createDefaultWell(G, rock);
      return;
   end

   
   if opt.DualPoro %This works but there is probably a better way to do this
       
       COMPDAT4=cell2mat(control.COMPDAT(:,4));
       COMPDAT5=cell2mat(control.COMPDAT(:,5));
       
        fractureindex=[];
        matrixindex=[];
        
       for i=1:numel(COMPDAT4)
                if COMPDAT4(i)>G.cartDims(3)
                    control.COMPDAT{i,4}=COMPDAT4(i)-G.cartDims(3);
                    control.COMPDAT{i,5}=COMPDAT5(i)-G.cartDims(3);
                 
                   
                    control.COMPDAT{i,15}='Fracture';
                    fractureindex=[fractureindex; i];
                    
                else
                    control.COMPDAT{i,15}='Matrix';
                    matrixindex=[matrixindex; i];

                end
       end    
   end

   if opt.DualPoro && ~opt.DualPoroSplit
       
       if ~isempty(fractureindex)
           %fracture perfs
            controlF=control;
            controlF.COMPDAT=control.COMPDAT(fractureindex,:);
            well_id      = enumerateWells(controlF);
            controlF      = insertDefaults(controlF, well_id, G.cartDims(3));
            [controlF, pf] = orderCompletions(controlF);
       end
        
       if ~isempty(matrixindex)
           %Matrix perfs
            controlM=control;
            controlM.COMPDAT=control.COMPDAT(matrixindex,:)
            well_id      = enumerateWells(controlM);
            controlM     = insertDefaults(controlM, well_id, G.cartDims(3));
            [controlM, pm] = orderCompletions(controlM);
       end
   else
       well_id      = enumerateWells(control);
       control     = insertDefaults(control, well_id, G.cartDims(3));
       [control, p] = orderCompletions(control);
       
   end
   


   inj  = struct('WCONINJ' , @process_wconinj , ...
                 'WCONINJE', @process_wconinje, ...
                 'WCONINJH', @process_wconinjh);

   prod = struct('WCONHIST', @process_wconhist, ...
                 'WCONPROD', @process_wconprod);

   post = struct('WPOLYMER', @process_wpolymer, ...
                 'WSURFACT', @process_wsurfact);

   % ----------------------------------------------------------------------
   % Well processing stages
   %
   stages = { inj, prod, post };

   W = [];

   if  opt.DualPoro && ~opt.DualPoroSplit
       
    %fracture wells
    if ~isempty(fractureindex)
       for stage = reshape(stages, 1, []),
          fn = reshape(sort(fieldnames(stage{1})), 1, []);

          for kw = fn(isfield(controlF, fn)),
             W = stage{1}.(kw{1}) (W, controlF, G, rock_fracture, well_id, pf, opt);
          end
       end

       if isempty(W)
          W = createDefaultWell(G, rock_fracture);
          return;
       end

       if opt.DepthReorder
          W = reorderWellPerforationsByDepth(W);
       end

       % for now, connection status can only be true if well status is true  
       for k = 1:numel(W)
           W(k).cstatus = and(W(k).status, W(k).cstatus);
           % if all completions are closed, well must be off
           W(k).status = and(W(k).status, any(W(k).cstatus));
       end
    end
       
             
    % for matrix wells
    if ~isempty(matrixindex)
       for stage = reshape(stages, 1, []),
          fn = reshape(sort(fieldnames(stage{1})), 1, []);

          for kw = fn(isfield(controlM, fn)),
             W = stage{1}.(kw{1}) (W, controlM, G, rock_matrix, well_id, pm, opt);
          end
       end

       if isempty(W)
          W = createDefaultWell(G, rock_matrix);
          return;
       end

       if opt.DepthReorder
          W = reorderWellPerforationsByDepth(W);
       end

       % for now, connection status can only be true if well status is true  
       for k = 1:numel(W)
           W(k).cstatus = and(W(k).status, W(k).cstatus);
           % if all completions are closed, well must be off
           W(k).status = and(W(k).status, any(W(k).cstatus));
       end
    end
       
       
       
   end
   
   
 if  opt.DualPoro && opt.DualPoroSplit
  %split the rock by flowing and stagnant porosity
    
  
      %Assume rock_fracture is the flowing rock porosity 
       for stage = reshape(stages, 1, []),
          fn = reshape(sort(fieldnames(stage{1})), 1, []);

          for kw = fn(isfield(control, fn)),
             W = stage{1}.(kw{1}) (W, control, G, rock_fracture, well_id, p, opt);
          end
       end

       if isempty(W)
          W = createDefaultWell(G, rock_fracture);
          return;
       end

       if opt.DepthReorder
          W = reorderWellPerforationsByDepth(W);
       end

       % for now, connection status can only be true if well status is true  
       for k = 1:numel(W)
           W(k).cstatus = and(W(k).status, W(k).cstatus);
           % if all completions are closed, well must be off
           W(k).status = and(W(k).status, any(W(k).cstatus));
       end
       

 end
 
end


    
    

%--------------------------------------------------------------------------
% Private helpers follow.
%--------------------------------------------------------------------------

function W = process_wconinj(W, control, G, rock, well_id, p, opt)
   for i = 1 : size(control.WCONINJ, 1),
      nm = control.WCONINJ{i,1};
      status = strcmp(control.WCONINJ{i,3}, 'OPEN');

      type = lower(control.WCONINJ{i,4});
      switch type,
        case 'rate', val = control.WCONINJ{i, 5};
        case 'resv', val = control.WCONINJ{i, 6};
        case 'bhp' , val = control.WCONINJ{i, 9};

        otherwise
          dispif(opt.Verbose, ...
                 ['Control mode ''%s'' unsupported for injector ', ...
                  '''%s''.  Ignored.\n'], upper(type), nm);
          continue
      end

      % Define injection fluid.  This is a hack.
      switch lower(control.WCONINJ{i,2}(1)),
        case 'w', compi = [1, 0, 0];  % Water, 1st phase
        case 'o', compi = [0, 1, 0];  % Oil, 2nd phase
        case 'g', compi = [0, 0, 1];  % Gas, 3rd phase
        otherwise,
          dispif(opt.Verbose, ...
                 ['Injection phase ''%s'' is unknown.  ', ...
                  'Well ignored.\n'], control.WCONINJ{i,2});
          continue
      end

      sizeW = numel(W);
      W = buildWell(W, G, rock, control, well_id(nm), p, type, val, ...
                    compi, opt.InnerProduct, 1);
      if numel(W) > sizeW
         W(end).lims.rate = control.WCONINJ{i, 5};
         W(end).lims.bhp  = control.WCONINJ{i, 9};
         W(end).status    = status;
      end
   end
end

%--------------------------------------------------------------------------

function W = process_wconinje(W, control, G, rock, well_id, p, opt)
   for i = 1 : size(control.WCONINJE, 1),
      nm = control.WCONINJE{i,1};
      status = strcmp(control.WCONINJE{i,3}, 'OPEN');

      type = lower(control.WCONINJE{i,4});
      switch type,
        case 'rate'
          val = control.WCONINJE{i, 5};
        case 'resv'
          dispif(opt.Verbose, ...
                 ['Control mode ''%s'' unsupported for injector ', ...
                  '''%s''.  Ignored.\n'], upper(type), nm);
          continue
          % val = control.WCONINJE{i, 6};
        case 'bhp' , val = control.WCONINJE{i, 7};
        otherwise
          dispif(opt.Verbose, ...
                 ['Control mode ''%s'' unsupported for injector ', ...
                  '''%s''.  Ignored.\n'], upper(type), nm);
          continue
      end


      % Define injection fluid.  This is a hack.
      switch lower(control.WCONINJE{i,2}(1)),
        case 'w', compi = [1, 0, 0];  % Water, 1st phase
        case 'o', compi = [0, 1, 0];  % Oil, 2nd phase
        case 'g', compi = [0, 0, 1];  % Gas, 3rd phase
        otherwise,
          dispif(opt.Verbose, ['Injection phase ''%s'' is unknown.  ', 'Well ignored.\n'], ...
                 control.WCONINJE{i,2});
          continue
      end

      sizeW = numel(W);
      W = buildWell(W, G, rock, control, well_id(nm), p, type, val, ...
                    compi, opt.InnerProduct, 1);

      if numel(W) > sizeW
         W(end).lims.rate = control.WCONINJE{i, 5};
         W(end).lims.bhp  = control.WCONINJE{i, 7};
         W(end).status    = status;
      end
   end
end

%--------------------------------------------------------------------------

function W = process_wconinjh(W, control, G, rock, well_id, p, opt)
   for i = 1 : size(control.WCONINJH, 1),
      nm = control.WCONINJH{i,1};
      status = strcmp(control.WCONINJH{i,3}, 'OPEN');

      type = 'rate'; % History matched well are always rate controlled
      val = control.WCONINJH{i, 4};

      % Define injection fluid.  This is a hack.
      switch lower(control.WCONINJH{i,2}(1)),
        case 'w', compi = [1, 0, 0];  % Water, 1st phase
        case 'o', compi = [0, 1, 0];  % Oil, 2nd phase
        case 'g', compi = [0, 0, 1];  % Gas, 3rd phase
        otherwise,
          dispif(opt.Verbose, ...
                 ['Injection phase ''%s'' is unknown.  ', ...
                  'Well ignored.\n'], control.WCONINJH{i,2});
          continue
      end

      sizeW = numel(W);
      W = buildWell(W, G, rock, control, well_id(nm), p, type, val, ...
                    compi, opt.InnerProduct, 1);
      if numel(W) > sizeW
         W(end).lims.rate = control.WCONINJH{i, 4};
         W(end).lims.bhp  = 6895 * barsa;  % \approx 1e5 psia, as for default wconinj
         W(end).status    = status;
      end
   end
end

%--------------------------------------------------------------------------

function W = process_wpolymer(W, control, G, rock, well_id, p, opt)
   npoly = size(control.WPOLYMER, 1);
   if npoly == 0
       return;
   end
   if ~isfield(W, 'c')
       [W.c] = deal([]);
   end
   for i = 1:numel(W)
       % Add zeros for polymer
       W(i).c = [W(i).c, 0];
   end
   for i = 1 : npoly
      for j = 1:size(W,1)
         if strcmp(W(j).name, control.WPOLYMER{i,1})
            W(j).c(end) = control.WPOLYMER{i,2};
         end
      end
   end
end

%--------------------------------------------------------------------------

function W = process_wsurfact(W, control, varargin)
   nsurf = size(control.WSURFACT, 1);
   if nsurf == 0
       return;
   end
   if ~isfield(W, 'c')
       [W.c] = deal([]);
   end
   for i = 1:numel(W)
       % Add zeros for surfactant
       W(i).c = [W(i).c, 0];
   end
   if ~isempty(W),
      Wn = { W.name };

      for i = 1 : nsurf
         j = find(strcmp(Wn, control.WSURFACT{i,1}));

         if ~isempty(j),
            [ W(j).c(end) ] = deal(control.WSURFACT{i,2});
         end
      end
   end
end

%--------------------------------------------------------------------------

function W = process_wconhist(W, control, G, rock, well_id, p, opt)
   for i = 1 : size(control.WCONHIST, 1),
      nm = control.WCONHIST{i,1};
      status = strcmp(control.WCONHIST{i,2}, 'OPEN');

      type = lower(control.WCONHIST{i,3});
      switch type
        case 'orat'
          val = - control.WCONHIST{i, 4};
          compi = [0, 1, 0];  % Oil, 2nd phase
        case 'wrat'
          val = - control.WCONHIST{i, 5};
          compi = [1, 0, 0];  % Water, 1st phase
        case 'grat'
          val = - control.WCONHIST{i, 6};
          compi = [0, 0, 1];  % Gas, 3rd phase
        case 'lrat'
          val = - sum([control.WCONHIST{i, 4:5}]);
          compi = [1, 1, 0];
        case 'resv'
          if opt.strictParsing
             error('RESV for WCONHIST not supported');
          else
             warning(['RESV for WCONHIST not supported. '...
                      'Inserting NaN values for control in well %s.'], nm);
          end
          val = nan;
          compi = [nan, nan, nan];
        case 'bhp'
          val = control.WCONHIST{i, 10};

        otherwise
          dispif(opt.Verbose, ...
                 ['Control mode ''%s'' unsupported for producer ', ...
                  '''%s''.  Ignored.\n'], upper(type), nm);
          continue
      end

      sizeW = numel(W);

      W = buildWell(W, G, rock, control, well_id(nm), p, type, val, ...
                    compi, opt.InnerProduct, -1);

      if numel(W) > sizeW
         switch type,
           case 'orat'
             W(end).lims.orat = -control.WCONHIST{i, 4};
             W(end).lims.wrat = -inf;
             W(end).lims.grat = -inf;
             W(end).lims.lrat = -inf;
             W(end).lims.bhp  = 1*atm;
           case 'wrat'
             W(end).lims.wrat = -control.WCONHIST{i, 5};
             W(end).lims.orat = -inf;
             W(end).lims.grat = -inf;
             W(end).lims.lrat = -inf;
             W(end).lims.bhp  = 1*atm;
           case 'grat'
             W(end).lims.grat = -control.WCONHIST{i, 6};
             W(end).lims.orat = -inf;
             W(end).lims.wrat = -inf;
             W(end).lims.lrat = -inf;
             W(end).lims.bhp  = 1*atm;
           case 'lrat'
             W(end).lims.lrat = -sum([control.WCONHIST{i, 4:5}]);
             W(end).lims.orat = -inf;
             W(end).lims.wrat = -inf;
             W(end).lims.grat = -inf;
             W(end).lims.bhp  = 1*atm;
           case 'bhp'
             W(end).lims.bhp  = control.WCONHIST{i, 10};
             W(end).lims.lrat = -inf;
             W(end).lims.orat = -inf;
             W(end).lims.wrat = -inf;
             W(end).lims.grat = -inf;
         end
         W(end).status    = status;
      end
      %W(end).lims  = [];

      %W(end).bhpLimit = 1*atm;
   end
end

%--------------------------------------------------------------------------

function W = process_wconprod(W, control, G, rock, well_id, p, opt)
   for i = 1 : size(control.WCONPROD,1),
      nm = control.WCONPROD{i,1};
      status = strcmp(control.WCONPROD{i,2}, 'OPEN');

      type = lower(control.WCONPROD{i,3});
      switch type,
        case 'orat',
          % 'val' is INJECTION rate.
          val = -control.WCONPROD{i, 4};
          compi = [ 0, 1, 0 ];

        case 'wrat',
          % 'val' is INJECTION rate.
          val = -control.WCONPROD{i, 5};
          compi = [ 1, 0, 0 ];

        case 'grat',
          % 'val' is INJECTION rate.
          val = -control.WCONPROD{i, 6};
          compi = [ 0, 0, 1 ];

        case 'lrat',
          % 'val' is INJECTION rate.
          val   = -control.WCONPROD{i, 7};
          compi = [ 1, 1, 0 ];  % LIQUID rate == water+oil rate (surface)

        case 'resv',
          % 'val' is INJECTION rate.
          val   = -control.WCONPROD{i, 8};
          compi = [ 1, 1, 1 ];

        case 'bhp',
          val   = control.WCONPROD{i, 9};
          compi = [0, 1, 0];  % Doesn't matter.

        otherwise
          dispif(opt.Verbose, ...
                 ['Control mode ''%s'' unsupported for producer ', ...
                  '''%s''.  Ignored.\n'], upper(type), nm);
          continue
      end

      sizeW = numel(W);

      W = buildWell(W, G, rock, control, well_id(nm), p, type, val, ...
                    compi, opt.InnerProduct, -1);

      if numel(W) > sizeW
         W(end).lims.orat = -control.WCONPROD{i, 4};
         W(end).lims.wrat = -control.WCONPROD{i, 5};
         W(end).lims.grat = -control.WCONPROD{i, 6};
         W(end).lims.lrat = -control.WCONPROD{i, 7};
         W(end).lims.bhp  = control.WCONPROD{i, 9};
         W(end).status    = status;
      end
   end
end

%--------------------------------------------------------------------------

function wid = enumerateWells(control)
   [wn, i] = sort(control.WELSPECS(:,1));
   assert (numel(unique(wn)) == numel(wn));

   try
      % Java's O(1) hash table string search support.
      ht = java.util.Hashtable;
      for n = 1 : numel(wn),
         ht.put(wn{n}, i(n));
      end
      wid = @(s) ht.get(s);
   catch %#ok
         % Fall back to (probably) linear structure field name search if Java
         % is unavailable.
      ht = struct();
      for n = 1 : numel(wn),
         ht.(regexprep(wn{n}, '\W', '_')) = i(n);
      end
      wid = @(s) ht.(regexprep(s, '\W', '_'));
   end
end

%--------------------------------------------------------------------------

function control = insertDefaults(control, id, nlayers)
   control.COMPDAT  = insertDefaultCOMPDAT (control, id, nlayers);
   control.WCONINJE = insertDefaultWCONINJE(control);
   control.WCONPROD = insertDefaultWCONPROD(control);
end

%--------------------------------------------------------------------------

function [control, p] = orderCompletions(control)
   comp = [vertcat(control.COMPDAT{:,1}), ...
           (1 : size(control.COMPDAT,1)).'];
   comp = sortrows(comp);

   [w, p] = rlencode(comp(:,1), 1);
   assert (all(diff(w) == 1));

   p = cumsum([1; p]);
   control.COMPDAT = control.COMPDAT(comp(:,2), :);  % Reorder completions.
end

%--------------------------------------------------------------------------

function W = buildWell(W, G, rock, control, i, p, ...
                       type, val, compi, ip, sgn)
   comp  = control.COMPDAT(p(i) : p(i + 1) - 1, :);
   perf  = arrayfun(@(i) active_perf(G, comp(i,:)), ...
                    (1 : size(comp,1)).', 'UniformOutput', false);
   nperf = reshape(cellfun(@numel, perf), [], 1);

   wdir  = rldecode([comp{:,13}], nperf, 2) .';
   Kh    = rldecode([comp{:,10}], nperf, 2) .';
   WI    = rldecode([comp{:, 8}], nperf, 2) .';
   Wdiam = rldecode([comp{:, 9}], nperf, 2) .';
   openShutFlag = rldecode(comp(:, 6), nperf);

   % Zero input value means to calculate quantities ourselves.
   % Set negative whence 'addWell' will DTRT[tm].
   %
   Kh(~ (Kh > 0)) = -1;
   WI(~ (WI > 0)) = -1;

   Wdiam(~ (Wdiam > 0)) = 1*ft;

   assert (all(Wdiam((Kh < 0) & (WI < 0)) > 0), ...
           ['Well bore diameter must be defined unless productivity ', ...
            'index is explicitly provided.']);

   RefDepth = control.WELSPECS{i,5};
   assert (isnumeric(RefDepth));

   if isnan(RefDepth),
      RefDepth = min(G.cells.centroids(vertcat(perf{:}), 3));
   end

   % Remove cells from connection which are shutdown
   perf = vertcat(perf{:});
   [new_perf, ia, ic] = unique(perf, 'last');
   ia = sort(ia); % don't want to sort perf accoring to cell-num
   cstatus = strcmp('OPEN', openShutFlag(ia));
   %ia = ia(ia_open);

   perf  =  perf(ia);
   wdir  =  wdir(ia);
   Kh    =  Kh(ia);
   WI    =  WI(ia);
   Wdiam =  Wdiam(ia);

   if sum(nperf) > 0,
      sizeW = numel(W);
      W = addWell(W, G, rock, perf,                       ...
                  'Type', type, 'Val', val, 'Dir', wdir, 'Kh', Kh,    ...
                  'Sign', sgn, 'Wi', WI, 'Comp_i', compi,             ...
                  'Name', control.WELSPECS{i, 1}, 'Radius', Wdiam./2, ...
                  'RefDepth', RefDepth, 'InnerProduct', ip);
      if numel(W) > sizeW
         W(end).lims = [];
         W(end).cstatus   = cstatus;
      end
   end
end

%--------------------------------------------------------------------------

function compdat = insertDefaultCOMPDAT(control, id, nlayers)
   compdat = control.COMPDAT;
   present = false([size(control.WELSPECS,1), 1]);

   if ~isempty(compdat),
      compdat(:,1) = cellfun(id, compdat(:,1), 'UniformOutput', false);
      present(vertcat(compdat{:,1})) = true;
   end

   % Wells for which no completions have been explicitly set using the
   % 'COMPDAT' input keyword get vertical completions in all layers.
   %
   tmpl = {-1, -1, 1, nlayers, 'OPEN', 0, 0.0, 0.1, -1, 0.0, ...
           'Default', 'Z', -1};
   for w = reshape(find(~present), 1, []),
      compdat = [compdat; [w, tmpl]]; %#ok
   end

   % Fill in defaulted (I,J) locations (of well heel).
   i = vertcat(compdat{:,2}) < 1;
   compdat(i,2) = control.WELSPECS(vertcat(compdat{i,1}), 3);

   i = vertcat(compdat{:,3}) < 1;
   compdat(i,3) = control.WELSPECS(vertcat(compdat{i,1}), 4);
end

%--------------------------------------------------------------------------

function wconinje = insertDefaultWCONINJE(control)
   wconinje = control.WCONINJE;

   if numel(wconinje) > 0,
      % Insert default pressure target (large value) where needed
      i = isnan(vertcat(wconinje{:,7}));
      wconinje(i,7) = { 6895 * barsa };  % \approx 1e5 psia
   end
end

%--------------------------------------------------------------------------

function wconprod = insertDefaultWCONPROD(control)
   wconprod = control.WCONPROD;

   if numel(wconprod) > 0,
      % Insert default pressure value (atmospheric pressure) where needed.
      i = isnan(vertcat(wconprod{:,9}));
      wconprod(i,9) = { 1 * atm };
   end
end

%--------------------------------------------------------------------------

function perf = active_perf(G, comp)
   [ijk{1:3}] = ndgrid(comp{2}, comp{3}, comp{4} : comp{5});
   i    = sub2ind(G.cartDims, ijk{1}(:), ijk{2}(:), ijk{3}(:));
   perf = cart2active(G, i);
end

function W = createDefaultWell(G, rock)
% default well will be handled as injector in checkLims
   rock.perm = ones(G.cells.num,1);
   W = addWell([], G, rock, 1, 'Val', 0, 'Type', 'rate', 'sign', 1, 'Comp_i', 1/3*[1, 1, 1], 'refDepth', G.cells.centroids(1,3));
   W.lims.rate = -inf;
   W.lims.bhp  = inf;

end
