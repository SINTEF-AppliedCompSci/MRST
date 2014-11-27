function [W, wells_shut] = updateSwitchedControls(sol, W, varargin)
   if isempty(W) || isempty(sol)
      wells_shut = false;
      return
   end
   
   if isfield(W, 'status')
       active = vertcat(W.status);
   else 
       active = true(size(W));
   end
   W_tmp = W(active);
   
   opt = struct('allowWellSignChange', false, ...
                'allowControlSwitching', true, ...
                'Verbose', mrstVerbose);
   opt = merge_options(opt, varargin{:});

   % Check if producers are becoming injectors and vice versa. The indexes
   % of such wells are stored in inx.
   wsg = vertcat(W_tmp(:).sign);
   ssg = sign(getTotalRate(sol));
   inx = wsg ~= ssg;

   % A well can be set to zero rate without beeing shut down. We update inx
   % to take into account this fact.
   wval = vertcat(W_tmp(:).val);
   wtype = {W_tmp(:).type}';
   inx = inx & ~strcmp('bhp', wtype) & (wval ~= 0);

   inx = find(inx);
   wells_shut = false;
   if ~opt.allowWellSignChange & opt.allowControlSwitching
      if any(inx)
         wells_shut = true;
         ostring = 'Wells ';

         for k = 1:numel(inx)
            W_tmp(inx(k)).status = false;
            ostring = [ostring,  W_tmp(inx(k)).name, ', '];
         end

         if opt.Verbose,
            fprintf([ostring(1:end - 2), ' shut down.\n']);
         end
      end
   end

   % Check if well-controls have been switch, if so, update W
   inx = find(~arrayfun(@(x,y)strcmp(x.type,y.type), W_tmp(:), sol(:)));
   for k = 1:numel(inx)
      fromTp = W_tmp(inx(k)).type;
      toTp   = sol(inx(k)).type;

      if opt.Verbose
         fprintf(['Well ', W_tmp(inx(k)).name,       ...
                  ' has switched from ', fromTp, ...
                  ' to ', toTp, '.\n']);
      end

      W_tmp(inx(k)).type = toTp;
      W_tmp(inx(k)).val  = sol(inx(k)).val;
   end
   W(active) = W_tmp;
end

%--------------------------------------------------------------------------

function qt = getTotalRate(sol)
   typelist = {'qWs', 'qOs', 'qGs'};
   types    = typelist(isfield(sol(1), typelist));
   qt       = zeros([numel(sol), 1]);

   for w = 1:numel(sol)
      for t = reshape(types, 1, []),
         x = sol(w).(t{1});
         if ~isempty(x),
            qt(w) = qt(w) + x;
         end
      end
   end
end
