function [W, wells_shut] = updateSwitchedControls(sol, W, varargin)

   if isempty(W) || isempty(sol)
      wells_shut = false;
      return
   end

   opt = struct('allowWellSignChange', false, 'allowControlSwitching', true);
   opt = merge_options(opt, varargin{:});

   % Check if producers are becoming injectors and vice versa. The indexes of
   % such wells are stored in inx.
   wsg = vertcat(W(:).sign);
   ssg = sign(getTotalRate(sol));
   inx = wsg ~= ssg;

   % A well can be set to zero rate without beeing shut down. We update inx to
   % take into account this fact.
   wval = vertcat(W(:).val);
   wtype = {W(:).type}';
   inx = inx & ~strcmp('bhp', wtype) & (wval ~= 0);

   inx = find(inx);
   wells_shut = false;
   if ~opt.allowWellSignChange & opt.allowControlSwitching
      if any(inx)
         wells_shut = true;
         ostring = 'Wells ';

         for k = 1:numel(inx)
            W(inx(k)).status = false;
            ostring = [ostring,  W(inx(k)).name, ', '];
         end

         fprintf([ostring(1:end - 2), ' shut down.\n']);
      end
   end

   % Check if well-controls have been switch, if so, update W
   inx = find(~arrayfun(@(x,y)strcmp(x.type,y.type), W(:), sol(:)));
   for k = 1:numel(inx)
      fromTp = W(inx(k)).type;
      toTp   = sol(inx(k)).type;
      fprintf(['Well ', W(inx(k)).name,       ...
               ' has switched from ', fromTp, ...
               ' to ', toTp, '.\n']);
      W(inx(k)).type = toTp;
      W(inx(k)).val  = sol(inx(k)).val;
   end
end


function qt = getTotalRate(sol)
   typelist = {'qWs', 'qOs', 'qGs'};
   typeindex = isfield(sol(1), typelist);
   types = typelist(typeindex);
   qt = zeros(numel(sol), 1);
   for w = 1:numel(sol)
      for i = 1:numel(types)
         qt(w) = qt(w) + sol(w).(types{i});
      end
   end
end
