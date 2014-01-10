function reginx = getRegMap(val, REGNUM, REGINX, varargin)
   opt = struct('cellInx', []);
   opt = merge_options(opt, varargin{:});
   nt  = numel(REGINX);

   if nt == 1,

      reginx = { ':' };

   else
      if isempty(opt.cellInx),
         % Entire domain.

         if numel(double(val)) ~= numel(REGNUM),
            % Do not use NUMEL in case of ADI.

            error('Region reference for input undefined');
         end

         reginx = REGINX;

      else
         % Reference to (small) subset of all cells

         cellInx = opt.cellInx;
         regnum  = REGNUM(cellInx);

         if numel(cellInx) > 1

            if size(val, 1) ~= numel(cellInx),
               % Do not use NUMEL in case of ADI.

               error('Number of cell indices must be same as input values');
            end

            reginx = arrayfun(@(x) find(x == regnum), 1 : nt, ...
                              'UniformOutput', false);

         elseif numel(cellInx) == 1,
            % Allow single input (for exploring single cell functions).

            reginx = { repmat(regnum, size(val)) };

         else

            error('Got empty cellInx input. This is not happening...');

         end
      end
   end
end
