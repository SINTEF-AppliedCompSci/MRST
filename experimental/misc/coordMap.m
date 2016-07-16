function cmap = coordMap(varargin)
   
   opt.hash_function = [];
   opt = merge_options(opt, varargin{:});
   
   %map = containers.Map;
   map = containers.Map('KeyType', 'double', 'ValueType', 'int32');
   hfun = @default_hash_function;
   if ~isempty(opt.hash_function)
      hfun = opt.hash_function;
   end
   
   count = 0;
   
   function res = impl(coord, varargin)

      key = hfun(coord);
      
      if map.isKey(key)
         res = map(key);
      else
         count = count + 1;
         map(key) = count;
         res = count;
      end
   end
   cmap = @impl;
end

% ----------------------------------------------------------------------------

function hash = default_hash_function(coord)
   
   hash = prod(coord); % worst hash function in the world
   %hash = CalcMD5(coord);
   
end
