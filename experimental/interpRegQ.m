function [yi, dyidxi] = interpRegQ(T, xi, reginx)
   
   if isempty(xi)
      yi     = [];
      dyidxi = [];
      return;
   end
   
   nreg = numel(reginx);
   
   if nreg > 0 && ischar(reginx{1}) && strcmp(reginx{1}, ':'),
      % Special case denoting entire domain in single region.
      
      yi = mex_nakeinterp1(T{1}(:,1), T{1}(:,2));
   
   elseif nreg > 0
      
      yi = zeros(size(xi));

      for k = 1:nreg,
         if ~isempty(reginx{k})
            yi(reginx{k}) = mex_nakeinterp1(T{k}(:,1), T{k}(:,2), ...
                xi(reginx{k}));
         end
      end
   
   end
   
   if nargout > 1,

      dyidxi = dinterpReg(T, xi, reginx);

   end
end
