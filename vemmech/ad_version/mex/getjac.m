function val = getjac(adivar)
   val = adivar.jac;
   
   for i = 1:numel(val)
      if ~issparse(val{i})
         val{i} = sparse(val{i});
      end
   end
   
end
