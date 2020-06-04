function adi = makeadi(vals, jac)

   % three next lines to workaround mex memory issue with sparse arrays
   if isscalar(jac)
      jac = sparse(numel(vals(:)), jac);
   end
   
   adi = ADI(vals(:), jac);
   
end
