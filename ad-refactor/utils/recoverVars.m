function x = recoverVars(eq, n, sol)
   % recover variables x at position n using solutions sol
   solInx = [1:(n-1), (n+1):(numel(sol)+1)];
   x = - eq.jac{n}\(eq.val);
   for k  = 1:numel(solInx)
       if(~isempty(x) && ~isempty(sol{k}))
           x = x - eq.jac{n}\(eq.jac{solInx(k)}*sol{k});
       end
   end
end
