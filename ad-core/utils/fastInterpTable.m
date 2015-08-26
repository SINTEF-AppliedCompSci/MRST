function yi = fastInterpTable(X, Y, xi)
   if isempty(xi)
       yi = [];
       return
   end
   F = griddedInterpolant(X, Y, 'linear', 'linear');
   yi = F(xi);
end
