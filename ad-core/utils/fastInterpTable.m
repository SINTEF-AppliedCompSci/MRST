function yi = fastInterpTable(X, Y, xi)
   F = griddedInterpolant(X, Y, 'linear', 'linear');
   yi = F(xi);
end
