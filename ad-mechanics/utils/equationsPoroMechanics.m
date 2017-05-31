function eqs = equationsPoroMechanics(x, fluidp,  G, rock, operators)

   s = operators.mech;
   alpha = rock.alpha;

   eqs{1} = s.A * x - s.gradP * (alpha .* fluidp) - s.rhs;

   fac =  1 / (1e6 * mean(G.cells.volumes));
   eqs{1} = eqs{1} * fac;

end