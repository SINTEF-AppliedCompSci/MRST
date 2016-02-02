function [s, krW, krO] = interpolateCorey()
   n      = 3;
   swc    = 0.2;
   sor    = 0.2;
   krOres = 0.5;
   krWres = 0.6;
   n      = 1.5;
   swc    = 0.05;
   sor    = 0.05;
   krOres = 1;
   krWres = 1;

   s = (0:0.01:1)';
   [krW, krO] = corey(s, n, swc, sor, krOres, krWres);

end

function [krW, krO] = corey(s, n, swc, sor, krOres, krWres)
  s = min(max(s, swc), 1 - sor);
  krW = krWres*((s - swc)./(1 - swc - sor)).^n;
  sO = 1 - s;
  krO = krOres*((sO - sor)./(1 - swc - sor)).^n;
end
