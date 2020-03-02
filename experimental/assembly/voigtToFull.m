function [constit,cfull] = voigtToFull(G, C)
% Convert  stiffness matrix from voigt format (as returned by Enu2C) to full format (as
% used in CC code)
   

   Nd = G.griddim;
   
   if Nd == 1
      nlin = 1;
      Voigt = 1; 
      Full = 1;
   elseif Nd == 2
      nlin = 3;
      warning('hack to get comaring with CC to work for 2D??')
      Voigt = [1 0 0 0; 0 0 0 1; 0 1 1 0];% probably wrong
      %Voigt =[ 1 0 0; 0 0 1;0 1 0;0 0 1];
      %Voigt = [1 0 0 0; 0 0 0 1; 0 .5 .5 0];
      Full = [1 0 0; 0 0 1; 0 0 1; 0 1 0];
   elseif Nd ==3
      nlin = 6;
      Voigt = [1 0 0 0 0 0 0 0 0; 0 0 0 0 1 0 0 0 0; 0 0 0 0 0 0 0 0 1; 0 0 0 0 0 .5 0 ...
               .5 0; 0 0 .5 0 0 0 .5 0 0; 0 .5 0 .5 0 0 0 0 0];
      Full = [1 0 0 0 0 0; 0 0 0 0 0 1; 0 0 0 0 1 0; ...
              0 0 0 0 0 1; 0 1 0 0 0 0; 0 0 0 1 0 0; ...
              0 0 0 0 1 0; 0 0 0 1 0 0; 0 0 1 0 0 0];
   end

   nc = G.cells.num;
   constit = cell(nc, 1);

   %[i,j] = blockDiagIndex(repmat(nlin,G.cellnum,1),repmat(nlin,G.cellnum,1))
   dims = G.griddim;
   cfull=nan(G.cells.num,(dims*dims).^2);
   for ic = 1 : G.cells.num
      c = reshape(C(ic, :), nlin, nlin);
      constit{ic} = Full*c*Voigt;
      cfull(ic,:)=constit{ic}(:);
   end
   
end
