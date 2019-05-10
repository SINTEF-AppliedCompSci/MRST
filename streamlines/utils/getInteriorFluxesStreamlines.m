function flux = getInteriorFluxesStreamlines(G, state, pvol, reverse)
   if nargin < 3
       reverse = false;
   end
   if size(state.flux, 2) > 1
      state.flux = sum(state.flux, 2);
   end

   if reverse
      state.flux = -state.flux;
   end
   % Make array face fluxes for each cell in grid (Not outer).
   cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
   cf     = G.cells.faces;
   flux   = accumarray([cellNo, cf(:,2)], state.flux(cf(:,1)));
   flux   = bsxfun(@rdivide, flux, pvol);
end