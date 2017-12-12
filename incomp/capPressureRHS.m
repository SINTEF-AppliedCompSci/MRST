function cc = capPressureRHS(g, mob, pc, pc_form)
% Compute capillary pressure contribution to system RHS
%
% SYNOPSIS:
%   cc = capPressureRHS(g, mob, pc, pc_form)
%
% DESCRIPTION:
%   Calculate halfface contribution to rhs of darcy equation from the
%   capillary pressure
%
%   g       - Grid data structure.
%
%   mob     - mobilities evaluated in all cells.
%
%   pc      - capillary pressure evaluated in all cells, = fluid.pc(state)
%
%   pc_form - capillary pressure formulation

%{
#COPYRIGHT#
%}

   if any(abs(pc) > 0)
      N = getNeighbourship(g, 'Topological', true);
      [cellno, cf] = getCellNoFaces(g);

      % Make an approximation of pc at faces.
      %
      % This has no effect in tpfa and mpfa (??) but is important for
      % consistency in mimetic discretisations.
      internal = ~ any(N == 0, 2);

      pcf           = zeros(size(N));
      pcf(N > 0)    = pc(N(N > 0));
      pcf           = sum(pcf, 2);
      pcf(internal) = pcf(internal) ./ 2;

      totmob = sum(mob, 2);

      dpc = pc(cellno) - pcf(cf(:,1));

      if strcmpi(pc_form, 'wetting')
         % for using wetting first phase as variable
         cc = (mob(cellno, 2) ./ totmob(cellno)) .* dpc;

      elseif strcmpi(pc_form, 'nonwetting')
         % for using non-wetting first phase as variable for pressure
         % pc = p_nw - p_w
         cc = - (mob(cellno, 1) ./ totmob(cellno)) .* dpc;

      elseif any(strcmpi(pc_form, {'symmetric', 'symetric'}))
         cc = (diff(mob(cellno,:), [], 2) ./ totmob(cellno)) .* dpc;

      else
         error('no such cappilary formulation implemented')
      end

   else
      cc = 0;
   end
end
