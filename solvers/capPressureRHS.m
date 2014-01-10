function cc = capPressureRHS(g, mob, pc, pc_form)
% Compute capillary pressure contribution to system RHS
%
%   SYNOPSIS:
%   cc = capPressureRHS(g, mob, pc, gpc, pc_form)
%
%  DESCRIPTION:
%   Calculate halfface contribution to rhs of darcy equation from the
%   capillary pressure
%
%
%   g       - Grid data structure.
%
%   mob     - mobilities evaluated in all cells.
%
%   pc      - cappilary pressure evaluated in all cells, = fluid.pc(state)
%
%   pc_form - capillary pressure formulation

   if ~all(pc==0),
      dim = g.griddim;%size(g.nodes.coords,2);
      %assert (1 < dim && dim < 4);
      cellno = rldecode(1 : g.cells.num, diff(g.cells.facePos), 2) .';
      % make an approximation of pc at faces.
      % This has no effect in tpfa and mpfa (??) but
      % is important for consitancy in mimetic
      pcf=zeros(g.faces.num,2);
      pcf(g.faces.neighbors>0)= pc(g.faces.neighbors(g.faces.neighbors>0));
      internal = sum(g.faces.neighbors>0,2)==2;
      pcf=sum(pcf,2);
      pcf(internal)=pcf(internal)/2;
      totmob=sum(mob,2);
      if(strcmp(pc_form,'wetting'))
         % for using wetting first phase as variable
         cc     = (mob(cellno,2)./totmob(cellno)) ...
                  .*(pc(cellno)-pcf(g.cells.faces(:,1)));
      elseif(strcmp(pc_form,'nonwetting'))
         % for using non-wetting first phase as variable for pressure
         % pc = p_nw -p_w
         cc     = -(mob(cellno,1)./totmob(cellno)) ...
                  .*(pc(cellno)-pcf(g.cells.faces(:,1)));
      elseif(strcmp(pc_form,'symetric'))
         cc     = ((mob(cellno,2)-mob(cellno,1))./totmob(cellno)) ...
                  .*(pc(cellno)-pcf(g.cells.faces(:,1)));
      else
         error('no such cappilary formulation implemented')
      end
   else
      cc = [];
   end
end


