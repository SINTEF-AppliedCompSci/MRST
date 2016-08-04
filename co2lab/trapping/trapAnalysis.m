function res = trapAnalysis(Gt, method, varargin)
% Compute and summarize the relevant trap analysis information
% SYNOPSIS:
%   function res = trapAnalysis(Gt, method)
%
% PARAMETERS:
%   Gt     - top surface grid to analyze
%   method - true : use cell-centroid-based implementation
%            false: use edge-based implementation
% OPTIONAL PARAMETERS:
%   varargin - pair of 'key'/'value', where currently supported key is:
%              * closed_boundary_edges - by default, boundary edges are
%                                       considered open, but the boundary
%                                       edges whose indices are found in the
%                                       vector 'closed_boundary_edges' will
%                                       be considered closed (no flow across)
%     
% DESCRIPTION:
%   The function computes and summarizes information that describes the
%   trapping structure using either a cell-centroid-based implementation
%   that requires the Matlab Boost Graph Library (from 3rd-party module
%   'matlab_bgl') or an edge-based implementation.  Regardless of
%   implementation, the resulting information pertains to cells (not
%   edges). 
%
% RETURNS:
%   res - structure describing the trap structure of Gt with the following
%         fields:
%         - traps        - one index per cell in Gt.  Gives index of trap
%                          cell belongs to, or '0' if the cell does not
%                          belong to a trap.
%         - trap_z       - One value per trap, giving the z-value at its
%                          bottom (spill point)
%         - trap_regions - one index per cell in Gt.  Gives the index of
%                          the trap that the cell 'spills into', or '0' if
%                          the cell spills out of the domain.  NB: A cell
%                          does not have to lie in a trap in order to spill
%                          into it. 
%         - trap_adj     - adjacency matrix describing connection between
%                          traps.  Values are thus 0 or 1).  A nonzero
%                          value at (i,j) indicates that region 'i' spills
%                          directly into region 'j'.  This matrix is stored
%                          on the sparse format.
%         - cell_lines   - One cell array per trap, conaining the 'rivers'
%                          exiting that trap.  A river is presented as a
%                          sequence of consecutive grid cells that lie
%                          geographically along the river.  A river starts
%                          in a trap and ends either in another trap or at
%                          the boundary of the domain.
%         - top          - indices for all cells that represents local
%                          maxima in the grid.  (NB: these are all trap
%                          cells, but there may be more than one local
%                          maxima per trap)
%
% EXAMPLE:
%
%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

assert(isscalar(method),...
      ['Input parameter <method> must be a logical. '...
      'Use true/false (or 1/0), not ''true''/''false''.'])
  
  opt.closed_boundary_edges = [];
  opt = merge_options(opt, varargin{:});
  
  if method
      % we will use the cell-based method
      if ~isempty(opt.closed_boundary_edges)
         warning(['Cell-based algorithm does not currently support closed ' ...
                  'boundary edges.  Use edge-based algorithm instead if needed. ' ...
                  ' Currently proceeding without closed boundaries.']);
      end
      mlist = mrstModule();
      moduleCheck('matlab_bgl','coarsegrid');
      res = cell_based_trap_analysis(Gt);
      mrstModule('reset', mlist{:})
  else
      % we will use the edge-based method
      res = edge_based_trap_analysis(Gt, opt.closed_boundary_edges);
  end
end

%===============================================================================
function res = edge_based_trap_analysis(Gt, closed_bedges)
   
    % Identifying closed boundary nodes, if any
    closed_bnodes = Gt.faces.nodes(mcolon(Gt.faces.nodePos(closed_bedges), ...
                                          Gt.faces.nodePos(closed_bedges+1)-1));

    % Computing essential trap information, based on geometry of edges
    ntraps = computeNodeTraps(Gt, closed_bnodes);

    % Projecting trap information onto cells (from edges)
    [ctraps, ctrap_zvals, ctrap_regions, csommets, cadj, crivers] = ...
        n2cTraps(Gt, ntraps.trap_regions, ntraps.trap_zvals, ntraps.dstr_neigh, ...
                 ntraps.connectivity, ntraps.traps, ntraps.rivers);  
    
    res.traps        = ctraps;
    res.trap_z       = ctrap_zvals;
    res.trap_regions = ctrap_regions;
    res.trap_adj     = cadj;
    res.cell_lines   = crivers';
    res.top          = csommets;
end

%===============================================================================
function res = cell_based_trap_analysis(Gt)
    trap_st = findTrappingStructure(Gt);
    res.top = trap_st.top;
    
    conn_st = findTrapConnections(trap_st.Gtop, trap_st.z_spill_loc);
    res.traps    = conn_st.traps;
    res.trap_adj = conn_st.trap_matrix;

    num_traps = max(res.traps);
    res.trap_z = zeros(num_traps,1);
    for i = 1:num_traps
        trapcells=find(res.traps == i);
        res.trap_z(i) = trap_st.z_spill_loc(trapcells(1)); % should be the
                                                           % same for all
                                                           % trap cells
        if(~isempty(conn_st.cell_line_traps))                                                
            trap_rivers = find(conn_st.cell_line_traps(:,1) == i);
            for j = 1:numel(trap_rivers)
                res.cell_lines{i}{j} = conn_st.cell_lines{trap_rivers(j)};
            end
        end
    end
    % remove spillpoint from trap definition
    %%{
    for i=1:numel(res.trap_z)
       tcells=find(res.traps==i);
       bcells=tcells((Gt.cells.z(tcells)==res.trap_z(i)));
       if(numel(bcells)>1)
          dispif(mrstVerbose, 'Warning, Trap boundary or trap.') 
       end
       z=max(Gt.cells.z(tcells));
       assert(z<res.trap_z(i));
       % remove b cell from trap
       res.traps(bcells)=0;       
    end
    %}
    % Computing trap regions (based on separate algorithm for now)
    %[C, N, CC] = maxTPFAGravityMatrix(Gt);
    C = maxTPFAGravityMatrix(Gt);
    b = components(C+C');
    % b = b(1:end-1); % removing 'out-of-region' index
    res.trap_regions = zeros(size(res.traps));
    for c = unique(b)'
        affected_cells = find(b == c);
        affected_region = nonzeros(unique(res.traps(affected_cells)));
        
        myswitch=numel(affected_region);
        if(myswitch>0)
            tcells=affected_cells(res.traps(affected_cells)>0);
            if(~(min(Gt.cells.z(tcells))==min(Gt.cells.z(affected_cells))))
               myswitch=0;
               dispif(mrstVerbose, ['Warning, top point of region higher than ' ...
                                   'trap in region.  Include in outside traps.\n']); 
            end            
        end
        switch myswitch;
          case 0
            % this is the 'outside' region
            res.trap_regions(affected_cells) = 0;
          case 1
            % this should be the normal case            
            res.trap_regions(affected_cells) = affected_region(1);
          otherwise
            dispif(mrstVerbose, ['Warning, ambiguous spill region detected.  ' ...
                                'Assign to highest trap\n'], numel(affected_region)); 
            [m,i]=min(Gt.cells.z(tcells));                  %#ok<ASGLU>
            res.trap_regions(affected_cells) = res.traps(tcells(i(1)));
        end
    end
end

