function [CG dual] = createPermutationMatrix(A, dual, CG, varargin)
% Create permutation matrices for the coarse orderings in the MsFV
%
% SYNOPSIS:
%    Creates permutation matrices
% DESCRIPTION:
%    The MsFV is dependent on reordering the system matrices. This takes a
%    dual grid and coarse grid and uses them to create permutation matrices
%    for both primal and inner, edge... orderings
%
% REQUIRED PARAMETERS:
%    A    -- Not in use
%    dual -- Dual coarse grid such as defined by dualPartition
%    CG   -- Coarse grid as defined by for example partitionUI
%
% OPTIONAL PARAMETERS:
%    wells   -- well data structure
%    SpeedUp -- Alternate formulation of the method requires alternate ordering matrices
%
% RETURNS:
%    CG   -- Coarse grid with updated partition structure to account for
%    wells
%    dual -- Updated dual grid with permutation matrices
% NOTE:
%    This function is meant for encapsulation within the MSFVM solver, and
%    is subject to change

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('wells', [], 'Speedup', false);
    opt = merge_options(opt, varargin{:});

    [N M] = size(A);
    assert(N==M,'Matrix should be NxN, not NxM');

    if opt.Speedup
        dual.ll = setdiff(dual.lineedge, dual.nn);
        dual.ee = setdiff(dual.ee, dual.ll);
    else
        dual.ll = [];
    end

    % Hack for unstructured grid test where some coarse cells are skipped
    % from categorization because they do not really have enough faces to
    % be 3d cells
    dual.nn = dual.nn(dual.nn~=0);

    % Normalize to row vectors...
    dual.nn = dual.nn(:)';
    dual.ee = dual.ee(:)';
    dual.ii = dual.ii(:)';
    dual.ll = dual.ll(:)';

    Nbhp = 0;
    Nw = numel(opt.wells);
    for i = 1:Nw
        w = opt.wells(i);
        if strcmp(w.type,'bhp')
            Nbhp = Nbhp + 1;
        end
    end

    % "Real" nodes are all nodes not corresponding to BHP well cells
    Nc = N - Nw;
    real_nodes = N - Nbhp;
    passed = false;
    for i = 1:Nw
        w = opt.wells(i);
        if strcmp(w.type,'rate')
            if passed
                error('Rate wells should come before BHP wells in the msfvm solver');
            end
            contact_nn = intersect(w.cells, dual.nn);
            contact_ee = intersect(w.cells, dual.ee);
            contact_ss = intersect(w.cells, dual.ll);
            if any(contact_nn)
                dual.nn = [dual.nn Nc+i];
            elseif any(contact_ss)
                dual.ll = [dual.ll Nc+i];
            elseif any(contact_ee)
                dual.ee = [dual.ee Nc+i];
            else
                dual.ii = [dual.ii Nc+i];
            end

            %add the well nodes to the same partition as the connected
            %cells (the first cell is chosen arbitrarily
            CG.partition = vertcat(CG.partition, CG.partition(w.cells(1)));
        else
            passed = true;
        end
    end
    if size(dual.ii,1) < size(dual.ii,2)
        ordering = double([dual.ii dual.ee dual.ll dual.nn]);
    else
        ordering = double([dual.ii; dual.ee; dual.ll; dual.nn]);
    end


    Nc = length(ordering);
    %store the actual number of cells
    dual.N = Nc;

    dual.P = sparse(1:Nc,...
                    ordering,...
                    1, Nc, Nc) > 0;
   %% Flux permutation matrix
   ordering = zeros(real_nodes,1);

   ind = 1;
   for i=1:CG.cells.num
       tmp = find(CG.partition == i);
       ordering(ind:(ind+length(tmp)-1)) = tmp;
       ind = ind + length(tmp);
   end

   dual.P_flux = sparse(1:real_nodes,...
                    ordering(1:real_nodes),...
                    1, Nc, Nc) > 0;
return
