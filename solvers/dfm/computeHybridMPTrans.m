function [g,S] = computeHybridMPTrans(g,S)
% Computes the hybrid-hybrid transmissibilities. Also update the grid and
% transmissibility structure accordingly.
%
% The call to computeMultiPointTrans with hybrid = true gives a
% transmissibility matrix where connections between hybrid cells are
% represented as boundaries (they are not connectiod). This function
% modifies the transmissibility matrix to produce connected faces.
% In the case that more than two hybrid cells meets (fracture
% intersection), we furthermore remove the small cell in the interaction.
%
% If the pressure values at the intersections between hybrid cells are denoted u_0,
% the discretized elliptic pressure equation can be splitt such that:
%
%( A B^T ; B D ) (u u0)^T = (f 0)^T,
%
% The pressure at the intersection can then be removed such that
%
% u0 = inv(D)*B
%
% Using the above expression the flux can be expressed%
% q = (T T0) (u u0)^T = (T - T0*inv(D)*B)u
%
%
% SYNOPSIS
%
%   [g, S] = computeHybridMPTrans(g,S)
%
%  PARAMETERS
%
%   g       - Grid data structure.
%   S       - Transmissibility data structure.
%
% Copyright 2011-2012 University of Bergen
%
% This file is licensed under the GNU General Public License v3.0.

% Pick out the transmissibility.
T = S.T;

% Assemble the discretization matrix. Cut and paste from incompMPFAlegacy
nc     = g.cells.num;
cellNo = rldecode(1:g.cells.num, diff(g.cells.facePos), 2) .';
C      = sparse(1:size(g.cells.faces, 1), cellNo, 1);
b      = any(g.faces.neighbors==0, 2);
I1     = [(1:g.cells.num)'; g.cells.num + find(b)];
D      = sparse(1:size(g.cells.faces,1), double(g.cells.faces(:,1)), 1);
A      = [C, -D(:,b)]' * T(:,I1);


% Pick out coefficients
I(b)    = g.cells.num + 1 : g.cells.num + sum(b);
num_hc  = length(g.hybridNeighbors.facePos)-1;
[i_b,j_b,s_b] = find(A(1:nc,I(g.hybridNeighbors.faces)));

% let u_k = u0_j for all hybrid cellfaces k connected to an intersection 0.
% j maps facenumber to intersection number.
j = rldecode(1:num_hc,diff(g.hybridNeighbors.facePos),2)';

% Assemble the hybrid-matrix coefficients.
B = sparse(i_b,j(j_b),s_b);

% pick out and assemble the hybrid-hybrid coefficients
[i_d,j_d,s_d] = find(A(I(g.hybridNeighbors.faces),I(g.hybridNeighbors.faces)));
D = sparse(j(i_d),j(j_d),s_d);

% pick out the hybrid transmissibilities.
I2 = g.hybridNeighbors.faces + nc;
[i_t,j_t,s_t] = find(T(:,I2));
T0 = sparse(i_t,j(j_t),s_t,length(T),num_hc);

% remove the inner boundaries
T(i_t,I2) = 0;

% update the transmissibility matrix
T(:,1:nc) = T(:,1:nc) - T0/D*B';

% remove the inner boundaries?
T(:,I2) = 0;

% update the transmissibility structure.
S.T = T;

% store cell2cell connections in the grid structure.
g.cells.neighbors = double(g.faces.neighbors(g.hybridNeighbors.faces(g.hybridNeighbors.neighbors)));

