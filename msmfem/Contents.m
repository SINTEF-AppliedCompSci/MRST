% Routines supporting the multiscale mixed FE method for the pressure equation.
%
% Files
%   basisMatrixHybrid            - Form hybrid versions of the multiscale basis function matrices.
%   basisMatrixMixed             - Form mixed versions of the multiscale basis function matrices.
%   compressPartition            - Renumber coarse block partitioning to remove any empty coarse blocks.
%   dynamicCoarseWeight          - Compute synthetic multiscale weighting function.
%   generateCoarseGrid           - Build coarse grid data structure from partion of existing fine grid.
%   generateCoarseSystem         - Construct coarse system component matrices from fine grid model.
%   generateCoarseWellSystem     - Construct coarse system component matrices for well contributions.
%   partitionCartGrid            - Partition a Cartesian grid.
%   partitionLayers              - Partition grid uniformly in logical (I,J) direction, non-uniformly in K.
%   partitionUI                  - Partition grid uniformly in logical space.
%   processPartition             - Split disconnected coarse blocks into new blocks.
%   solveCoarsePsysBO            - Solve coarsened fine-scale well system (for Black Oil).
%   solveIncompFlowMS            - Solve coarse (multiscale) pressure system.
%   solveIncompFlowMSSpeedUp     - Solve (multiscale) pressure system assuming that some values are precomputed.
%   speedUpMS                    - preassemble basis
%   subFaces                     - Extract fine-grid faces constituting individual coarse grid faces.
%   unpackWellSystemComponentsMS - Extract coarse linear system components from wells.
