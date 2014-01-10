% NUC Support for non-uniform coarsening/agglomeration method
%
% Files
%   mergeBlocks      - Merge blocks in a partitioning that are smaller than
%                      the given limit.
%   mergeBlocks2     - Alternative implementation of merge primitive that
%                      seeks to uphold the upper bound on the flow
%                      indicator and optionally on the number of cells in a
%                      block.
%   refineBlocks     - Refine blocks for which indicator value exceeds
%                      a given upper limit using either a uniform or a
%                      greedy refinement algorithm.
%   refineGreedy     - Refine blocks in a partition using a greedy
%                      algorithm that grows a new block inward from the
%                      block boundary, adding one ring of neighbors in each
%                      iteration.
%   refineGreedy2    - Alternative greedy algorithm that may grow only
%                      parts of the neighboring ring to ensure that upper
%                      bound is not violated.
%   refineGreedy3    - Same as refineGreedy2, except that cells are sorted
%                      by increasing numbers of faces shared with cells in
%                      the block
%   refineGreedy4    - Same as refineGreedy2, except that cells are sorted
%                      by difference in flow indicator from that of the
%                      expanding block.
%   refineUniform    - Refine blocks in a partition by uniform partitioning
%   segmentIndicator - Segments a fine grid into blocks according to
%                      indicator.
%
% Tutorials (in subdirectory "examples")
%   example1         - Discussion of different flow indicators combined
%                      with uniform refinement of high-flow zones or the
%                      NUC algorithm.
%   example2         - Examples of hybrid grids that combined flow adaption
%                      with a regular partition.
%   example3         - Examples of contrained coarsening.
%
%   adaptiveRefinement - Discusses grids that dynamically adapt to an
%                      advancing saturation front.
%   runSPE10_NUC     - The NUC algorithm applied to individual layers of
%                      the SPE10 dataset.
%   simpleNUC        - Discusses the four basic steps of the NUC algorithm.
