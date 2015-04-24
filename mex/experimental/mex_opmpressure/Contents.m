% Files
%   mex_cfs_tpfa.m               - Perform single fixed point iteration using compiled C code.
%   mex_compute_coarse_BI.m      - Compute coarse-scale inverse inner product from fine-scale contributions
%   mex_compute_coarse_contrib.m - Generate contributions to coarse-scale linsys using compiled C code.
%   mex_compute_fs_flux.m        - Derive fine-scale contact fluxes from basis-function projection
%   mex_compute_press_flux.m     - Derive pressure and flux from hybrid system using compiled C code.
%   mex_compute_trans.m          - Compute one-sided transmissibilities using compiled C code.
%   mex_generate_coarsegrid.m    - Build coarse grid data structure from partition of existing grid (MEX)
%   mex_generate_coarsesystem.m  - Build coarse system data structure from partition of existing grid (MEX)
%   mex_ifsh.m                   - Discretise and solve flow equation using compiled C code.
%   mex_ifsh_ms.m                - Discretise and solve flow equation using compiled C code.
%   mex_ip_simple.m              - Compute 'ip_simple' inner product values using compiled C code.
%   mex_partition_compress.m     - Remove empty blocks/bins using compiled C code.
%   mex_partition_invert.m       - Invert cell-to-block map (creating block-to-cell) using compiled C code.
%   mex_partition_process.m      - Split disconnected blocks into new blocks using compiled C code.
%   mex_partition_ui.m           - Partition grid uniformly in logical space using compiled C code.
%   mex_schur_comp.m             - Compute hybrid system component matrices using compiled C code.
%   mex_schur_comp_symm.m        - Compute hybrid system component matrices using compiled C code.
%   qfs_mex.m                    - {
%   test_mex_ifsh.m              - {
%   test_mex_schur_comp_symm.m   - run ../../startup

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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
