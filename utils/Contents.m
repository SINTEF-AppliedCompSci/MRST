% UTILS
%   Supporting utilities for MATLAB Reservoir Simulation Toolbox (MRST).
%
% Files
%   ADI                         - class: simple implementation of automatic differentiation for easy construction of jacobian matrices.
%   blockDiagIndex              - Compute subscript or linear index to nonzeros of block-diagonal matrix
%   buildmex                    - Wrapper around MEX which abstracts away details of pathname generation.
%   cellDims                    - cellDims -- Compute physical dimensions of all cells in single well
%   cellDimsCG                  - cellDims -- Compute physical dimensions of all cells in single well
%   cellFlux2faceFlux           - Transform cell-based flux field to face-based.
%   computeWellIndex            - Undocumented Utility Function for Connection Transmissibility Factors
%   dinterpq1                   - Compute derivative of piecewise linear interpolant.
%   dinterpTable                - Compute derivative of one-dimensional interpolant, possibly using splines.
%   dispif                      - Produce textual output contingent upon predicate.
%   evalFunStandalone           - Utility for performing stand-alone function evaluation in new matlab session
%   evalFunWrapper              - Utility to launch (mrst-)function evaluation in seperate matlab session 
%   faceFlux2cellFlux           - Transform face-based flux field to cell-based.
%   faceFlux2cellVelocity       - Transform face-based flux field to one constant velocity per cell.
%   findFilesSubfolders         - Find all files in a directory hierarchy
%   formatTimeRange             - Small utility which returns a human readable string from seconds.
%   geomspace                   - Geometrically spaced vector.
%   getSortedCellNodes          - Construct n x 2 table of cell edges with edges oriented the same
%   githubDownload              - Download objects from GitHub (.ZIP or collection of files)
%   incompHydrostaticPressure   - Undocumented Utility Function
%   initVariablesADI            - Initialize a set of automatic differentiation variables
%   interpTable                 - Interpolate a one-dimensional table, possibly using splines.
%   interpTableMEX              - Interpolate a one-dimensional table with MEX acceleration
%   isCoarseGrid                - Check if a grid is a coarse grid or a fine grid
%   mcolon                      - Compute vector of consecutive indices from separate lower/upper bounds.
%   md5sum_fallback             - Alternative implementation of md5sum for systems without C compiler.
%   merge_options               - Override default control options.
%   merge_options_relaxed       - A less general version of merge_options focused on specific choices:
%   moduleGUI                   - Interactive user interface for activation/deactivation of known mrst modules
%   mrstConfigureMETISLocations - Configure Installed Location of METIS Library and Associate Header Files
%   mrstDataDirectory           - Set or retrieve the current canonical data directory for MRST
%   mrstDebug                   - Globally control default settings for MRST debugging information.
%   mrstDefaultMexFlags         - Define Common Compiler and Linker Flags/Libraries for MRST's MEX Functions
%   mrstExamples                - Discover Example M-Files Pertaining to One or More MRST Modules
%   mrstExtraDirs               - Get List of Directories Added to MATLAB's Search PATH by MRST
%   mrstIsLiveEditorDir         - Detect if Script is Run From Live Editor or in Cell Mode
%   mrstModule                  - Query or modify list of activated add-on MRST modules
%   mrstNargInCheck             - Check number of input arguments to function
%   mrstOutputDirectory         - Set or retrieve the current canonical data directory for MRST
%   mrstPath                    - Establish and maintain mapping from module names to system directory paths
%   mrstStartupMessage          - Print a welcome message with helpful commands for new MRST users
%   mrstVerbose                 - Globally control default settings for MRST verbose information.
%   mrstWebSave                 - Get Call-Back for Downloading Online Resources Specified by URLs
%   msgid                       - Construct Error/Warning message ID by prepending function name.
%   multiEig                    - Solve sequence of general (unsymmetric) eigenvalue problems using LAPACK
%   multiSymmEig                - Solve sequence of symmetric eigenvalue problems using LAPACK
%   propertynames               - Undocumented Utility Function
%   reduceToDouble              - Reduce ADI variable to double.
%   require                     - Announce and enforce module dependency.
%   rldecode                    - Decompress run length encoding of array `A` along dimension `dim`.
%   rlencode                    - Compute run length encoding of array A along dimension dim.
%   ROOTDIR                     - Retrieve full path of Toolbox installation directory.
%   subsetMinus                 - Undocumented Utility Function
%   subsetPlus                  - Undocumented Utility Function
%   ternaryAxis                 - Create a ternary axis and mappings to ternary space
%   tetrahedralAxis             - Create a ternary axis and mappings to ternary space
%   ticif                       - Evaluate function TIC if input is true.
%   tocif                       - Evaluate function TOC if input is true.
%   uniqueStable                - Support `unique(A, 'stable')` in all versions of MATLAB
%   value                       - Remove AD state and compact 1 by n cell arrays to matrices

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
