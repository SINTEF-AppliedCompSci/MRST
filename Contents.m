% Near-wellbore modeling (nwm) module
%
%  -----------------------------------------------------------------------
% |     Copyright (C) 2019 Lin Zhao. See COPYRIGHT.TXT for details.       |
%  ----------------------------------------------------------------------- 
%
% Routines for the implementation and manipulation of near-wellbore 
% modeling (NWM) methodology
%
% AUTHOR:      Lin Zhao (zhaolin_9@126.com)
% INSTITUTION: Reservoir Engineering group, China University of Petroleum
%              (Beijing), Supervisors: Hanqiao Jiang, JunJian Li, and Jie Li 
%              Computational Geosciences group, SINTEF Digital (2018-2019),
%              Supervisor: Knut-Andreas Lie
% ACKNOWLEDGEMENT:
%   Sincerely thanks for the technical supports of Knut-Andreas Lie, 
%   Atgeirr F. Rasmussen, Stein Krogstad, Xavier Raynaud, and Kai Bao from
%   SINTEF Digital; Hanqiao Jiang, JunJian Li, and Jie Li from China 
%   University of Petroleum (Beijing); Xin Sun from China University of 
%   Petroleum (HuaDong); Shuai Ma from China University of Geosciences 
%   (Beijing).
%   Sincerely thanks for the financial support of China Scholarship
%   Council (CSC).
%
%
%
%  COPYRIGHT.txt     - Copyright file
%  license.txt       - GNU License file
%
% /data
%   MSW.data   
%      - Eclipse data file for example 'nearWellBoreModelingMultiSegWell'
%   NWM.data
%      - Eclipse data file for examples 'nearWellBoreModelingGrids' 
%        and 'nearWellBoreModelingWorkFlow'
%   actnum.inc
%      - INCLUDE file of 'NWM.data' for keyword ACTNUM
%   coord.inc
%      - INCLUDE file of 'NWM.data' for keyword COORD
%   perm.inc
%      - INCLUDE file of 'NWM.data' for keyword PERMX
%   zcorn.inc
%      - INCLUDE file of 'NWM.data' for keyword ZCORN
%   trajectory.mat
%      - 3D discrete well-trajectory points for examples 
%        'nearWellBoreModelingGrids' and 'nearWellBoreModelingWorkFlow'
%
% /example
%   nearWellBoreModelingGrids  
%      - Demonstration of the grids in NWM 
%   nearWellBoreModelingMultiSegWell
%      - Demonstration of the coupled model of NWM and multi-segment well
%        (MSW)
%   nearWellBoreModelingSim
%      - Demonstration of the generations of necessary data structures 
%        passed to the mrst AD simulators for the hybrid grid
%   nearWellBoreModelingWorkFlow
%      - Demonstration of the NWM workflow
%
% /example/gridding examples
%   connWRCartGrids 
%      - Example of connecting the well-region grid to Cartesian grid by 
%        the Voronoi (pebi) grid
%   glueRadCartGrids
%      - Demonstration of gluing the 2D radial grid to Cartesian grid
%   layeredGrids
%      - Examples of creating layered grids by function 'makeLayeredGridNWM'
%   simRadCartGrids
%      - Example of simulation on a radial-Cartesian hybrid grid
%   standalonelRadialGrid
%      - Demonstration of building an individual radial grid
%
% /example/gridding examples/demo functions
%   demoComputeAuxPts
%      - An excerpt from 'generateVOIGridNodes' for demonstration purpose
%   demoHandleVoronoiDiag
%      - An excerpt from 'generateVOIGridNodes' for demonstration purpose
%   demoPlotLine
%      - Plot a 3D line
%   demoPlotPoly
%      - Plot a closed polygon
%
% /gridding
%   assembleGrids 
%      - Assemble multiple grids   
%   buildRadialGrid 
%      - Build the 2D radial grid from point and dimension definitions
%   distmesh_2d_nwm
%      - The modified distmesh2d in 'distmesh'
%   extractBdyNodesCells
%      - Extract the sorted boundary nodes and cells (in counterclockwise) 
%        of a inner continuous region
%   generateHWGridNodes
%      - Generate 3D points of all radial HW grid surfaces and 2D planar 
%        points for horizontal well (HW) grid
%   generateVOIGridNodes
%      - Generate 3D points of all surfaces of volume of interest (VOI)   
%        grid and connectivity list corresponding to the 2D planar points
%   getConnListAndBdyNodeWR2D
%      - Get connectivity list and boundary nodes of 2D well region
%   makeConnListFromMat
%      - Make the connectivity list from the node distribution matrix for a
%        structured grid
%   makeLayeredGridNWM
%      - Extrude 2D grid to layered 3D grid according the topology of 2D 
%        grid and provided surface point set
%   passToDistmesh
%      - Generate parameters passed to 'distmesh_2d_nwm'
%   pointsSingleWellNode
%      - Generate the 2D well region points corresponding to single well
%        node
%   radCartHybridGrid
%      - Build the hybrid grid by gluing the radial grid in the near-well 
%        region to the Cartesian grid elsewhere in the reservoir
%
% /models
%   HorWellRegion
%      - Class for HW region in VOI grid
%   MultiSegWellNWM
%      - Derived class for generating necessary data structures for coupled
%        NWM and MSW model
%   NearWellboreModel
%      - Class for generating necessary data structures passed to the mrst 
%        AD simulators for the hybrid grid 
%   VolumeOfInterest
%      - Class for VOI in the Corner-point grid (CPG) or Cartesian grid
%
% /trans
%   computeRadTransFactor
%      - Compute the radial half transmissibility factor (ft) for the 2D 
%        radial grid
%   handleMatchingFaces
%      - Compute intersection relation between layered boundaries of 
%        subgrids
%   handleNonMatchingFaces
%      - Compute the intersection relation of a surface shared two grids
%   
% /utils
%      - Contains some utility functions
