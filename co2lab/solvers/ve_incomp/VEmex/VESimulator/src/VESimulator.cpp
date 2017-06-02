/*
-------------------------------------------------------------------------------
File Name :

Purpose :

Creation Date : 2010-11-08

Last Modified : Wed Feb  2 15:54:25 2011

Created By :  Martin Ertsaas (martiert@student.matnat.uio.no)
-------------------------------------------------------------------------------
*/
/*! \mainpage CPU vertical average simulator
 *
 * \section intro_sec Introduction
 *
 * This is the introduction.
 *
 * \section install_sec Installation and Compiling
 *
 *     Download the hg close ???
 *     Read the README file in the root directory.
 *     > mkdir build
 *     > cd build
 *     > cmake ..
 *     > make VESimulatorCPU
 *     > make VETransportCPU
 *
 *  This will generate build/mex/VETransportCPU.mexa64 which is the mex interface to the transport code.
 *
 * \subsection Matlab test cases
 *
 * If success full compilation of the code is done. There is example of using VETransportCPU together with MRST
 * in ROOT/matlab directory. Compeared to the standard vertical average matlab module in MRST the difference is that
 * VETransportCPU has to be called with out arguments like
 * > VETransportCPU()
 * if you need to change some parameters, except the state of the simulation or
 * anything in the opts_t struct (See State.h and Options.h). This is to avoid copying data between matlab and the c++ code.
 *
 */

#include "VESimulator.h"
