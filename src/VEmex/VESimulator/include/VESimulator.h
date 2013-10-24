#ifndef VESIMULATOR_H
#define VESIMULATOR_H

#include <vector>
#include <algorithm>
#include <cstring>

#include <assert.h>
#include <math.h>

#include <omp.h>
#include <cmath>
#include "State.h"
#include "Grid2D.h"
#include "Rock.h"
#include "Fluid.h"
#include "Well.h"
#include "BC.h"
#include "Options.h"
#include <stdio.h>
#include "configure.h"

/**
 *  The main simulator class. Controlls the simulation.
 */
template<class Real>
class VESimulator
{
    public:
        /**
         *  Makes an empty simulator. No structure present.
         *  \todo
         *      Implement!
         */
        VESimulator();

        /**
         *  Makes a simulator with given options.
         *  \param
         *      state   - The current state of the simulator. See State.h
         *  \param
         *      grid    - The 2D grid of the simulator. See Grid2D.h
         *  \param
         *      rock    - The rock structure of the reservoir. See Rock.h
         *  \param
         *      fluid   - The fluid structure for this problem. See Fluid.h
         *  \param
         *      opts    - Optional parameters for the simulator. See Options.h
         */
        VESimulator(State<Real> *state, Grid2D<Real> *grid, Rock<Real> *rock,
                Fluid<Real> *fluid, opt_t *opts);

        /**
         *  Makes a new simulator with grid and rock structure defined by the
         *  files given.
         *  \param
         *      grdecl          - Name of the grid file.
         *  \param
         *      porosity        - Name of the file containing porosity
         *      information.
         *  \param
         *      permeability    - Name of the file containing permeability
         *      information.
         *  \todo
         *      Implement!
         */
        VESimulator(char *grdecl, char *porosity, char *permeability);

        /**
         *  Destroys the simulator.
         */
        ~VESimulator();

        /**
         *  Does one explicit timestep. The state objekt is updated after this
         *  function have terminated to contain the new state of the reservoir.
         *  \param
         *      dT  - The timestep to do.
         *  \note
         *      Not fully implemented.
         */
        void doTimeStep(Real& dT);

        /**
         *  Computes a source array containing contributions from the boundary
         *  conditions, the wells and the src.
         *  \return
         *      The collected contributions from the sources.
         *  \note
         *      Not fully implemented.
         */
        Real* sources();

        /**
         *  Gives a pointer to the current state of the simulation.
         *  \return
         *      The current state of the simulation.
         */
        State<Real>* getState() const;

        /**
         *  Change the current state
         *  \param
         *      state - The new state.
         */
        void changeState(State<Real> *state);

        /**
         *  Change the current options.
         *  \param
         *      opts - The new options.
         */
        void changeOptions(opt_t *opts);
    private:
        State<Real>   *m_state;
        Grid2D<Real>  *m_grid;
        Rock<Real>  *m_rock;
        Fluid<Real> *m_fluid;
        Well<Real>  *m_wells;
        BC<Real>    *m_bc;

        bool    verbose;
        bool    computeDt;
        bool    intVert;
        bool    intVert_poro;
        bool    semi_implicit;
        bool    no_dif;
        bool    central;
        bool    gravity_upwind;
        Real    timestep;
        Real    heightWarn;
        Real    gravity;

        /// GPU memory. Only used on the GPU.
        /*
        float *poro;
        float *kr_H;
        float *dtkr_H;
        float *n;
        float *grav;
        float dfgrav, dfflux;
        */
        void precomputeGPUStuff();
        Real computeTimeStep(const Real& dT,
                             const Real* dflux,
                             const Real* q,
			     const Real& rho_diff,
                             const std::vector<Real>& poro,
			     const std::vector<Real>& grav,
                             const std::vector<Real>& dtkr_H);
};

#include "tmpl/VESimulator_tmpl.hpp"

#endif // VESIMULATOR_H
