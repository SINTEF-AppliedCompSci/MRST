#ifndef FLUID_H
#define FLUID_H

#include <vector>

#include "configure.h"

/**
 *  Structure for incompressible two-phase fluid.
 */
template<class Real>
class Fluid
{
    public:
        /**
         *  Makes a new fluid object for two phases.
         *  \param
         *      mu_arg      - The viscosities for the two fluids.
         *  \param
         *      rho_arg     - The densities for the two fluids.
         *  \param
         *      sr_arg      - Residual phase saturation for CO2.
         *  \param
         *      sw_arg      - Residual phase saturation for water.
         *  \param
         *      flux_arg    - Whether mobility is computed from table.
         *  \param
         *      kwm_arg     - Phase relative permeability at sr.
         */
        Fluid(Real mu_arg[2], Real rho_arg[2], Real sr_arg, Real sw_arg,
                Real flux_arg, Real kwm_arg[2]);

        /**
         *  Destructor.
         */
        ~Fluid();
    public:
        /// The viscosities for the two fluids.
        Real  mu[2];

        /// The densities for the two fluids.
        Real  rho[2];

        /// Residual phase saturation for CO2.
        Real  sr;

        /// Residual phase saturation for water.
        Real  sw;

        /// Flux interpeter.
        Real  fluxInterp;

        /// Phase relative permeability.
        Real  kwm[2];
};

#include "tmpl/Fluid.tmpl"

#endif
