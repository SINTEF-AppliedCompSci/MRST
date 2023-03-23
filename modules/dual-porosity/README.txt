====================================================================================================================
MRST's Dual-Porosity Module
====================================================================================================================
This module was developed by Dr. Rafael March and Victoria Spooner, of the Carbonates Reservoir Group of
Heriot-Watt University (https://carbonates.hw.ac.uk/). Any inquiries, sugestions or feedback should be 
sent to r.march@hw.ac.uk, s.geiger@hw.ac.uk or f.doster@hw.ac.uk.

A good starting point is the single phase example in the examples folder.

Description of folders and files:

- ad_models:
--- DualPorosityReservoirModel: Base class for all reservoirs that show dual-porosity behaviour.
--- ThreePhaseBlackOilModelDP:  Three-phase dual-porosity black-oil model.
--- TwoPhaseOilWaterModelDP: Two-phase dual-porosity oil-water model.
--- WaterModelDP: single-phase dual-porosity water model.

- ad_models\equations:
--- equationsOilWaterDP: equations for the TwoPhaseOilWaterModelDP.
--- equationsWaterDP: equations for the WaterModelDP.

- transfer_models:
--- TransferFunction: this is a base class for all the transfer models. All the transfer models should extend this class.
					  The other files ending with "...TransferFunction" are special implementations of transfer functions
					  available in the literature. The most traditional one is KazemiOilWaterGasTransferFunction (see
					  references).
					  
- transfer_models\shape_factor_models:
--- ShapeFactor: this is a base class for all the shape factors. All the shape factors should extend this class.
					  The other files ending with "...ShapeFactor" are special implementations of shape factors
					  available in the literature. The most traditional one is KazemiShapeFactor (see
					  references).					  
		
- examples:
--- example_1ph: oil production by depletion in a dual-porosity reservoir.
---	example_2ph_imbibition: water injection and oil production by imbibition in a dual-porosity reservoir.

References:
[1] Kazemi et al. Numerical Simulation of Water-Oil Flow in Naturally Fractured Reservoirs, SPE Journal, 1976
[2] Quandalle and Sabathier. Typical Features of a Multipurpose Reservoir Simulator, SPE Journal, 1989
[3] Lu et al. General Transfer Functions for Multiphase Flow in Fractured Reservoirs, SPE Journal, 2008

