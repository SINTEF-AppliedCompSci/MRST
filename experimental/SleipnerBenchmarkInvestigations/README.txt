This directory investigates the Sleipner Benchmark model.

Certain files are required to create the Sleipner model grids.

'inspectSleipnerGridModels.m' is a script intended to make comparisons
between the models.

'studySleipnerBenchmarkFUN.m' is the main function used to run an injection
scenario. Grids are generated using the function 'makeSleipnerModelGrid.m'.
Once a particular grid is generated, it is written to a .mat file to avoid
re-generation each time grid is desired for injection scenario.

Simulation results saved to .mat file.

Call 'studySleipnerBenchmarkFUN.m' individually, or run a script such as
'variousInjectScenarios.m' to study several different injection scenarios.

'testSleipnerSensFUN.m' and 'testSleipnerSensitivities.m' are used to test
how sensitive the match between observed and simulated CO2 migration is to
a model's caprock topography.