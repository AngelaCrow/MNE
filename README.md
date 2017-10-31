# ModelosDarwinUICN

This project is develop by Juan M. Barrios <j.m.barrios@gmail.com> and Angela P. 
Cuervo-Robayo <acuervo@gmail.com> to process data for the Darwin Initiative.

The main goal is process species occurrence data to generate Ecological Niche Models,
the `MNE.R` script is the responsible of this task. As an input the 
`MNE.R` script expects a csv file with the occurrences. 

# Tranferring species niche in time
Please use `DataFormating.R`, `Calibration.R`, `Projections.R`, `ENMEvaluation.R` and `RangeShift.R` for transferring models between geographic regions or periods in time.

You must modify:
- `covarAOIDataFolder_fc45`
- `covarAOIDataFolder_fc85` 
- `covarAOIDataFolder_fl45` 
- `covarAOIDataFolder_fl85` 

