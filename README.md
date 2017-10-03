# ModelosDarwinUICN

This project is develop by Juan M. Barrios <j.m.barrios@gmail.com> and Angela P. 
Cuervo-Robayo <acuervo@gmail.com> to process data for the Darwin Initiative.

The main goal is process species occurrence data to generate Specie Distribution
Models, the `MNE.R` script is the responsible of this task. As an input the 
`MNE.R` script expects a csv file with the occurrences. 

There are some aditional parameters that must be edited to use this code with your data. 
In `DataFormating.R` edit the following routes:
- `baseDataFolder` where is the spatial data:
  `shapePath`
  `shapeLayer` 
  `covarDataFolder`

# Tranferring species niche
If you are modeling species niche under diferente conditions, you must modify:
- `covarAOIDataFolder_fc45`
- `covarAOIDataFolder_fc85` 
- `covarAOIDataFolder_fl45` 
- `covarAOIDataFolder_fl85` 

