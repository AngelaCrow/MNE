# ModelosDarwinUICN

This project is develop by Juan M. Barrios <j.m.barrios@gmail.com> and Angela P. 
Cuervo-Robayo <acuervo@gmail.com> to process data for the Darwin Initiative.

The main goal is process species occurrence data to generate Specie Distribution
Models, the `MNE.R` script is the responsible of this task. As an input the 
`MNE.R` script expects a csv file with the occurrences. 
`DataFormating.R` will arranged ocurrence and enviromental data, for `MNE.R`. 
There are some aditional parameters that must be edited to use this code with your data. In `MNE.R` edit
the following routes:

- `baseDataFolder` where is the base 
