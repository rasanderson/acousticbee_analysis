# acousticbee_analysis

Baseline analyses of acousticbee data from Ben

The **concat_regr** branch contains a simple "regression-style" deep learning model with two branches: 1) meta-data, e.g. temperature, rain, wind, east, north etc. 2) acoustic-data which are 28 indices that describe the 'soundscape' in the colony. The file `data/All_Data_Master_Shortened.csv` contains a separate row for each timestep.

The data are hierarchical in structure, with 3 to 4 colonies at each site. **Note**: the current deep learning model ignores both the hierarchical structure of the data, and its temporal characteristics. So each row (time-point) in the CSV file is treated as an entirely independent data point. This is of course unrealistic, and some sort of temporal model is needed.
