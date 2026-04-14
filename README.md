# Actively updating and modifying
# Version 
Original Version Parameters:
- Multigap Configuration (total thickness 0.81 cm)
- Uses Gas Table Method
- 15 kV
- 4 ns time window
- 200 timebin
- 0.1 AvalancheMicroscopic Time Window
- Pion (100 GeV, this value is just to follow the value of the Garfield Website Guide.)

1. `Original` - From `Example/RPC`, changed to positive signal and increase pion energy to 100 GeV
2. `O_NoGasTable` - Change to no gas table method, comparison with `Original` 
3. `O_1Gap` - Change from `Original` multigap configuration to single gap.

Above are all versions with only signal plots. Starting from below version, an extra Drift Line Plot is included. 

1. `O_ViewDrift` -  From `Original`, added view drift line feature.
2. `O_ViewDrift_Vneg_1Gap`- From `O_1Gap`, added drift line plot, and flip the voltage by a negative sign.
3. 
