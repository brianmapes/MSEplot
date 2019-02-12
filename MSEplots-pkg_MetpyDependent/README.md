### MSEplots project
#### pip install MSEplots-pkg
![Generic badge](https://img.shields.io/badge/Python3.0-<COLOR>.svg)
![Generic badge](https://img.shields.io/badge/MetPy-<COLOR>.svg)
------
A Python package built for the moist static energy (MSE) analysis of sounding data/ model output which provides required vertical profiles of thermodynamic parameters. 

```python
from MSEplots import plots as mpt
:
mpt.msed_plots(pressure,Temp,q,altitude,ent_rate=np.arange(0,2,0.05),entrain=True)
```
<img src="https://github.com/weiming9115/Working-Space/blob/master/MSEplots_metpy/demo.png" width="550" height="400">

1. Required paramters: Air temperature, Mixing ratio, Pressure, Altitude [optional]. NOT specifically for sounding data!
2. Functions are provided for deriving thermodynamic variables eg. potential tmeperature and static energy. All calculations included depend on the metpy.calc.thermo module.
(https://unidata.github.io/MetPy/latest/_modules/metpy/calc/thermo.html)
3. Plotting options: thermo.plots, theta.plots, and mesd_plots
