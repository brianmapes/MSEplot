# MSEplots-pkg package and its parts live here 

Moist static energy plots (with dry static energy and saturation moist statice energy curves), [a better view of soundings](https://github.com/brianmapes/MSEplot/blob/master/StaticEnergyPlots-IWM-VI_extended_abs.pdf) in terms of seeing humidity and the prospects for moist convection. 

The package for PyPI is called MSEplots-pkg. It has MetPy and Siphon as dependencies. 



-------
How to publish to pypi: 

0. cd to package directory on disk corresponding to https://github.com/brianmapes/MSEplot/tree/master/MSEplots-pkg_MetpyDependent

1. Edit setup.py with new release increment. Delete old dist/ contents. Delete egg_info also. 

2. python setup.py sdist

3. twine upload dist/* which was created in step 2.
