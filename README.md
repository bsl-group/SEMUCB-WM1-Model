# SEMUCB-WM1-Model
- Download SEMUCB_WM1 model and codes to read it and plot it
- There are two codes, depending on the specific purpose:
  1. To obtain values of the model (Vs and xi) at a particular (radius, lon, lat) or to obtain values at a particular radius on an N x N grid - useful for plotting purposes and direct comparison with other models, download: (tar and gzipped files). Note: the values of the model in the crust are not interpretable. See French and Romanowicz (2014) for details of the crustal model.
  2. To predict waveforms, dispersion curves, or travel times of particular phases using SEMUCB_WM1, here is a code that allows you to construct 1D depth profiles of the model at a particular (lat, lon) location: (tar and gzipped files)


__You must use this code for these purposes__. Indeed, it is important to remember that the mantle model needs to be combined with the crustal model that was constructed at the same time. Using our mantle model with another crustal model will not lead to a fair assessment. Our crustal model is not a standard model like crust2.0 or crust1.0, it is a smooth crustal model designed to fit a global dispersion dataset (Shapiro and Ritzwoller, 2002) as well as our waveform dataset, and it is adjusted at each iteration of the inversion. To construct a 3D description of the model on a grid, you should pay close attention to the Moho depths - which are not realistic in the oceans (saturated to 30 km) but are designed to be equivalent to a real crustal models for seismic waves of periods > 30 s.

It is straightforward to construct a 3D grid model from the 1D profiles obtained with the code provided.
