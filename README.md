# Computed Tomography Imaging Spectrometer simulator
This repository details a computed tomography imaging spectrometer (CTIS) simulator, which is designed to simulate the generation of a center zeroth-order and surrounding first-order diffraction spots of a diffractive optical element (DOE). 
## Overview and use
The CTIS simulator is explained in the supplementary material of our paper _The hybrid approach - Convolutional Neural Networks and Expectation Maximization Algorithm - for Tomographic Reconstruction of Hyperspectral Images_ [[1]](#1).

This repository contains the following MATLAB scripts

- _em.m_:        Expectation Maximization (EM) algorithm (see function for details, inputs and output)
- _generateH.m_: Generate a system matrix $\boldsymbol{H}$ for a computed tomography imaging spectrometer (see function for details, inputs and output)
- _demo.m_:      Demonstration file showing generation of system matrix $\boldsymbol{H}$, simulation of CTIS image $\boldsymbol{g}$ and reconstruction of hyperpsectral cube $\boldsymbol{f}$ using the EM algorithm.

and .mat files:

- cube_HSI_colorchecker.mat: A $400\times200\times216$ hyperspectral cube captured by our pushbroom hyperspectral imaging system.
- halogen25_mat:  A $1\times 25$ column vector containing the spectrum of the halogen lamps used as illumination for the wavelength range 400 nm to 740 nm.
- wavelength25_mat: A $1\times 25$ column vector containing the wavelength axis for 25 spectral bands (400 nm to 740 nm)
- sensitivity25.mat: A $9 \times 25$ matrix containing the diffraction sensitivity (diffraction efficiency of the DOE, transmission of optical system and sensor response) for the 9 diffraction order spots.


## Brief introduction to CTIS
A CTIS system is described by the linear imaging equation:
$$\boldsymbol{g} = \boldsymbol{H}\boldsymbol{f}+ \boldsymbol{n},$$

where $\boldsymbol{g}$ is a column vector with $q^2$ elements and $\boldsymbol{f}$ is the vectorized hyperspectral cube with $r = x \cdot y \cdot z$ voxels, where $x$, $y$ and $z$ denote the two spatial dimensions and the number of spectral channels, respectively, while $\boldsymbol{n}$ corresponds to a random noise vector.
The $q^2 \times r$ system matrix $\boldsymbol{H}$ describes the projection of the $i$-th voxel in $\boldsymbol{f}$ to the $j$-th pixels in $\boldsymbol{g}$ - equivalent to the nine projections shown below, which consist of a central zeroth order, surrounded by eight first orders.

![CTIS_sim_fig1](https://user-images.githubusercontent.com/25078549/159441650-dad683ce-b5ed-4f01-be8a-174402e091c7.png)


## References
<a id="1">[1]</a> 
M. J. Ahlebæk, M. S. Peters, W.-C. Huang, M. T. Frandsen, R. L. Eriksen, and B. Jørgensen, “The hybrid approach - Convolutional Neural Networks and Expectation Maximization Algorithm - for Tomographic Reconstruction of Hyperspectral Images,” 2022, DOI: [10.48550/ARXIV.2205.15772](https://doi.org/10.48550/arXiv.2205.15772).
