Components for the MuMMI software release. LANL Copyright No. C21095.

The Department of Energy and the National Cancer Institute have developed new software for conducting multi-scale simulations of complex systems. This software, called the Multiscale Machine-Learned Modeling Infrastructure (MuMMI), couples simulations on three spatial scales to study slow, large-scale reorganizations of biomolecular systems with the speed of continuum and coarse-grained models while revealing selected interactions at full atomic precision. In these simulations, coarse-to-fine model conversions are used to spawn relevant fine-scale simulations along chosen order parameters, and fine-to-coarse feedback is used to iteratively improve the accuracy and multi-scale consistency of coarse-scale and continuum simulations. Many parts of the MuMMI framework will be released as open-source software by Lawrence Livermore National Laboratory and the entire MuMMI package is necessary for full functionality. This release covers a subset of the MuMMI components that were developed exclusively at the Los Alamos National Laboratory.

Software components included here have been developed to (a) parameterize cellular signaling proteins for molecular simulation, (b) convert coarse-grained models to all-atom models with multi-scale consistency, (c) quantify biases in machine-learning-guided resampling of existing data, and (d) Sobol resampling for optimal data coverage.

The software will be useful to study systems with interaction potentials that depend on precise arrangements that can only be revealed at high resolution, and whose fundamental dynamics only emerge at large scales and/or long timescales.

We include modified code from: https://github.com/Tsjerk/Backward , which is described in the following publication: https://pubs.acs.org/doi/10.1021/ct400617g We also include a patch that modifies open-source gromacs software version 2019.6, which can be downloaded here: https://manual.gromacs.org/documentation/2019.6/download.html

The entire MuMMI code relies on software beyond that listed here. The code listed here relies on: 
- open-source gromacs software version 2019.6, which can be downloaded here: https://manual.gromacs.org/documentation/2019.6/download.html 
- unmodified MDAnalysis python libraries: https://github.com/MDAnalysis/mdanalysis

