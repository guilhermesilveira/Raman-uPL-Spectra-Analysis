# Raman-μPL-Spectra-Analysis
Raman and Micro-Photoluminescence (μPL) spectroscopy are powerful characterization technique providing many useful informations on the sample, they are non-destructives and have a spatial resolution around 100nm, which makes them widely used in nanoscience. This repository contains the source code to analyze raw measurements performed by Raman and μPL spectroscopy.

## Characterized samples
<img align="right" src="https://raw.githubusercontent.com/Aurelien-Pelissier/Raman-uPL-Spectra-Analysis/master/img/nb.png" width=400>

The samples used in our experiments are "nanobeam", it consist of a one dimensional photonic cristal for photons confinement and a quantum well for electrons confinement, both particles can be coupled to create a Quantum Electro Dynamic system, and our goal is to have a better understanding of the physical phenomenons involved in nanobeams by performing space and time dependent measurement.


### Micro-Photoluminescence
<img align="left" src="https://raw.githubusercontent.com/Aurelien-Pelissier/Raman-uPL-Spectra-Analysis/master/img/PL.png" width=300>


μPL measurements provide information about the optical and electronic properties of the material. Particularly, in nanobeam μPL, the background light is from the Quatum well emission and can be fitted with a Gaussian, while the 4 observed peaks refers to the different confinement mode of the photonic cristal and are fitted with a Lorentzian. By anylizing both the position and the quality factor of these curve, we can extract informations about the quality and the temperature of the nanobeam. The code for this sample curve is available in the folder `src/μPL/sample-spectrum`.



&nbsp;

### Raman spectroscopy
<img align="left" src="https://raw.githubusercontent.com/Aurelien-Pelissier/Raman-uPL-Spectra-Analysis/master/img/Raman.png" width=300>
Raman spectroscopy indicate the strain state of the measured sample. The shidt in the spectrum is linear with th deformation.
Three different modes are visibles on the Raman spectra, one corrsepond to silicon and the other two are GaN
spectrum analysis with lorentzian fit. By anylizing both the position and the quality factor of these curve, we can extract informations about the strain state of the nanobeam. The code for this sample curve can be found in the folder `src/Raman/sample-spectrum`.


&nbsp;


&nbsp;


&nbsp;


&nbsp;


&nbsp;


## Measurements & Results

#### Time dependant μPL

We have performed time dependent μPL analysis while changing the potential to see how the nanobeam isa affected by a change in the external electric field.
Spectrum analysis with lorentzian and gaussian fit 
time evolution of QW and PhC mode

#### Raman mapping
Raman maesurement has been performed at different position on the nanobeam. Each spectra are then fited with 3 lorentzian to extract the position of the modes.
mapping of raman shift
spectrum analysis with lorentzian fit


## Report
For more details about the methods used in our measurements, or more detailed results, you can check the report, available in the `report/` folder.
