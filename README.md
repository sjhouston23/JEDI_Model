# JEDI_Model

This directory has all of the programs used to produce X-ray spectra from oxygen and sulfur precipitation based on synthetic spectra from Hui et al., 2019.

To run this program, you use the following command:
gfortran -o jedi.x JEDIModel.f08 JEDIOxyPrecip.f08 JEDISulPrecip.f08 formattingOxy.o formattingSul.o InterpFlux.f08 XRayCX.f08 XRayDE.f08 XRayCXAtm2.f08 XRayDEAtm2.f08 Opacity.f08 OpacityAtm2.f08 Absorption.f08 spline.f08 splineinterp.f08 combCXDE.f08 combCXDEAtm2.f08 && ./jedi.x

JEDIModel.f08 is the main program which stitches together all of the other routines. It begins by asking what JEDI ion spectrum you want to open. The file type it expects is .ds2, which is the file type sent to me originally from the JEDI team (i.e. PJ1.d2s, what I have declared as "version" in the code). If you have a new oxygen and/or sulfur file it needs to be placed in ./JunoData/Spectra/Oxy and/or ./JunoData/Spectra/Sul, respectively.

It takes these files, reads them in and interpolates them with the subroutines JEDIOxyPrecip.f08 and JEDISulPrecip.f08 - it also renormalizes them into units of counts/cm^2/s. When doing this, it creates new files that show the interpolated data in ./Output/Oxy/ and ./Output/Sul/. It then takes the interpolated data and multiplies it by the output files from the original ion precipitation code to create results from an ion precipitation with the fluxes seen by that JEDI measurement. These files are located in ./Output/Oxy/[version]/.

If output files have already been created from a previous use of the code, that'll be detected and it will ask if you want to recalculate the files again with a simple [Y]es/[N]o question. The sulfur to oxygen ratio (S:O) can only be calculated in this program if you allow the JEDI...Precip.f08 codes to run by answering with an "N" to the questions.

It then generates multiple spectra for different viewing angles from the photon output files by using the subroutines XRayCX.f08, XRayDE.f08, XRayCXAtm2.f08, XRayDEAtm2.f08. These subroutines calculate the X-ray production from charge exchange (CX) and direct excitation (DE) separately. It also considers the X-rays being generated in a second atmosphere (Atmosphere 2), which is a well-mixed atmosphere. This matters when taking into account opacity effects with the subroutine Opacity.f08 and Absorption.f08.

The X-ray output is placed into ./Output/XRay/ and ./Output/XRayAtm2/. In these directories are directories for different JEDI passes. Going into one of those, you'll find the output from the X-ray subroutines which include the line emission in the directory "Points", synthetic spectra with Gaussian distributions in "Spectra", and finally the line emission for a given photon energy vs. altitude in "SpecVsAlt".

The line emission in the ./Points/ directory is useful for comparison with XMM-Newton data. This is the data I would send to Will Dunn that he would process using the XMM-Newton instrument response functions to generate a comparable spectrum to observation. However, I split everything based on CX and DE. The X-rays themselves have no idea whether they came from CX or DE, so when comparing to observation it doesn't matter whether it came from CX or DE. Thus, the final part of the JEDIModel.f08 code combines the CX and DE output within the ./Points/ directory so it's all in one file.

----------------------------------------------------------------

When doing a run with new JEDI data a few new directories need to be created. Firstly, within both ./Output/Oxy/ and ./Output/Sul/, a directory with the JEDI data name needs to be created. E.g., if you get a new data set titled PJ14.ds2, then go into both of those output directories and create directories titled "PJ14". This will allow the output data to be placed in there. If this isn't done, the code will fail. Next, in both ./Output/XRay/ and ./Output/XRayAtm2/ the same thing needs to be done. However, within the newly created directory, you need to include the directories for Points, Spectra, and SpecVsAlt. Within these, the different viewing angle directories need to be created. To simplify this, I have created a directory template, per se. Located in ./Output/Directory Template/ are the directories titled ./Points, ./Spectra, and ./SpecVsAlt. These have the necessary viewing angle directories within them so you don't have to go creating a bunch everywhere. Instead, just copying and pasting all 3 of these directories into the new ./Output/XRay/PJ14/ and ./Output/XRayAtm2/PJ14/ directories is much easier.

If you want to do a different viewing angle than the default 0°, 60°, and 80°, then you first need to change the variable "Angles" in JEDIModel.f08. This variable is fed into the subroutines, so it only has to be set once. Right now I only allow it to do 3 different angles at a time. If you want to add more, you'll need to go into all the subroutines and account for this. After changing that variable, you need to create a directories in the outputs. Say for example you want to do 55° for PJ14, then you'll need to add the following directories: ./Output/XRay/PJ14/Points/55/ and ./Output/XRayAtm2/PJ14/55_deg/.

Below this are notes on some of the individual subroutines for more information about them.

----------------------------------------------------------------

XRaySpectraCX.f08 is the main program which produces synthetic spectra from charge exchange collisions for various monoenergetic ion beams. It takes line emission points and creates a Gaussian distribution around them, where the variance and sigma can be set.

Opacity.f08 is a subroutine which takes into account the depth effects due to the atmosphere of Jupiter. This requires absorption cross-sections for Hydrogen (Multiplied by 2 for H2), Helium, and Carbon (for CH4) that have been taken from Nataly Ozak's old code. It allows you to account fro multiple viewing angles using the Chapman function. The angle is an input taken from XRaySpectraCX.f08

Absorption.f08 is a subroutine called within Opacity.f08. It requires and input of wavelength (in nanometers) and outputs absorption cross-sections for H, He, and C. It interpolates the original cross-sections to get a value for any wavelength. It uses the spline.f08 and splineinterp.f08 subroutines for the interpolation. H and He are both entirely interpolated with spline, while C has it's tails interpolated with spline, but a very narrow midrange energy linearly interpolated because of a discontinuity of the absorption cross-section due to the n=2 to n=1 principal quantum number change.
