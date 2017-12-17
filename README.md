# I. Protocol for Running Pareto Scans and Plotting Results

## The following instructions use the noble metal-oxide system as an example (see system 3 under II. Notes on Different Systems)

- Change directories to the MetalOxide folder

`cd MetalOxide`

- To view the contents of this folder, use the 'ls' command (ls for list)

`ls`

- At a minimum, you should see the following files/folders after typing ls:

`DIEL  input.txt  plot_pareto.gnu  plot_spectra.gnu  Scan.c  Scan.exe  SPECTRA`

- Scan.exe is a program that will scan through a number of structural and material parameters, identify the set of Pareto optimal structures, and compute the spectral efficiency, useful power density, and thermal emission spectrum for all Pareto optimal structures

- input.txt is an input file that specifies the material and structural parameters that Scan.exe will use.  Before running Scan.exe, open input.txt and verify it contains parameters relevant to the system you want to study.  You can use the text editor 'vim' to open input.txt

`vim input.txt`

    Nmin_Nmax

    5 28

    d1min_d1max

    0.1  0.3

    d2min_d2max

    0.1  0.3

    vfmin_vfmax

    0.0 1.0

    Tmin_Tmax

    1000 1700

    FilePrefix

    Rh_Alumina

    AbsorberFileName

    DIEL/Rh_Spline.txt


# II. Notes on Different Systems
(1) oxide-oxide layers (ex. Fe3O4-Al2O3, or Fe3O4-Fe2O3) on top of BR

    Fe3O4 is absorbing and has a high melting point.
    Oxidation states has to be well controlled (Fe2O3: transparent).
    There is an intermediate growth condition where both Fe2O3 and Fe3O4 formed.

(2) nitride-nitride layer buried under BR (ex. HfN-AlN)

    AlN, TiN, and HfN are resistive to oxidation when buried under thick oxide layers (hundreds of nanometer). 
    Nitrides may be a better barrier layer for W diffusion outward to oxide BR layers from W substrate.
    However, the effect of coupling to oxide BR is less dramatic when the weak absorber is below BR compared to when it’s above the BR.
    Precise control over stoichiometry at high temperature is needed. (Metal rich: absorbing, N rich: transparent)

(3) noble metal-oxide layers (ex. Ru-Al2O3) on top of BR

    Ru has a high melting point and similar optical properties to W
    RuO2 also has a higher melting point than other noble metal oxides (not sufficient though).

Attached slides contain my preliminary literature surveys on thermal stability and optical properties of the materials mentioned above. FYI, I also attached the refractive indices data obtained from the website: https://refractiveindex.info. So, my main question to you is which materials systems you expect to work best in terms of optical properties. For oxide BR, we agreed to use Al2O3 for low n, and ZrO2 or HfO2 for high n for better thermal stability. I'll measure the refractive indices of these oxides grown by PEALD for next runs of simulations.
