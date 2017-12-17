# I. Protocol for Running Pareto Scans and Plotting Results

## The following instructions use the noble metal-oxide system as an example (see system 3 under II. Notes on Different Systems)

- Change directories to the MetalOxide folder

`cd MetalOxide`

- To view the contents of this folder, use the 'ls' command (ls for list)

`ls`

- At a minimum, you should see the following files/folders after typing ls:

`DIEL  input.txt  plot_pareto.gnu  plot_spectra.gnu  Scan.c  Scan.exe  SPECTRA`

- Scan.exe is a program that will scan through a number of structural and material parameters, identify the set of Pareto optimal structures, and compute the spectral efficiency, useful power density, and thermal emission spectrum for all Pareto optimal structures

- input.txt is an input file that specifies the material and structural parameters that Scan.exe will use.  Before running Scan.exe, open input.txt and verify it contains parameters relevant to the system you want to study.  You can use the text editor 'vim' to open input.txt:

`vim input.txt`

- Once input.txt is open, you will see a variety of parameter labels (words) and parameter values (numbers).  An example of input.txt is provided in Section <b> II. Input File </b> for reference, and the explanation of the labels and values are commented in this readme for clarity.  Comments are in parenthesis and should <b> NOT </b> appear in the actual file input.txt  

## II. Input File
- Note: Comments about the meaning of the input parameters are in parenthesis and should <b> NOT </b> appear in the actual file input.txt

>Nmin_Nmax   (minimum and maximum number of layers in the structure) 
>
>5 28        (scan will consider structures with number of layers between 5 and 28)

>d1min_d1max (minimum and maximum thickness of low refractive index layer) 
>
>0.1  0.3    (scan will consider low refractive index layers between 0.1 and 0.3 microns thick)
 
>d2min_d2max  (minimum and maximum thickness of high refractive index layer) 
>
>0.1  0.3     (scan will consider high refractive index layers between 0.1 and 0.3 microns thick)

>vfmin_vfmax  (minimum and maximum volume fraction of alloy layer)
>
>0.0 1.0      (scan will consider alloys with volume fractions between 0 and 100% of the metal in the oxide)

>Tmin_Tmax    (minimum and maximum temperature of the structures)
>
>1000 1700    (Scan will consider performance of structures when they are between 1000 and 1700 K)

>FilePrefix   (Controls what the name of various output files will be)
>
>Rh_Alumina   (All output files created by Scan.exe will begin with the word 'Rh_Alumina')

>AbsorberFileName (Specifies where Scan.exe should read the data that defines the metal material)
>
>DIEL/Rh_Spline.txt (Scan.exe will read the metal material data from the file Rh_Spline.txt in the DIEL folder... this is data for the metal Rhodium (Rh))


# III. Notes on Different Systems
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
