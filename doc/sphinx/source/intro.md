# Introduction

Seismic wave propagation stands at the center of numerous endeavours in science and engineering, ranging from deep-earth tomography to seismic hazard assessment, nuclear monitoring, resource exploration, dynamic earthquake rupture and tsunami generation, volcano monitoring, earthquake early warning, non-destructive testing and extraterrestrial geophysics. Solving the underlying equations accurately in realistic 3D media at resolutions where observed data exists is a formidable challenge, ranking amongst the most demanding supercomputing problems. 

AxiSEM3D is a versatile solver for wave propagation in elastic, viscoelastic, acoustic, fully anisotropic media in whole planets, spherical sections, or local-scale Cartesian domains. Its characteristic feature is a flexible and automated adaptivity to the complexity of underlying structures and wavefields, such that it can run up to 5 orders of magnitude faster in 1D or axisymmetric structures, and up to 2-3 orders of magnitude faster in smooth 3D structures compared to fully discretised 3D models with conventional methods (such as SEM). We call this Azimuthal Complexity Adaptation (ACA), such that AxiSEM3D is an ACA-SEM. The code simulates singular point sources (e.g. single forces, moment tensors, explosions) inside elastic and fluid media (including the atmosphere), but can also be used for simultaneous multiple sources. It boasts an extensive array of output options, including singular seismograms, surfaces, wavefields in 2D and 3D, each for a range of quantities such as displacement, acceleration, pressure, curl, stress tensor, deformation. 

It has extra functionality for a wavefield injection approach [Leng et al., 2020] to further boost computational efficiency for localised heterogeneities, as well as a framework for computing discrete and continuous adjoint wavefields which underpin the sensitivity kernels (Frechet derivatives) for full-waveform inversion with the same speedup [Szenzier et al., 2020]. While this capability to compute kernels is the most crucial computational aspect of full-waveform inversion, an actual iterative gradient scheme or other processing steps are not included here. 

<span style="color: red;"><b>*** TO DO: Missing Figure</b> </span>

![Wavefield snapshot for a simulation of the Virginia earthquake](https://www.earth.ox.ac.uk/wp-content/uploads/2020/07/scattering-768x538.jpg)

Future extensions (interested?)
Localised bases for azimuthal expansion; local time-stepping; GPU; improved load balancing 

> **Note:**
> AxiSEM3D is a community project. As such, we encourage contributions from the community to improve this code and tue manual over time.

 

## License

AxiSEM3D is published under the MIT license and open for any non-commercial usage. 

 ## Acknowledgements

The development of AxiSEM has been funded through a wide variety of grants to the authors. Initial development was supported under ...  
Contineud development

AxiSEM3D depends on the Salvus Mesher developed by Martin van Driel.


 ## Citing AxiSEM3D

In demonstrating continued relevance of this project to sponsors, we ask for you to cite the appropriate references if you publish results that were obtained to some part using AxiSEM3D. Acknwoledgement to the many who have contributed to the 
development of AxiSEM3D confers much desserved credit and helps the project receive funding.

In citing AxiSEM3D, please cite BOTH the code and the relevant published work.

CIte the code as:

    @article {

    }

Add the following to your list of References:

    Leng, Nissen-Meyer, van Driel, 2016, Efficient global wave propagation adapted to 3-D structural complexity: a pseudo- spectral/spectral-element approach, Geophysical Journal International, 207, 1700-1721. https://doi.org/10.1093/gji/ggw363
   
    Leng, Nissen-Meyer, van Driel, Hosseini, Al-Attar, 2019. AxiSEM3D: broad-band seismic wavefields in 3-D global earth models with undulating discontinuities, Geophysical Journal International, 217, 2125–2146, https://doi.org/10.1093/gji/ggz092

Bibtek:

    @article{Leng2016,
    title = {Efficient global wave propagation adapted to 3-D structural complexity: a pseudospectral/spectral-element approach},
    volume = {207},
    ISSN = {1365-246X},
    url = {http://dx.doi.org/10.1093/gji/ggw363},
    DOI = {10.1093/gji/ggw363},
    number = {3},
    journal = {Geophysical Journal International},
    publisher = {Oxford University Press (OUP)},
    author = {Leng,  Kuangdai and Nissen-Meyer,  Tarje and van Driel,  Martin},
    year = {2016},
    month = sep,
    pages = {1700–1721}
    }

    @article{Leng2019,
    title = {AxiSEM3D: broad-band seismic wavefields in 3-D global earth models with undulating discontinuities},
    volume = {217},
    ISSN = {1365-246X},
    url = {http://dx.doi.org/10.1093/gji/ggz092},
    DOI = {10.1093/gji/ggz092},
    number = {3},
    journal = {Geophysical Journal International},
    publisher = {Oxford University Press (OUP)},
    author = {Leng,  Kuangdai and Nissen-Meyer,  Tarje and van Driel,  Martin and Hosseini,  Kasra and Al-Attar,  David},
    year = {2019},
    month = feb,
    pages = {2125–2146}
    }


If relevant cite one or more of the following as well:

<span style="color: red;"><b>*** TO DO:ADD BIBTEKS</b></span>

 

	Leng, Korenaga, Nissen-Meyer, 2020. Three-dimensional scattering of elastic waves by small-scale heterogeneities in the Earth’s mantle, Geophysical Journal International, 223, 1, 502–525, https://doi.org/10.1093/gji/ggaa331

    Szenicer, Leng, Nissen-Meyer, 2020. A complexity-driven approach towards global waveform tomography, Geophysical Journal International, https://doi.org/10.1093/gji/ggaa349


    Haindl, Leng, Nissen-Meyer, 2021. A 3D Complexity-Adaptive Approach to Explore Sparsity in Visco-Elastic Wave Propagation, Geophysics, 86, 1, T331-T335, https://doi.org/10.1190/geo2020-0490.1

        @article{Haindl2021,
        title = {A 3D complexity-adaptive approach to explore sparsity in elastic wave propagation},
        volume = {86},
        ISSN = {1942-2156},
        url = {http://dx.doi.org/10.1190/geo2020-0490.1},
        DOI = {10.1190/geo2020-0490.1},
        number = {5},
        journal = {Geophysics},
        publisher = {Society of Exploration Geophysicists},
        author = {Haindl,  Claudia and Leng,  Kuangdai and Nissen-Meyer,  Tarje},
        year = {2021},
        month = may,
        pages = {T321–T335}
    }

    Tesoniero, Leng, Longs, Nissen-Meyer. 2020. Full wave sensitivity of SK(K)S phases to arbitrary anisotropy in the upper and lower mantle, Geophysical Journal International, 222, 412–435, https://doi.org/10.1093/gji/ggaa171

    Fernando, Leng, Nissen-Meyer, 2020. Oceanic high-frequency global seismic wave propagation with realistic bathymetry, Geophysical Journal International, 222,  1178–1194, https://doi.org/10.1093/gji/ggaa248 


### Data Availability
We strongly recommend making your data available for reproducibilty and replicability. Consider depositing your data (code, parameter files, data, log files, ...) in an approved repository, which will assign an identifier (e.g. a DOI) and enable
citation of your data. See [geodynamics.org software publishing guidance](https://geodynamics.org/software/software-bp/software-publishing).
Then add the following to your data availability statement:

    The code modifications, parameter, data, and log files used for the models in the study are available at DOI (Authors X, Y, Z) under CC BY-NC-SA 4.0.

    AxiSEM version X.X.X, (ADD CITATION) used in these computations is freely available under the MIT license through its software landing page https://geodynamics.org/resources/axisem3d and is being actively developed on GitHub and can be accessed via https://github.com/AxiSEMunity/AxiSEM3D.

### Acknowledgements
Please consider using the following text in your Acknowledgements section:

    We thank ....