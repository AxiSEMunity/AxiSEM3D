# Introduction

Seismic wave propagation stands at the centre of numerous endeavours in science and engineering, ranging from deep-earth tomography to seismic hazard assessment, nuclear monitoring, resource exploration, dynamic earthquake rupture and tsunami generation, volcano monitoring, earthquake early warning, non-destructive testing and extraterrestrial geophysics. Solving the underlying equations accurately in realistic 3D media at resolutions where observed data exists is a formidable challenge, ranking amongst the most demanding supercomputing problems.

AxiSEM3D is a versatile solver for wave propagation in elastic, viscoelastic, acoustic, fully anisotropic media in whole planets, spherical sections, or local-scale Cartesian domains. Its characteristic feature is a flexible and automated adaptivity to the complexity of underlying structures and wavefields, such that it can run up to 5 orders of magnitude faster in 1D or axisymmetric structures, and up to 2-3 orders of magnitude faster in smooth 3D structures compared to fully discretised 3D models with conventional methods (such as SEM). We call this Azimuthal Complexity Adaptation (ACA), such that AxiSEM3D is an ACA-SEM. The code simulates singular point sources (e.g. single forces, moment tensors, explosions) inside elastic and fluid media (including the atmosphere), but can also be used for simultaneous multiple sources. It boasts an extensive array of output options, including singular seismograms, surfaces, wavefields in 2D and 3D, each for a range of quantities such as displacement, acceleration, pressure, curl, stress tensor, deformation.

It has extra functionality for a wavefield injection approach 
{cite:t}`Leng2020`
to further boost computational efficiency for localised heterogeneities, as well as a framework for computing discrete and continuous adjoint wavefields which underpin the sensitivity kernels (Frechet derivatives) for full-waveform inversion with the same speedup
{cite:t}`Szenicer2020`. While this capability to compute kernels is the most crucial computational aspect of full-waveform inversion, an actual iterative gradient scheme or other processing steps are not included here.

Wavefield snapshot for a simulation of the Virginia earthquake (click to go to the YouTube video to watch the simulation):

[<img src="_static/wavefield_screenshot.png" width="500">](https://www.youtube.com/watch?v=v7_HqSzaBEg)

Interested in contributing? Here are some areas for future development—this is, of course, non-exhaustive: localized bases for azimuthal expansion; local time-stepping; GPU acceleration; and improved load balancing.

> **Note:**
> AxiSEM3D is a community project. As such, we encourage contributions from the community to improve this code and the manual over time.



## License

AxiSEM3D is published under the MIT license and open for any non-commercial usage.

 ## Acknowledgements

The development of AxiSEM has been funded through a wide variety of grants to the authors. Initial development was supported under ...
Continued development

AxiSEM3D depends on the Salvus Mesher developed by Martin van Driel.


 ## Citing AxiSEM3D

In demonstrating continued relevance of this project to sponsors, we ask for you to cite the appropriate references if you publish results that were obtained to some part using AxiSEM3D. Acknowledgment to the many who have contributed to the development of AxiSEM3D confers much deserved credit and helps the project receive funding.

In citing AxiSEM3D, please cite BOTH the code and the relevant published work.

Cite the code as:

    @software{leng_2026_20142089,
    author       = {Leng, Kuangdai and
                  Fernando, Benjamin and
                  Wolf, Jonathan and
                  Heister, Timo and
                  Hwang, Lorraine and
                  Nissen-Meyer, Tarje},
    title        = {AxiSEM3D v2.1.0},
    month        = may,
    year         = 2026,
    publisher    = {Zenodo},
    version      = {2.1.0},
    doi          = {10.5281/zenodo.20142089},
    url          = {https://doi.org/10.5281/zenodo.20142089},
    }

Add the following articles {cite:p}`Leng2016` and {cite:p}`Leng2019` to your list of References:

    Leng, Nissen-Meyer, van Driel, 2016, Efficient global wave propagation adapted to 3-D structural complexity: a pseudo- spectral/spectral-element approach, Geophysical Journal International, 207, 1700-1721. https://doi.org/10.1093/gji/ggw363

    Leng, Nissen-Meyer, van Driel, Hosseini, Al-Attar, 2019. AxiSEM3D: broad-band seismic wavefields in 3-D global earth models with undulating discontinuities, Geophysical Journal International, 217, 2125–2146, https://doi.org/10.1093/gji/ggz092

Bibtex:

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


If relevant, cite one or more of the following as well:

	@article{Leng2020,
    author = {Leng, K and Korenaga, J and Nissen-Meyer, T},
    title = {3-D scattering of elastic waves by small-scale heterogeneities in the Earth's mantle},
    journal = {Geophysical Journal International},
    volume = {223},
    number = {1},
    pages = {502--525},
    year = {2020},
    month = jul,
    issn = {0956-540X},
    doi = {10.1093/gji/ggaa331},
    url = {https://doi.org/10.1093/gji/ggaa331},
    eprint = {https://academic.oup.com/gji/article-pdf/223/1/502/33518142/ggaa331.pdf}
    }

    @article{Szenicer2020,
    author = {Szenicer, Alexandre and Leng, Kuangdai and Nissen-Meyer, Tarje},
    title = {A complexity-driven framework for waveform tomography with discrete adjoints},
    journal = {Geophysical Journal International},
    volume = {223},
    number = {2},
    pages = {1247--1264},
    year = {2020},
    month = jul,
    issn = {0956-540X},
    doi = {10.1093/gji/ggaa349},
    url = {https://doi.org/10.1093/gji/ggaa349},
    eprint = {https://academic.oup.com/gji/article-pdf/223/2/1247/33763352/ggaa349.pdf}
    }
          
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

    @article{Tesoniero2020,
    author = {Tesoniero, Andrea and Leng, Kuangdai and D. Long, Maureen and Nissen-Meyer, Tarje},
    title = {Full wave sensitivity of SK(K)S phases to arbitrary anisotropy in the upper and lower mantle},
    journal = {Geophysical Journal International},
    volume = {222},
    number = {1},
    pages = {412-435},
    year = {2020},
    month = {07},
    issn = {0956-540X},
    doi = {10.1093/gji/ggaa171},
    url = {https://doi.org/10.1093/gji/ggaa171},
    eprint = {https://academic.oup.com/gji/article-pdf/222/1/412/33181794/ggaa171.pdf},
    }

    @article{Fernando2020,
    author = {Fernando, Benjamin and Leng, Kuangdai and Nissen-Meyer, Tarje},
    title = {Oceanic high-frequency global seismic wave propagation with realistic bathymetry},
    journal = {Geophysical Journal International},
    volume = {222},
    number = {2},
    pages = {1178-1194},
    year = {2020},
    month = {08},
    issn = {0956-540X},
    doi = {10.1093/gji/ggaa248},
    url = {https://doi.org/10.1093/gji/ggaa248},
    eprint = {https://academic.oup.com/gji/article-pdf/222/2/1178/33341061/ggaa248.pdf},
    }   

### Data Availability
We strongly recommend making your data available for reproducibility and replicability. Consider depositing your data (code, parameter files, data, log files, ...) in an approved repository, which will assign an identifier (e.g., a DOI) and enable citation of your data. See [geodynamics.org software publishing guidance](https://geodynamics.org/software/software-bp/software-publishing).
Then add the following to your data availability statement:

    The code modifications, parameter, data, and log files used for the models in the study are available at DOI (Authors X, Y, Z) under CC BY-NC-SA 4.0.

    AxiSEM version X.X.X, (ADD CITATION) used in these computations is freely available under the MIT license through its software landing page https://geodynamics.org/resources/axisem3d and is being actively developed on GitHub and can be accessed via https://github.com/AxiSEMunity/AxiSEM3D.

### Acknowledgments
Please consider using the following text in your Acknowledgments section:

    We thank ....

### Getting Help

For questions about AxiSEM3D on all levels, please use the [AxiSEM3D forum](https://community.geodynamics.org/c/axisem).

See also [CONTRIBUTING.md](../../../CONTRIBUTING.md).

### Development
AxiSEM3D is being developed as a community code project. As such, we rely on the community to help maintain and contiribute new features to the code.

Development priorities arise from community consensus. The development priorites are maintained as separate [GitHub issues](https://github.com/AxiSEMunity/AxiSEM3D/issues?q=is%3Aissue%20state%3Aopen%20label%3ADEV%3Ahigh%2CDEV%3Alow%2CDEV%3Amedium) and maybe assigned high, medium, or low priority according to community ranking.
