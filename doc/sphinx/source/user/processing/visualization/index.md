# Wavefield Visualization

Axisem3D can be configured to output seismic wavefields in two main ways: by stations or by elements. The stations options is perhaps the more traditional approach, in which a list of seismic station coordinates is provided, and various seismic wavefields are dumped at those stations. The elements option enables more general output geometries, such as dumping various seismic wavefields on slices through the planet. Outputting wavefields in this way allows for downstream applications such as wavefield animations or kernel computations.

In thissection we will go through both options in detail, and in particular explain the various parameters that affect the station and element outputs.

---

```{toctree}
---
maxdepth: 1
---
station.md
element.md
```
