#!/bin/bash
python -m salvus_mesh_lite.interface AxiSEMCartesian --basic.model "sfba_2018.bm" \
  --basic.period 2.0 --attenuation.frequencies 0.01 10.0 \
  --cartesian2Daxisem.x 50.0 --output_file AxiSEMCartesian_sfba_2018_2s.e
