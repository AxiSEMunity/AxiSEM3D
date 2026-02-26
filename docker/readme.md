# Building local docker images

This directory contains Dockerfiles which can be used to build ubuntu containers containing AxiSEM3D.


To build a docker image of AxiSEM3D, run this command from the top level
of your Axisem3D directory:
docker build -f docker/Dockerfile -t axisem3d --platform linux/amd64 .

# Using the Axisem3D image hosted on Dockerhub

The built docker image is also available on Dockerhub at geodynamics/axisem3d.
To use this prebuilt docker image, run these commands:

```
  docker pull geodynamics/axisem3d
  docker run --rm -it geodynamics/axisem3d
```
