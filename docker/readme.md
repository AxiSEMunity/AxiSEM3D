# Building local docker images

This directory contains Dockerfiles which can be used to build ubuntu containers containing AxiSEM3D.
These containers might be helpful if you are struggling to compile AxiSEM3D yourself.

The docker_copy directory contains a Dockerfile which copies the contains of your local Axisem3D directory
into the container, which might be useful if you are modifying or contributing new code to Axisem3D.

The docker_git directory contains a Dockerfile which clones the official Axisem3D repository,
which might be useful if you just want to use Axisem3D.

To build a docker image from the docker_git Dockerfile, run this command from the top level
of your Axisem3D directory:
docker build -f docker/docker_git/Dockerfile -t axisem3d --platform linux/amd64 .

To build a docker image from the docker_copy Dockerfile, run this command from the top level
of your Axisem3D directory:
docker build -f docker/docker_copy/Dockerfile -t axisem3d --platform linux/amd64 .

When building an image from the docker_copy Dockerfile it is especially important to run the build
command from the top level directory, because otherwise docker will not be able to copy all of your
Axisem3D files into the container.

# Using the Axisem3D image hosted on Dockerhub

The built docker image is also available on Dockerhub at geodynamics/axisem3d.
To use this prebuilt docker image, run these commands:

docker pull geodynamics/axisem3d;
docker run --rm -it geodynamics/axisem3d;
