# Requirements

To ensure a stable working environment, we recommend using [Docker](https://www.docker.com). To follow this recommendation, 
docker needs to be installed and running on your machine. We provide a Dockerfile based on CUDA 11.5 and Ubuntu 20.4.

If you want to use permGWAS2 without Docker, you need to install all packages mentioned in the 
[requirements file](../Docker/requirements.txt). 

# Installation Guide

1. Clone the repository into the directory where you want to set up the project

```shell
git clone https://github.com/grimmlab/permGWAS.git
```

2. To use permGWAS2 within a Docker environment, navigate to `Docker` and build a Docker image using the provided Dockerfile.

```shell
cd permGWAS/Docker
docker build -t IMAGENAME .
```

3. Run an interactive Docker container based on the created image.\
You have to mount the directory where the repository is located on your machine in the Docker container. 
If you want to work on GPU, specify the GPUs to mount.

```shell
docker run -it -v PATH_TO_REPO_FOLDER:/NAME_OF_DIRECTORY_IN_CONTAINER --gpus device=DEVICE_NUMBER --name CONTAINERNAME IMAGENAME
```

### Example

1. Assume our repository is located in a folder called `/myhome` and we want to name our image `permGWAS_image`

```shell
cd /myhome/permGWAS/Docker
docker build -t permGWAS_image .
```

2. Further, assume that we want to call our container `permGWAS_container`, our data is located in (subfolders of)
`/myhome` (i.e. we only need to mount one directory) and we want to use GPU 1. Then we have to run the following command:

```shell
docker run -it -v /myhome/:/myhome_in_container/ --gpus device=1 --name permGWAS_container permGWAS_image
```

3. If we need to mount a second directory (e.g. we want to save our results in a different folder called `/results`), 
we can run the following:

```shell
docker run -it -v /myhome/:/myhome_in_container/ -v /results/:/results/ --gpus device=1 --name permGWAS_container permGWAS_image
```

With this the setup is finished. For details on how to run permGWAS, see our [Quickstart Guide](./QUICKSTART.md).