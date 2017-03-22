# deep-photo-styletransfer
Code for paper "Deep Photo Style Transfer": [arXiv link coming soon]


## Setup
This code is based on torch. It has been tested on Ubuntu 14.04 LTS. 

Dependencies:
* [Torch](https://github.com/torch/torch7) (with [matio-ffi](https://github.com/soumith/matio-ffi.torch) and [loadcaffe](https://github.com/szagoruyko/loadcaffe))
* [Matlab](https://www.mathworks.com/)

CUDA backend:
* [CUDA](https://developer.nvidia.com/cuda-downloads)
* [cudnn](https://developer.nvidia.com/cudnn)

Download VGG-19:
```
sh models/download_models.sh
```

Compile ``cuda_utils.cu``
```
make clean && make
```

## Usage
To generate all results (in ``/examples`` folder) using the provided scripts, simply run 
```
run('gen_laplacian/gen_laplacian.m')
```
in Matlab and then
```
python gen_all.py
```
in Python. The final output will be in ``/examples/final_results`` folder.

### Basic usage
1. Given input and style images with semantic segmentation masks, put them in ``/examples`` folder respectively. They will have the filenames like: ``/examples/in/in<id>.png``, ``/examples/tar/tar<id>.png`` and ``/examples/segmentation/in<id>.png``, ``/examples/segmentation/tar<id>.png``;
2. 

## Examples












