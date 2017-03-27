# deep-photo-styletransfer
Code and data for paper "[Deep Photo Style Transfer](https://arxiv.org/abs/1703.07511)"


## Setup
This code is based on torch. It has been tested on Ubuntu 14.04 LTS.

Dependencies:
* [Torch](https://github.com/torch/torch7) (with [matio-ffi](https://github.com/soumith/matio-ffi.torch) and [loadcaffe](https://github.com/szagoruyko/loadcaffe))
* [Matlab](https://www.mathworks.com/) or [Octave](https://www.gnu.org/software/octave/)

CUDA backend:
* [CUDA](https://developer.nvidia.com/cuda-downloads)
* [cudnn](https://developer.nvidia.com/cudnn)

Download VGG-19:
```
sh models/download_models.sh
```

Compile ``cuda_utils.cu`` (Adjust ``PREFIX`` and ``NVCC_PREFIX`` in ``makefile`` for your machine):
```
make clean && make
```

## Usage
### Quick start
To generate all results (in ``examples/``) using the provided scripts, simply run
```
run('gen_laplacian/gen_laplacian.m')
```
in Matlab or Octave and then
```
python gen_all.py
```
in Python. The final output will be in ``examples/final_results/``.

### Basic usage
1. Given input and style images with semantic segmentation masks, put them in ``examples/`` respectively. They will have the following filename form: ``examples/input/in<id>.png``, ``examples/style/tar<id>.png`` and ``examples/segmentation/in<id>.png``, ``examples/segmentation/tar<id>.png``;
2. Compute the matting Laplacian matrix using ``gen_laplacian/gen_laplacian.m`` in Matlab. The output matrix will have the following filename form: ``gen_laplacian/Input_Laplacian_3x3_1e-7_CSR<id>.mat``;
3. Run the following script to generate segmented intermediate result:
```
th neuralstyle_seg.lua -content_image <input> -style_image <style> -content_seg <inputMask> -style_seg <styleMask> -index <id> -serial <intermediate_folder>
```
4. Run the following script to generate final result:
```
th deepmatting_seg.lua -content_image <input> -style_image <style> -content_seg <inputMask> -style_seg <styleMask> -index <id> -init_image <intermediate_folder/out<id>_t_1000.png> -serial <final_folder> -f_radius 15 -f_edge 0.01
```

### Image segmentation

Note: In the main paper we generate all comparison results using automatic scene segmenation algorithm modified from [DilatedNet](https://arxiv.org/abs/1606.00915). Manual segmentation enables more diverse tasks hence we provide the masks in ``examples/segmentation/``.

Here are some automatic or manual tools for creating a segmentation mask for a photo image:

#### Automatic:
MIT Scene Parsing: http://sceneparsing.csail.mit.edu/
SuperParsing: http://www.cs.unc.edu/~jtighe/Papers/ECCV10/
Nonparametric Scene Parsing: http://people.csail.mit.edu/celiu/LabelTransfer/
Berkeley: https://www2.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/resources.html
CRF-RNN for Semantic Image Segmentation: https://github.com/torrvision/crfasrnn
Selective Search: https://github.com/belltailjp/selective_search_py

#### Manual:
Photoshop Quick Selection Tool: https://helpx.adobe.com/photoshop/using/making-quick-selections.html

## Examples
Here are some results from our algorithm (from left to right are input, style and our output):
<p align='center'>
  <img src='examples/input/in3.png' height='194' width='290'/>
  <img src='examples/style/tar3.png' height='194' width='290'/>
  <img src='examples/final_results/best3_t_1000.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in4.png' height='194' width='290'/>
  <img src='examples/style/tar4.png' height='194' width='290'/>
  <img src='examples/final_results/best4_t_1000.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in13.png' height='194' width='290'/>
  <img src='examples/style/tar13.png' height='194' width='290'/>
  <img src='examples/final_results/best13_t_1000.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in9.png' height='194' width='290'/>
  <img src='examples/style/tar9.png' height='194' width='290'/>
  <img src='examples/final_results/best9_t_1000.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in20.png' height='194' width='290'/>
  <img src='examples/style/tar20.png' height='194' width='290'/>
  <img src='examples/final_results/best20_t_1000.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in1.png' height='194' width='290'/>
  <img src='examples/style/tar1.png' height='194' width='290'/>
  <img src='examples/final_results/best1_t_1000.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in39.png' height='194' width='290'/>
  <img src='examples/style/tar39.png' height='194' width='290'/>
  <img src='examples/final_results/best39_t_1000.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in57.png' height='194' width='290'/>
  <img src='examples/style/tar57.png' height='194' width='290'/>
  <img src='examples/final_results/best57_t_1000.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in47.png' height='194' width='290'/>
  <img src='examples/style/tar47.png' height='194' width='290'/>
  <img src='examples/final_results/best47_t_1000.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in58.png' height='194' width='290'/>
  <img src='examples/style/tar58.png' height='194' width='290'/>
  <img src='examples/final_results/best58_t_1000.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in51.png' height='194' width='290'/>
  <img src='examples/style/tar51.png' height='194' width='290'/>
  <img src='examples/final_results/best51_t_1000.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in7.png' height='194' width='290'/>
  <img src='examples/style/tar7.png' height='194' width='290'/>
  <img src='examples/final_results/best7_t_1000.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in23.png' width='290'/>
  <img src='examples/input/in23.png' width='290'/>
  <img src='examples/final_results/best23_t_1000.png' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in16.png' height='194' width='290'/>
  <img src='examples/style/tar16.png' height='194' width='290'/>
  <img src='examples/final_results/best16_t_1000.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in30.png' height='194' width='290'/>
  <img src='examples/style/tar30.png' height='194' width='290'/>
  <img src='examples/final_results/best30_t_1000.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in2.png' height='194' width='290'/>
  <img src='examples/style/tar2.png' height='194' width='290'/>
  <img src='examples/final_results/best2_t_1000.png' height='194' width='290'/>
</p>




## Acknowledgement
* Our torch implementation is based on Justin Johnson's [code](https://github.com/jcjohnson/neural-style);
* We use Anat Levin's Matlab [code](http://www.wisdom.weizmann.ac.il/~levina/matting.tar.gz) to compute the matting Laplacian matrix.

## Contact
Feel free to contact me if there is any question (Fujun Luan fl356@cornell.edu).

