# deep-photo-styletransfer
Code and data for paper "[Deep Photo Style Transfer](https://arxiv.org/abs/1703.07511)"

## Disclaimer 
**This software is published for academic and non-commercial use only.**

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

**Note: Please make sure that the content image resolution is consistent for Matting Laplacian computation in Matlab and style transfer in Torch, otherwise the result won't be correct.**

3. Run the following script to generate segmented intermediate result:
```
th neuralstyle_seg.lua -content_image <input> -style_image <style> -content_seg <inputMask> -style_seg <styleMask> -index <id> -serial <intermediate_folder>
```
4. Run the following script to generate final result:
```
th deepmatting_seg.lua -content_image <input> -style_image <style> -content_seg <inputMask> -style_seg <styleMask> -index <id> -init_image <intermediate_folder/out<id>_t_1000.png> -serial <final_folder> -f_radius 15 -f_edge 0.01
```

You can pass `-backend cudnn` and `-cudnn_autotune` to both Lua scripts (step 3.
and 4.) to potentially improve speed and memory usage. `libcudnn.so` must be in
your `LD_LIBRARY_PATH`. This requires [cudnn.torch](https://github.com/soumith/cudnn.torch).

### Image segmentation

Note: In the main paper we generate all comparison results using automatic scene segmentation algorithm modified from [DilatedNet](https://arxiv.org/abs/1606.00915). Manual segmentation enables more diverse tasks hence we provide the masks in ``examples/segmentation/``.

The mask colors we used (you could add more colors in `ExtractMask` function in two `*.lua` files):

| Color variable  | RGB Value | Hex Value | 
| ------------- | ------------- | ------------- |
| `blue`  | `0 0 255`  | `0000ff`  | 
| `green`  | `0 255 0`  | `00ff00`  |
| `black`  | `0 0 0`  | `000000`  |
| `white`  | `255 255 255`  | `ffffff`  |
| `red`  | `255 0 0`  | `ff0000`  |
| `yellow`  | `255 255 0`  | `ffff00`  |
| `grey`  | `128 128 128`  | `808080`  |
| `lightblue`  | `0 255 255`  | `00ffff`  |
| `purple`  | `255 0 255`  | `ff00ff ` |

Here are some automatic and manual tools for creating a segmentation mask for a photo image:

#### Automatic:
* [MIT Scene Parsing](http://sceneparsing.csail.mit.edu/)
* [SuperParsing](http://www.cs.unc.edu/~jtighe/Papers/ECCV10/)
* [Nonparametric Scene Parsing](http://people.csail.mit.edu/celiu/LabelTransfer/)
* [Berkeley Contour Detection and Image Segmentation Resources](https://www2.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/resources.html)
* [CRF-RNN for Semantic Image Segmentation](https://github.com/torrvision/crfasrnn)
* [Selective Search](https://github.com/belltailjp/selective_search_py)
* [DeepLab-TensorFlow](https://github.com/DrSleep/tensorflow-deeplab-lfov)

#### Manual:
* [Photoshop Quick Selection Tool](https://helpx.adobe.com/photoshop/using/making-quick-selections.html)
* [GIMP Selection Tool](https://docs.gimp.org/en/gimp-tools-selection.html)
* [GIMP G'MIC Interactive Foreground Extraction tool](http://gmic.eu/gimp.shtml)

## Examples
Here are some results from our algorithm (from left to right are input, style and our output):
<p align='center'>
  <img src='examples/input/in3.png' height='194' width='290'/>
  <img src='examples/style/tar3.png' height='194' width='290'/>
  <img src='examples/refine_posterization/refine_3.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in4.png' height='194' width='290'/>
  <img src='examples/style/tar4.png' height='194' width='290'/>
  <img src='examples/refine_posterization/refine_4.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in13.png' height='194' width='290'/>
  <img src='examples/style/tar13.png' height='194' width='290'/>
  <img src='examples/refine_posterization/refine_13.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in9.png' height='194' width='290'/>
  <img src='examples/style/tar9.png' height='194' width='290'/>
  <img src='examples/refine_posterization/refine_9.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in20.png' height='194' width='290'/>
  <img src='examples/style/tar20.png' height='194' width='290'/>
  <img src='examples/refine_posterization/refine_20.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in1.png' height='194' width='290'/>
  <img src='examples/style/tar1.png' height='194' width='290'/>
  <img src='examples/refine_posterization/refine_1.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in39.png' height='194' width='290'/>
  <img src='examples/style/tar39.png' height='194' width='290'/>
  <img src='examples/refine_posterization/refine_39.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in57.png' height='194' width='290'/>
  <img src='examples/style/tar57.png' height='194' width='290'/>
  <img src='examples/refine_posterization/refine_57.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in47.png' height='194' width='290'/>
  <img src='examples/style/tar47.png' height='194' width='290'/>
  <img src='examples/refine_posterization/refine_47.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in58.png' height='194' width='290'/>
  <img src='examples/style/tar58.png' height='194' width='290'/>
  <img src='examples/refine_posterization/refine_58.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in51.png' height='194' width='290'/>
  <img src='examples/style/tar51.png' height='194' width='290'/>
  <img src='examples/refine_posterization/refine_51.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in7.png' height='194' width='290'/>
  <img src='examples/style/tar7.png' height='194' width='290'/>
  <img src='examples/refine_posterization/refine_7.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in23.png' width='290'/>
  <img src='examples/input/in23.png' width='290'/>
  <img src='examples/final_results/best23_t_1000.png' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in16.png' height='194' width='290'/>
  <img src='examples/style/tar16.png' height='194' width='290'/>
  <img src='examples/refine_posterization/refine_16.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in30.png' height='194' width='290'/>
  <img src='examples/style/tar30.png' height='194' width='290'/>
  <img src='examples/refine_posterization/refine_30.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in2.png' height='194' width='290'/>
  <img src='examples/style/tar2.png' height='194' width='290'/>
  <img src='examples/final_results/best2_t_1000.png' height='194' width='290'/>
</p>

<p align='center'>
  <img src='examples/input/in11.png'  width='290'/>
  <img src='examples/style/tar11.png' width='290'/>
  <img src='examples/refine_posterization/refine_11.png'  width='290'/>
</p>


## Acknowledgement
* Our torch implementation is based on Justin Johnson's [code](https://github.com/jcjohnson/neural-style);
* We use Anat Levin's Matlab [code](http://www.wisdom.weizmann.ac.il/~levina/matting.tar.gz) to compute the matting Laplacian matrix.

## Citation
If you find this work useful for your research, please cite:
```
@article{luan2017deep,
  title={Deep Photo Style Transfer},
  author={Luan, Fujun and Paris, Sylvain and Shechtman, Eli and Bala, Kavita},
  journal={arXiv preprint arXiv:1703.07511},
  year={2017}
}
```

## Contact
Feel free to contact me if there is any question (Fujun Luan fl356@cornell.edu).

