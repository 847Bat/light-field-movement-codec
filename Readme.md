# Light field translation codec

Table of content:
- Description & Terms of Use
- Installation
- Usage
- Remarks

## Description & Terms of Use
Light field technology has recently been gaining traction in the research community. Several acquisition technologies have been demonstrated to properly capture light field information, and portable devices have been commercialized to the general public. However, new and efficient compression algorithms are needed to sensibly reduce the amount of data that needs to be stored and transmitted, while maintaining an adequate level of perceptual quality. 
In this paper, we propose a novel light field compression scheme that uses view estimation to recover the entire light field from a small subset of encoded views.
Experimental results on a widely used light field dataset show that our method achieves good coding efficiency with average rate savings of 54.83% with respect to HEVC.

The codec is implemented in MATLAB. If you use this code in your research, please cite our paper:

B. Hériard-Dubreuil, I. Viola and T. Ebrahimi. "Light field compression using translation-assisted view estimation". 2019. 34th Picture Coding Symposium (PCS), Ningbo, China, November 12-15, 2019. 

## Installation

To begin with, clone this repository in a working folder.
```git clone insert_path_here```

You will need a complementary library: the [Light Field Toolbox v0.4](https://ch.mathworks.com/matlabcentral/fileexchange/49683-light-field-toolbox-v0-4). Please download it and move the folder **light-field-toolbox** to **./modules**.

You will also need a HEVC encoder. Please download [HEVC 15.0 + Rext 8.1](https://hevc.hhi.fraunhofer.de/svn/svn_HEVCSoftware/tags/HM-15.0+RExt-8.1/) and move the folder **HM-15.0+RExt-8.1** to the location **../**

Finally, you need to add some light field images. To do so, please add such a picture in the folder **./refs/image_name/**, view by view, each view being a 10 bits RGB PPM picture of the name **%03d_%03d.ppm** (example: 001_004.ppm for the view (1,4)).

## Usage

You will find a file **test_global.m**, divided in cells to load the data, select the parameters, encode, visualize the bitrate, decode and compute the quality metrics. Change the img name and the parameters to your will. You will find a description of the parameters usage in the header of **encoder.m**.

The code is supposed to check in your repository if you already compressed this set of references with the same parameters, and if not, compress them with HEVC.

## Remarks

The code has not been prepared for a different image size (434x626). If you want to change that, you will need to modify **refs_compression.m**.
