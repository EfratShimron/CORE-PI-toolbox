# CORE-PI-toolbox
-----------------------------------------------------------------------------------
A Matlab package for CORE-PI - a parameter-free parallel MRI reconstruction method
----------------------------------------------------------------------------------

CORE-PI is a general reconstrution method, suitable for image reconstruction
from multi-coil (parallel imaging) acquisition of 2D Cartesian k-space
data. This method was published in:
     Shimron, Webb, Azhari, "CORE-PI: Non-iterative convolution-based 
     reconstruction for parallel MRI in the wavelet domain." 
     Medical Physics 46.1 (2019):199-214

CORE-PI is a parameter-free method, so you don't need to calibrate any params
It also enables flexible 1D undersampling of a 2D Cartesian k-space.
The toolbox includes demos with various undersampling schemes - periodic / 
varying-period / variable-density / random 

## Getting Started
Clone or download the CORE-PI code. 

### Prerequisites
To use the CORE-PI toolbox, you need a liscence for Matlab. The code was tested with Matlab2017R. 

## Running the examples
Open the "main.m" function in Matlab, choose one example from the list, and run the code.

There are 9 reconstruction examples, divided to 3 groups:
I. Recosntructions of an analyltical brain phantom using different subsampling schemes - periodic / 
   varying-period / variable-density / random subsampling - all with a reduction factor of R=6, and 
   with the default wavelet type 'db2'. 
II. Reconstructions of the same brain phantom using different wavelet types - haar / coif1 / sym4.
III. Reconstruction examples with two in-vivo 7t brain scans data, subsampled periodically with R=4. 

Results for examples in group I:

![examples with different subsampling schemes](https://github.com/EfratShimron/CORE-PI-toolbox/blob/master/README_figures/phantom_examples.png)


## Acknowledgments
The in-vivo data is courtesy of Prof. Andrew G. Webb from Leiden University Medical Center (LUMC). 

The Realistic Analytical Brain Phantom (resterized) data was utilized here with permission from
the authors of:
    Guerquin-Kern, Matthieu, et al. "Realistic analytical phantoms for parallel 
    magnetic resonance imaging." IEEE Transactions on Medical Imaging 31.3
    (2011): 626-636.
If you use that data in your publications, please cite this paper. 
