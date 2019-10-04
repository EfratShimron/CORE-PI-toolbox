# ![alt text](https://github.com/EfratShimron/CORE-PI-toolbox/blob/master/README_figures/CORE_PI%20logo%20v2.png)


<img src="https://github.com/EfratShimron/CORE-PI-toolbox/blob/master/README_figures/CORE_PI%20logo%20v2.png" width=200 align=left>


-----------------------------------------------------------------------------------
A Matlab package for CORE-PI - a parameter-free parallel MRI reconstruction method
----------------------------------------------------------------------------------

CORE-PI is a general reconstrution method, suitable for image reconstruction
from multi-coil (parallel imaging) acquisition of 2D Cartesian k-space
data. This method was published in:
     [Shimron, Webb, Azhari, "CORE-PI: Non-iterative convolution-based 
     reconstruction for parallel MRI in the wavelet domain." 
     Medical Physics 46.1 (2019):199-214](https://aapm.onlinelibrary.wiley.com/doi/full/10.1002/mp.13260)


CORE-PI is a **parameter-free method**, so users do not need to calibrate any params!


CORE-PI is enables flexible 1D undersampling of a 2D Cartesian k-space.
The toolbox includes demos with various undersampling schemes - periodic, varying-period, variable-density and random schemes.

## Getting Started
Clone or download the CORE-PI code. 

## Prerequisites
A liscence for Matlab is required. The code was tested with Matlab2017R. 

## Running the examples
Open the "main.m" function in Matlab, choose one example from the list, and run the code.

There are 9 reconstruction examples, divided to 3 groups:


1. **Analyltical brain phantom demos** - reconstructions were performed for **different subsampling schemes**, 
   all with a reduction factor (sub-sampling rate) of R=6:
   - Periodic 
   - Varying-period
   - Variable-density 
   - Random subsampling 
   In all these demos CORE-PI was impelmented with wavelet 'db2'.
   

2. **Analyltical brain phantom demos** in which CORE-PI was implemented using **different wavelet types**:
    - haar 
    - coif1 
    - sym4.


3. **In-vivo 7t brain scans demos** - data was retrospectivly subsampled with R=4 (using periodic subsampling),
     and CORE-PI was impelmented with wavelet 'db2'.


## Analaytical Brain Phantom demos - different subsampling schemes

Results for examples in group 1:

![examples with different subsampling schemes](https://github.com/EfratShimron/CORE-PI-toolbox/blob/master/README_figures/phantom_examples.png)


## In-vivo 7T Brain Scans Demos 

Example from group 3:

![examples with different subsampling schemes](https://github.com/EfratShimron/CORE-PI-toolbox/blob/master/README_figures/CORE-PI_in_vivo_fig.png)


## Acknowledgments
The in-vivo data is courtesy of Prof. Andrew G. Webb from Leiden University Medical Center (LUMC). 

The Realistic Analytical Brain Phantom (resterized) data was utilized here with permission from
the authors of:
    Guerquin-Kern, Matthieu, et al. "Realistic analytical phantoms for parallel 
    magnetic resonance imaging." IEEE Transactions on Medical Imaging 31.3
    (2011): 626-636.
If you use that data in your publications, please cite this paper. 
