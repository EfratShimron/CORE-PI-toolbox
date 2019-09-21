# CORE-PI-toolbox
-----------------------------------------------------------------------------------
A Matlab package for CORE-PI - a parameter-free parallel MRI reconstruction method
----------------------------------------------------------------------------------

CORE-PI is a general reconstrution method, suitable for image reconstruction
from multi-coil (parallel imaging) acquisition of 2D Cartesian k-space
data. 

CORE-PI was published in:
     Shimron, Webb, Azhari, "CORE-PI: Non-iterative convolution-based 
     reconstruction for parallel MRI in the wavelet domain." 
     Medical Physics 46.1 (2019):199-214

CORE-PI is user-friendly: it is a parameter-free linear (non-iterative) method. 
Thus, users do not need to calibrate any params!

CORE-PI enables flexible undersampling of a 2D Cartesian k-space with 1D undersampling 
schemes.

Exmaples include: 
(1) Recosntructions of a realistic analyltical brain phantom using different subsampling schemes.
(2) Examples of in-vivo 7t brain images reconstruction



![brain phantom - four examples with different subsampling schemes]./images/CORE-PI phantom examples.jpg



Note - the Realistic Analytical Brain Phantom data was reproduced with
    permission from the authors of this paper:
    Guerquin-Kern, Matthieu, et al. "Realistic analytical phantoms for parallel 
    magnetic resonance imaging." IEEE Transactions on Medical Imaging 31.3
    (2011): 626-636.
