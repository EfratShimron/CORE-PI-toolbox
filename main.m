clear all;
close all;

% This Matlab code is an implementation of the CORE-PI method published in:
% Shimron, Efrat, Andrew G. Webb, and Haim Azhari. "CORE?PI: Non?iterative 
% convolution?based reconstruction for parallel MRI in the wavelet domain." 
% Medical physics 46.1 (2019): 199-214.?

% Also, this toolbox contains data that was produced using the Realistic
% Analytical Brain Phantom, which was published in:
% Guerquin-Kern, Matthieu, et al. "Realistic analytical phantoms for parallel 
% magnetic resonance imaging." IEEE Transactions on Medical Imaging 31.3 (2011): 626-636.?

% Please kindly include the aforementioned citations whenever you present or 
% publish results that are based on this toolbox.

% (c) E. Shimron, H. Azhari, 2019

% ====================================================
%      CHOOSE ONE EXAMPLE FROM THE FOLLOWING LISTS 
% ====================================================
% ---- examples with different subsampling schemes ----
demo = 'brain_phantom_example';  sampling_scheme='periodic';          wavelet_type = 'db2';
%demo = 'brain_phantom_example';  sampling_scheme='variying-period';   wavelet_type = 'db2';
%demo = 'brain_phantom_example';  sampling_scheme='variable-density';  wavelet_type = 'db2';
%demo = 'brain_phantom_example';  sampling_scheme='random';            wavelet_type = 'db2';

% --- examples with different wavelet types ---
 demo = 'brain_phantom_example';  sampling_scheme='periodic';   wavelet_type = 'haar';  %  Try different wavelet types: 'haar' / 'db5' / 'sym4' / 'coif1' (see fig. 5)
%demo = 'brain_phantom_example';  sampling_scheme='periodic';   wavelet_type = 'coif1';  %  Try different wavelet types: 'haar' / 'db5' / 'sym4' / 'coif1' (see fig. 5)
%demo = 'brain_phantom_example';  sampling_scheme='periodic';   wavelet_type = 'sym4';  %  Try different wavelet types: 'haar' / 'db5' / 'sym4' / 'coif1' (see fig. 5)

% ---- examples with in-vivo data --------
% demo = 'In_vivo_example_1';      sampling_scheme='periodic';          wavelet_type = 'db2';
% demo = 'In_vivo_example_2';      sampling_scheme='periodic';          wavelet_type = 'db2';


% NOTE: this toolbox currently supports various types of under-sampling
% for the brain phantom data, and only periodic under-sampling for
% the in-vivo data. 

% ================ preparations load k-space data & sensitivity maps  ================
D = DataProcess(demo,sampling_scheme,wavelet_type);

% ================ display sampling mask  ================
figure;
imshow(D.KspaceSampPattern_DC_in_center); axis equal; axis tight; axis off;
title_str = [sampling_scheme,' Sampling, R=',num2str(D.R)];
title(title_str,'FontSize',12);  colormap (gray);

% ================ display gold standard image ================
figure; imagesc(D.GoldStandard4display); title(['Gold Standard']); caxis([D.cmin D.cmax]); axis off; colormap (gray); axis image;

% ============================================
%                     CORE-PI               
% ============================================
% compute the CORE-PI reconstruction
D = CORE_PI(D);

% ======== display wavelet-domain coeffs ======

% concatenate matrices for visualizing
SWT_Rec_MAT_CORE_PI = [D.conv_image_LP_channel_4display  ones(D.N,5) D.conv_image_HP_channel_4display  ];

figure; imagesc(abs([SWT_Rec_MAT_CORE_PI])); axis off; axis image; colormap gray; caxis([0 D.cmax]);
title(['Low-Pass (approximation)      High-Pass (details)      '])
suptitle('Reconstructed SWT decomposition')

% ========= Calc error image & NRMSE ========
err_mat = abs(abs(D.GoldStandard4display)- abs(D.CORE_PI_Rec4display));
NRMSE = calc_NRMSE(D.GoldStandard4display,D.CORE_PI_Rec4display);

% ======== display Gold Standard + Rec + Error ======
MAT = [D.GoldStandard4display   ones(D.N,5) D.CORE_PI_Rec4display ; ones(2,5+2*D.N); ones(D.N,D.N)  ones(D.N,5) err_mat*4];

figure; imagesc(abs(MAT)); axis off; axis image; colormap gray; caxis([0 D.cmax]);
text(10,10,'Gold Standard','Color','w')
text(10+D.N,10,'CORE-PI','Color','w')
text(10+D.N,D.N+2+10,'Error magnified x4','Color','w');
text(10+D.N,2*D.N-10,sprintf('NRMSE=%.5f',NRMSE),'Color','w');

