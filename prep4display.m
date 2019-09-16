function im_final = prep4display(im,extra_fftshift_flag,rotangle,mask4display)
% This function prepares a reconstructed image for display.
% It performs shifting, rotation, masking and taking absolute values.

% Compute image magnitude
im_final = abs(im);

% fftshift (such a 1D fftshift is needed for our data, obtained from the Phillips scanner at LUMC, Leiden University)
if extra_fftshift_flag==1
    im_final = fftshift(im_final,2);
end

% rotation
if rotangle~=0
    im_final = imrotate(im_final,rotangle);
end

% masking
im_final = im_final.*mask4display;

