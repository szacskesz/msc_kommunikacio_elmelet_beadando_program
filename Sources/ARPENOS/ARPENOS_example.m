% ARPENOS example script
% code for the paper 
% "Automated removal of quasiperiodic noise using
% frequency domain statistics"
% IS&T / SPIE Journal of Electronic Imaging
% vol. 24, no. 1, pages 013003/1-19
% Frederic Sur, Michel Grediac
% 2015


% to recreate Fig. 4 and 5 of the paper:
clear all, close all
im_gt=double(imread('mandril_gray.tif'));
[x,y]=meshgrid(1:size(im_gt,2),1:size(im_gt,1));
im=im_gt+(50*sin(2*pi*200/size(im_gt,1)*x).*sin(2*pi*2/size(im_gt,2)*y))+0*randn(size(im_gt));
im_noise_remov=ARPENOS(im,128);
diff=im_noise_remov(128:end-128,128:end-128)-im_gt(128:end-128,128:end-128);
disp(['RMSE: ',num2str(std(diff(:)))])


% to recreate Fig. 7 and 8 of the paper:
clear all, close all
im_gt=double(imread('boat.png'));
[x,y]=meshgrid(1:size(im_gt,2),1:size(im_gt,1));
im=im_gt+(50*sin(2*pi*100/size(im_gt,1)*x).*sin(2*pi*40/size(im_gt,2)*y))+0*randn(size(im_gt));
im_noise_remov=ARPENOS(im,128);
diff=im_noise_remov(128:end-128,128:end-128)-im_gt(128:end-128,128:end-128);
disp(['RMSE: ',num2str(std(diff(:)))]) 


% spectrum interpolation by TV minimization:
im_noise_remov=ARPENOS(im,128,1);
% note that ringing artifacts have disappeared
