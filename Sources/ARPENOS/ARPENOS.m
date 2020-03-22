function im_noise_remov = ARPENOS( im , size_patches, reg)
% code for the paper 
% "Automated removal of quasiperiodic noise using
% frequency domain statistics"
% IS&T / SPIE Journal of Electronic Imaging
% vol. 24, no. 1, pages 013003/1-19
% Frederic Sur, Michel Grediac
% 2015
% usage:
% im: input image (double)
% size_patches: size of the patches (L in the paper)
% reg (optional argument): 0= no TV regularization, 1= constrained TV regularization

disp(' ')
disp('ARPENOS v. 1.2')
disp('Please cite:')
disp('Frederic Sur and Michel Grediac')
disp('Automated removal of quasiperiodic noise using frequency domain statistics')
disp('IS&T / SPIE Journal of Electronic Imaging')
disp('vol. 24, no. 1, pages 013003/1-19')
disp('2015.')
drawnow('update')

switch nargin
    case 2
        reg=0;
end;

% main constants (cf. paper)
step_patches=round(size_patches/8); % translation between patches (in pixels)
f2=8/size_patches; 
f0=f2/4;
f1=0.2;
upperlimit=3; % average spectrum coefficient considered as an outlier if distance to the regression line above upperlimit*sigma


% display power spectrum (log scale) and distribution
Fim=fft2(im);
Specim=abs(fftshift(Fim)).^2;
figure; imagesc(log(Specim)); colormap(gray), colorbar
title('power spectrum of the initial image (log scale)')
colorbarspectrum=[min(min(2*log(abs(fftshift(Fim))))), max(max(2*log(abs(fftshift(Fim)))))];
[rows, cols] = size(im);
[fx, fy] = meshgrid([-cols/2:cols/2-1]/cols, [-rows/2:rows/2-1]/rows);
Freq=sqrt(fx.^2+fy.^2); 
figure; loglog(Freq(:),Specim(:),'x','LineWidth',2)
xlabel('Frequency f (cycles per pixel)'), ylabel('|a_n(f)|^2')
grid on
title(num2str('Initial power spectrum distribution'));


% calculate average power spectrum
N=numel(1:step_patches:(rows-size_patches))*numel(1:step_patches:(cols-size_patches));
disp(' ')
disp(['Number of patches: ',num2str(N)])
H=hanning(size_patches)*hanning(size_patches)';
logspecmean=zeros(size_patches,size_patches);
for i=1:step_patches:(rows-size_patches),
  for j=1:step_patches:(cols-size_patches),
  imagette=H.*im(i:i+size_patches-1,j:j+size_patches-1);
  logspecmean=logspecmean+2/N*log(abs(fft2(imagette)));
  end;
end;
specmean=fftshift(exp(logspecmean));
logspecmean=fftshift(logspecmean);


% display average power spectrum 
figure; imagesc(logspecmean); colormap(gray), colorbar
title('average power spectrum on patches (log scale)')


% calculate power law parameters and display distribution
[fx, fy] = meshgrid([-size_patches/2:size_patches/2-1]/size_patches, [-size_patches/2:size_patches/2-1]/size_patches);
Freq=sqrt(fx.^2+fy.^2); 
J=find((Freq>f0)&(Freq<f1));
[coefb,stats]=robustfit(log(Freq(J)), log(specmean(J)));
std_b=stats.s;
Freqplot=[f0:0.01:f1];
Freqplotb=[f2:0.01:0.9];
figure, loglog(Freq(:),specmean(:),'x', Freqplot, exp(coefb(1))*Freqplot.^coefb(2),'r','LineWidth',2)
hold on 
loglog(Freqplotb,exp(coefb(1)+upperlimit*std_b)*Freqplotb.^coefb(2),'g','LineWidth',2)
xlabel('Frequency f (cycles per pixel)'), ylabel('|a_n(f)|^2')
grid on
title(num2str('Power spectrum distribution'));
disp(' ')
disp(['A = ',num2str(coefb(1))])
disp(['alpha = ',num2str(-coefb(2))])
disp(' ')


% Normalization of the average power spectrum
B1=coefb(1)*ones(size(Freq));
B2=coefb(2)*ones(size(Freq));
Model=B1+B2.*log(Freq);
figure; imagesc((logspecmean-Model)./std_b); colorbar, colormap(gray)
title('Normalized power spectrum')


% Outlier map 
Peaks=specmean.*(((logspecmean-Model)./std_b)>upperlimit);
Peaks(Freq<f2)=0;
figure; imagesc(Peaks>0); colormap(gray), colorbar;
title('Power spectrum outliers')


% Interpolation of the outlier map to the original image dimension
Peaks_b=zeros(size(Peaks)+2);
Peaks_b(2:end-1,2:end-1)=Peaks;
[fx, fy] = meshgrid([-size_patches/2-1:size_patches/2]/size_patches, [-size_patches/2-1:size_patches/2]/size_patches);
[x, y] = meshgrid([-cols/2:cols/2-1]/cols, [-rows/2:rows/2-1]/rows);
Out=interp2(fx,fy,Peaks_b,x,y);

% Notch filter desing
H = fspecial('gaussian',[12 12],2);
Notch_filter=sym_conv2(double(Out>0),H);
Notch_filter= fftshift(fft2(real(ifft2(ifftshift(Notch_filter))))); % (symmetrization to get rid of rounding errors)
M=max(Notch_filter(:));
if (M>0) Notch_filter=Notch_filter/M; end;


% Method noise estimation
Noise=real(ifft2(Fim.*ifftshift(Notch_filter)));


% Final display
figure; imagesc(2*log(abs(fftshift(Fim).*(1-Notch_filter))),colorbarspectrum); colormap(gray), colorbar;
title('Corrected power spectrum')
figure; imagesc(im,[0 255]); colormap(gray), colorbar
title('Original image')
im_noise_remov=im-Noise;
figure; imagesc(im_noise_remov,[0 255]); colormap(gray), colorbar;
title('Denoised image')
figure; imagesc(Noise), colormap(gray); colorbar
title('Method noise')

% Constrained TV regularization
if (reg)
    disp('Constrained TV regularization')
    disp(' ')
    niter=500; % max. iterations
    epsilon=0.1; % regularization parameter;
    Ima=im_noise_remov; % steepest descent initialization
    Gstep=1;
    Gx=diff(Ima([end 1:end],:),1,1); 
    Gy=diff(Ima(:,[end 1:end]),1,2); 
    regnorm=sqrt(epsilon^2+Gx.^2+Gy.^2);
    TV=mean(regnorm(:));
    disp(['Mean Total Variation before regularization: TV = ',num2str(TV)])
    for i=1:niter,
        tempx=Gx./regnorm; tempy=Gy./regnorm;
        Gxx=diff(tempx([1:end 1],:),1,1); 
        Gyy=diff(tempy(:,[1:end 1]),1,2); 
        divergence=Gxx+Gyy;
        GradProj=-ifft2(fft2(divergence).*ifftshift(Notch_filter),'symmetric'); % projected gradient
        Ima2=Ima-Gstep*GradProj; % constant step size
        Gx2=diff(Ima2([end 1:end],:),1,1); 
        Gy2=diff(Ima2(:,[end 1:end]),1,2); 
        regnorm2=sqrt(epsilon^2+Gx2.^2+Gy2.^2);
        TV2=mean(regnorm2(:)); % new value of the mean TV
        while (TV2>TV) 
            Gstep=Gstep/2; 
            disp([' Gstep = ',num2str(Gstep)]); 
            Ima2=Ima-Gstep*GradProj; 
            Gx2=diff(Ima2([end 1:end],:),1,1); 
            Gy2=diff(Ima2(:,[end 1:end]),1,2); 
            regnorm2=sqrt(epsilon^2+Gx2.^2+Gy2.^2);
            TV2=mean(regnorm2(:)); 
        end;
        if (TV-TV2<1e-4)
            break;
        else
            TV=TV2;
            regnorm=regnorm2;
            Ima=Ima2;
            Gx=Gx2;
            Gy=Gy2;
        end;
    end;
    disp(['Mean Total Variation after ',num2str(i),' iterations of gradient descent: TV = ',num2str(TV)])
    disp(' ')
    im_noise_remov_regul=Ima;
    figure; imagesc(2*log(abs(fftshift(fft2(Ima)))),colorbarspectrum); colormap(gray), colorbar;
    title('power spectrum of the regularized denoised image')
    figure; imagesc(im_noise_remov_regul,[0 255]); colormap(gray), colorbar;
    title('regularized denoised image')
    figure; imagesc(im-im_noise_remov_regul), colormap(gray); colorbar
    title('Difference between original image and regularized denoised image')
    figure; imagesc(im_noise_remov-im_noise_remov_regul), colormap(gray); colorbar
    title('Difference between denoised images before and after regularization')
    
    im_noise_remov=im_noise_remov_regul;
    
end

end

