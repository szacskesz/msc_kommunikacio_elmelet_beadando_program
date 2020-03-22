function [convertedImg, nLevel] = CustomARPENOS(im)

% main constants (cf. paper)
size_patches = 128;
step_patches=round(size_patches/8); % translation between patches (in pixels)
f2=8/size_patches; 
f0=f2/4;
f1=0.2;
upperlimit=3; % average spectrum coefficient considered as an outlier if distance to the regression line above upperlimit*sigma
[rows, cols] = size(im);


Fim=fft2(im);


% calculate average power spectrum
N=numel(1:step_patches:(rows-size_patches))*numel(1:step_patches:(cols-size_patches));
H=hanning(size_patches)*hanning(size_patches)';
logspecmean=zeros(size_patches,size_patches);
for i=1:step_patches:(rows-size_patches)
  for j=1:step_patches:(cols-size_patches)
  imagette=H.*im(i:i+size_patches-1,j:j+size_patches-1);
  logspecmean=logspecmean+2/N*log(abs(fft2(imagette)));
  end
end
specmean=fftshift(exp(logspecmean));
logspecmean=fftshift(logspecmean);


% calculate power law parameters and display distribution
[fx, fy] = meshgrid([-size_patches/2:size_patches/2-1]/size_patches, [-size_patches/2:size_patches/2-1]/size_patches);
Freq=sqrt(fx.^2+fy.^2); 
J=find((Freq>f0)&(Freq<f1));
[coefb,stats]=robustfit(log(Freq(J)), log(specmean(J)));
std_b=stats.s;


% Normalization of the average power spectrum
B1=coefb(1)*ones(size(Freq));
B2=coefb(2)*ones(size(Freq));
Model=B1+B2.*log(Freq);


% Outlier map 
Peaks=specmean.*(((logspecmean-Model)./std_b)>upperlimit);
Peaks(Freq<f2)=0;


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
if (M>0) 
    Notch_filter=Notch_filter/M; 
end


% Noise removal
Noise=real(ifft2(Fim.*ifftshift(Notch_filter)));
convertedImg = im-Noise;

% Noise level estimation
original = 2*log(abs(fftshift(Fim)));
original(original<0) = 0;
new = 2*log(abs(fftshift(Fim).*(1-Notch_filter)));
new(new<0) = 0;
nLevel = sum(abs(original - new), 'all');

end
