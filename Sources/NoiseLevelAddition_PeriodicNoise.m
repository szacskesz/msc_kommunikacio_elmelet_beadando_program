function [convertedImg] = NoiseLevelAddition_PeriodicNoise(img, value)
    convertedImg = img;
    alpha = 40; % 100
    beta = 40;   % 40
    
    value = value / 100 * 0.25;

    for i = 1:size(img,3)
        [x,y] = meshgrid(1:size(convertedImg, 2), 1:size(convertedImg,1));
        convertedImg(:,:,i) = convertedImg(:,:,i) + value * sin(2*pi*alpha/size(convertedImg,1)*x);
    end
    
    convertedImg(convertedImg > 1) = 1;
    convertedImg(convertedImg < 0) = 0;
end
