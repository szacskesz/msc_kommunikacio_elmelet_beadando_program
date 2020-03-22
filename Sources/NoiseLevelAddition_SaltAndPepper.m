function [convertedImg] = NoiseLevelAddition_SaltAndPepper(img, value)
    convertedImg = img;
    convertedImg = imnoise(convertedImg, 'salt & pepper', 0.005 * value);
    
    convertedImg(convertedImg > 1) = 1;
    convertedImg(convertedImg < 0) = 0;
end
