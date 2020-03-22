function [convertedImg] = NoiseLevelAddition_Speckle(img, value)
    convertedImg = img;
    convertedImg = imnoise(convertedImg, 'speckle', 0.005 * value);
    
    convertedImg(convertedImg > 1) = 1;
    convertedImg(convertedImg < 0) = 0;
end
