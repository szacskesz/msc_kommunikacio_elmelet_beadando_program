function [convertedImg] = NoiseLevelAddition_Gauss(img, value)
    convertedImg = img;
    convertedImg = convertedImg + 0.003 * value * randn(size(convertedImg));
    
    convertedImg(convertedImg > 1) = 1;
    convertedImg(convertedImg < 0) = 0;
end
