function [convertedImg] = ConvertImageToGrayscale(img)
    convertedImg = img;
    convertedImg = rgb2gray(convertedImg);
end
