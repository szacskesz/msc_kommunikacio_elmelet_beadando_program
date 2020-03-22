function [convertedImg] = NoiseLevelRemoval_MedianMethod(img, patchSize)
    convertedImg = img;
    for i = 1:size(img,3)
        channel = img(:,:,i);
        c = medfilt2(channel, [patchSize, patchSize]);
        convertedImg(:,:,i) = c;
    end
end
