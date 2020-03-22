function [convertedImg] = NoiseLevelRemoval_WienerMethod(img, patchSize)
    convertedImg = img;
    for i = 1:size(img,3)
        channel = img(:,:,i);
        c = wiener2(channel, [patchSize, patchSize]);
        convertedImg(:,:,i) = c;
    end
end
