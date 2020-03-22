function [convertedImg] = NoiseLevelRemoval_MeanMethod(img, patchSize)
    convertedImg = img;
    for i = 1:size(img,3)
        channel = img(:,:,i);
        c = filter2(fspecial('average', patchSize), channel);
        convertedImg(:,:,i) = c;
    end
end
