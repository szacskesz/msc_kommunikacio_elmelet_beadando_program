function [convertedImg] = NoiseLevelRemoval_FrequencyOutlinersMethod(img)
    convertedImg = img;
    for i = 1:size(img,3)
        channel = img(:,:,i);
        [c, level] = CustomARPENOS(channel);
        convertedImg(:,:,i) = c;
    end
end
