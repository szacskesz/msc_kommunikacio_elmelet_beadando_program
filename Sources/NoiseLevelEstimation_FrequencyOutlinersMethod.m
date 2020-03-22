function [nLevels] = NoiseLevelEstimation_FrequencyOutlinersMethod(img)  
    nLevels = zeros(1, size(img,3));
    for i = 1:size(img,3)
        channel = img(:,:,i);
        [c, level] = CustomARPENOS(channel);
        nLevels(i) = level;
    end
end
