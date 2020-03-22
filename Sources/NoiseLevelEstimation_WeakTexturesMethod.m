function [noise] = NoiseLevelEstimation_WeakTexturesMethod(img, patchSize)
    noise = NoiseLevel(double(img), patchSize);
end
