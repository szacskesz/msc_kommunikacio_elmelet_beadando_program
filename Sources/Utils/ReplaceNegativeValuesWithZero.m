function [arr] = ReplaceNegativeValuesWithZero(arr)
    arr(arr < 0) = 0;
end

