function [indexValues, values] = FindMinIndex(y1,y2,y3,y4,y5)
    indexValues = zeros(1,length(y1));
    values = zeros(1,length(y1));
    for j = 1:length(y1)
        [minValue, minIndex] = min([y1(j) y2(j) y3(j) y4(j) y5(j)]);
        indexValues(j) = minIndex;
        values(j) = minValue;
    end
end