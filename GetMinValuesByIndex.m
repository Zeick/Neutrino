function values = GetMinValuesByIndex(indexValues, y1, y2, y3, y4, y5)
    yall = [y1; y2; y3; y4; y5];
    values = zeros(1,length(indexValues));
    for j = 1:length(indexValues)
        values(j) = yall(indexValues(j), j);
    end
end