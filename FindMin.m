% (C) Timo Karkkainen 2017
% Helper function.
% Compares two sets of values currentMin and newStuff, which must have
% same lengths. Returns a vector of same length where a smaller value is
% chosen.
function values = FindMin(currentMin, newStuff)

values = currentMin;
for j=1:length(values)
    if(values(j) > newStuff(j))
       values(j) = newStuff(j); 
    end
end