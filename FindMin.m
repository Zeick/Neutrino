% (C) Timo Karkkainen 2017
% Helper function.

function values = FindMin(currentMin, newStuff)

values = currentMin;
for j=1:length(values)
    if(values(j) > newStuff(j))
       values(j) = newStuff(j); 
    end
end