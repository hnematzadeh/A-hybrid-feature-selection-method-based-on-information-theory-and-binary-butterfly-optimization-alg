

function ig = computegains( Tr, yr )


    classEnt = shannonEntropy(yr);

    % compute entropy for attributes:
    attents = zeros(size(Tr,2), 1);
    gains = zeros(size(Tr,2), 1);

    for i = 1:size(Tr,2)
        attents(i) = attributeEntropy(Tr(:,i), yr);
        gains(i) = classEnt - attents(i);
    end
   
    ig = mean(gains);
end

