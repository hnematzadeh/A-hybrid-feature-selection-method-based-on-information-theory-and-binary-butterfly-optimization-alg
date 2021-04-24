function [ ent ] = attributeEntropy( x, y )
% entropy of attribute x with classes in y
    ent = 0;
    
    if isa(x, 'cell')
        for v = unique(x)'
            idx = strcmp(x, v);
            prob = sum(idx) / length(x);
            if prob > 0
                ent = ent + prob * shannonEntropy(y(idx));
            end
        end
    elseif isa(x, 'numeric')
        [hc, edges, bin] = histcounts(x, 'BinMethod', 'sturges');
        for v = unique(bin)'
            idx = bin == v;
            prob = sum(idx) / length(x);
            if prob > 0
                ent = ent + prob * shannonEntropy(y(idx));
            end
        end
    else
        error(['cant handle attribute of class ', class(x)]);
    end
end

