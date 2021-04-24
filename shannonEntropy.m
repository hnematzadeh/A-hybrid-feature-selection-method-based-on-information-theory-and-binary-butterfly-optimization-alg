% Author: Cássio M. M. Pereira <cassiomartini@gmail.com>
% 01/10/2015

function [ ent ] = shannonEntropy( x )

    if isa(x, 'numeric')
        hc = histcounts(x, 'BinMethod', 'sturges'); % estimate number of bins via sturges
    elseif isa(x, 'cell')
        T = tabulate(x);
        hc = cell2mat(T(:,2)); % get table count
    else
        error(['i dont know how to treat var of class ', class(x)]);
    end
    
    ps = hc ./ sum(hc); % probabilities
    ps = ps(ps > 0); % ignore bins with 0 prob, as 0*log(0)=0
    ent = -sum(ps .* log2(ps)); 
    
end
