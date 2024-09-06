function [dist] = hammingDistance(p,q)
    dist = 0;
    L = length(p);
    for i = 1 : L
        if p(i) == q(i)
            continue 
        else 
            dist = dist + 1;
        end
    end
end