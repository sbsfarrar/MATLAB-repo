%function [patternRC] = reverseComp(pattern)
pattern = 'GCTAGCT';
b = '';
    for i = 1 : length(pattern)
        if pattern(i) == 'A'
            a = 'T';
        else
            if pattern(i) == 'T'
                a = 'A';
            else
                if pattern(i) == 'C'
                a = 'G';
                else
                    if pattern(i) == 'G'
                    a = 'C';
                    end
                end 
            end
        end
        b = [a b]; 
        patternRC = b; 
    end
%   end
disp(patternRC)
 