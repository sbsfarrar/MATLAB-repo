function [posStart,printPosStart] = approxPatternMatchingProblem(pattern,text,d)
%for the subsequent COUNTd(TEXT,PATTERN) problem:
        % simply find length(posStart)
posStart = [];
    for i = 1 : (length(text)-length(pattern)+1) %was +1, omitted for postion zero
        if hammingDistance(pattern,text(i:i+length(pattern)-1)) <= d
            posStart = [posStart, i-1]; %accomodate indexing starting at position zero
        end
    end
L = length(posStart);
printPosStart = "";
for i = 1 : L
        if i==1
            printPosStart = num2str(posStart(i));
        else
           printPosStart = printPosStart + " " + num2str(posStart(i)); 
        end
end