function [count] = approxPatternCount(text, pattern, d)
count = 0;
for i = 1 : (length(text)-length(pattern)+1)
    kmer = text(i:i+length(pattern)-1);
    if hammingDistance(pattern,kmer) <= d
        count = count + 1;
    end
end
end