clear all;
text = ['TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT'];
k = 3;
kmer = {};
freq = {};
n = length(text);
x = 0;
for i = 1 : (n-k+1)
    pattern = text(i:i+(k-1));
    if any(strcmp(kmer,pattern))
        %kmer{x} = pattern;
    else
        x = x + 1;
        kmer{x} = pattern;
    end
end 
for i = 1 : length(kmer)
    count = patternCount(text, cell2mat(kmer(i)));
    freq{i} = count; 
end
index = find([freq{:}] == max(cell2mat(freq)));
output = {};
for i = 1:length(index)
    output{i} = kmer{index(:,i)};
end
disp(output);