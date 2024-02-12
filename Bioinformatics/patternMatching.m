pattern = 'ATAT';
genome = 'GATATATGCATATACTT';
position = [];
for i = 1 : (length(genome)-length(pattern)+1)
    if genome(i:i+length(pattern)) == pattern
        position.append(i)
    end 
end