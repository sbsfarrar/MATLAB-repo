function [skewedGenome,minSkewed,locMinSkewed] = skewGenome(genome) %use '' for character vector
    %for text files, use 
    %file = fileread('/Users/Prog_ReS.S/Downloads/TEXTFILE.txt')%
    L = length(genome);
    skewed = (0);
    skewedGenome = "";
    for i = 1 : L
       if genome(i)=="A" || genome(i)=="T"
           skewed = [skewed, skewed(i)];
       else 
           if genome(i)=="G"
               skewed = [skewed, (skewed(i)+1)];
           else 
               if genome(i)=="C"
                   skewed = [skewed,(skewed(i)-1)];
               end
           end
       end
       disp((L-i) + " iterations remaining");
    end
    skewed(1) = [];
   L2 = length(skewed);
    disp("Generating Genome");
    for i = 1 : L2
        if i==1
            skewedGenome = num2str(skewed(i));
        else
           skewedGenome = skewedGenome + " " + num2str(skewed(i)); 
        end
        disp(i+" of "+L2);
    end
    %minSkewed%
    locMinSkewed = find(skewed(:)==min(skewed));
    %locMaxSkewed = find(skewed(:)==max(skewed));
    minSkewed = min(skewed);
    %maxSkewed = max(skewed);
end