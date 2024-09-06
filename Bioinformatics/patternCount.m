function [count] = patternCount(Text, Pattern)
  count = 0;
  for i = 0 : (length(Text)-length(Pattern))
    if Text(i+1:(length(Pattern)+i)) == Pattern
      count = count + 1;
    end 
  end 