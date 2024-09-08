%
% Word_Find.m
%
% Finds scrambled English words of specified length in a set of letters.
%
% ==> Word_List.mat file must be loaded or be in the Current Folder. <==
%
% The list of 370k words is from: https://github.com/dwyl/english-words
%
% Kurt Oldenburg - 09/04/17 
%
%% Check for Word_List then input Letters and Word_Length
%
%  The inputs are set up to repeat the last value by just hitting Enter.

if ~exist('Word_List','var')
    load Word_List;
end
 
if ~exist('Letters','var')||isempty(Letters)
    Letters='abcdef';
end
Old_Letters=Letters;

if ~exist('Word_Length','var')||isempty(Word_Length)
   Word_Length=3;
end
Old_Length=Word_Length;

Letters=input(strcat('Letters (', Old_Letters(1,:), '):  '),'s');

 if isempty(Letters)
     Letters=Old_Letters;
 end 

if ~exist('Word_Length','var')||isempty(Word_Length)
    Word_Length=3;
end

Word_Length =input(strcat('Word Length (', num2str(Old_Length), '):  '));

if isempty(Word_Length)
   Word_Length=Old_Length;
end
Old_Length=Word_Length;

%% Scramble the Letters then find the words.
%
%  After creating all possible permutations of the letters, nchoosek() 
%  pulls out the all of the combos of the requested length. There will
%  be duplicates which are removed with unique()
% 
%  The intersect() function finds the combos that are in the list of
%  English words.

AllCombos=perms(Letters);
ComboPerms=[];

for I=1:size(AllCombos,1)
    ComboPerms=[ComboPerms;nchoosek(AllCombos(I,:),Word_Length)];
end

UniqueCombos=cellstr(unique(ComboPerms,'rows'));
Found_Words=intersect(UniqueCombos,Word_List);

clear AllCombos I ComboPerms UniqueCombos Old_Letters Old_Length

%%
%New List with seven character size exclusion
A = string(Word_List);
count = 0;
for i = 1: size(A,1)
    WordLength = strlength(A(i));
    if WordLength == 7
        count = count + 1;
        newList(count) = A(i);
    end
end

%%%---------------------------------------%%%
%SEARCH FOR FOUR-LETTER and THREE-LETTER WORDS FIRST
%Letter Iteration 1: R + U + D + E + 2 unknowns
clear newString_7
unknownLetters = ["a";"b";"c";"f";"g";"h";"i";"j";"k";"l";"m";"n";"o";"p";...
    "q";"s";"t";"v";"w";"x";"y";"z"];
patString = "rudee";
iter = 0;
for i = 1:size(unknownLetters,1)
    newString_6 = append(patString,unknownLetters(i));
    for t = 1:size(unknownLetters,1)
        if unknownLetters(t) == unknownLetters(i)
        else
            iter = iter + 1;
            newString_7(iter,1) = append(newString_6,unknownLetters(t));
        end
    end
end
%% BEER
unknownLetters = ["a";"c";"f";"g";"h";"i";"j";"k";"l";"m";"n";"o";"p";...
    "q";"s";"t";"v";"w";"x";"y";"z"]; %b omitted
patString = "ue";
iter = 0;
for i = 1:size(unknownLetters,1)
    newString_2 = append(patString,unknownLetters(i));
    for t = 1:size(unknownLetters,1)
        if unknownLetters(t) == unknownLetters(i)
        else
            iter = iter + 1;
            newString_3(iter,1) = append(newString_2,unknownLetters(t));
        end
    end
end
%%
%NOW FIND ALL THREE LETTER WORDS CREATED BY EACH 
Word_Length = 3;
Found_Words_Store = [];
for k = 1:size(newString_3,1)

    AllCombos=perms(newString_3{k});
    ComboPerms=[];

    for I=1:size(AllCombos,1)
        ComboPerms=[ComboPerms;nchoosek(AllCombos(I,:),Word_Length)];
    end

    UniqueCombos=cellstr(unique(ComboPerms,'rows'));
    Found_Words=intersect(UniqueCombos,Word_List);
    Found_Words_Store = [Found_Words_Store; Found_Words];
%     [s,~,t] = unique(Found_Words);
%     if isempty(t)
%     else
%         modeString(k,1) = string(s{mode(t)});
%     end
    clear AllCombos I ComboPerms UniqueCombos Old_Letters Old_Length
    fprintf('%d of %d \n',k,size(newString_3,1));
end
[s,~,t] = unique(Found_Words_Store);
tally = accumarray(t,1);
freqOccurrence = [cellstr(s), tally];
% pattern = "r";
% Index = contains(newList,pattern);
%% ORIGINAL
%NOW FIND ALL FOUR LETTER WORDS CREATED BY EACH 
Found_Words_Store = [];
Word_Length = 4;
for k = 1:size(newString_7,1)

    AllCombos=perms(newString_7{k});
    ComboPerms=[];

    for I=1:size(AllCombos,1)
        ComboPerms=[ComboPerms;nchoosek(AllCombos(I,:),Word_Length)];
    end

    UniqueCombos=cellstr(unique(ComboPerms,'rows'));
    Found_Words=intersect(UniqueCombos,Word_List);
    Found_Words_Store = [Found_Words_Store; Found_Words];
    clear AllCombos I ComboPerms UniqueCombos Old_Letters Old_Length
    fprintf('%d of %d \n',k,size(newString_7,1));
end
% pattern = "r";
% Index = contains(newList,pattern);
[s,~,t] = unique(Found_Words_Store);
%write for loop eliminating modes until top 20 words are selected