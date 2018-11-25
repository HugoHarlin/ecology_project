
len = length(death);
death_clean = [];
index = 1;
for i=1:len
    if(death(i) ~= -1)
        death(i)
        death_clean(index) = death(i);
        index = index +1;
    end
end
histfit(death_clean);