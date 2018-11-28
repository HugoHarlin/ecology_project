
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
figure(1)
histfit(death_clean);

figure(2)
data = hist(death_clean,30);
lambdahat = poissfit(data(2:end));
y = poisspdf(data(2:end),lambdahat);
plot(data(2:end),y,'+')


