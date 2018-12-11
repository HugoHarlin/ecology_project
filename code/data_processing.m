close all;


len = length(death);
death_clean = [];
index = 1;
for i=1:len
    if(death(i) ~= -1)
        death_clean(index) = death(i);
        index = index +1;
    end
end

%% --------------- removing early collapses ---------------------
len = length(death_clean);
death_cleanClean = [];
index = 1;
for i=1:len
    if(death_clean(i) > 100)
        death_cleanClean(index) = death_clean(i);
        index = index +1;
    end
end

%% ------------ FIGURES --------------

figure(1)
h = histfit(death_cleanClean, 30,'gamma');
h(1).FaceColor = [0 100/256 208/256];
xlabel('Generations, t','FontWeight','bold')
ylabel('Collapsed Ecosystems','FontWeight','bold');
set(gca,'fontsize', 20)
phat = gamfit(death_cleanClean);
hold on
legend("data","gamma distribution");
%
%hist(death_cleanClean,30)
%h = findobj(gca,'Type','patch');
%h.FaceColor = [0 100/256 208/256];

xlabel('Generations, t','FontWeight','bold')
ylabel('Collapsed Ecosystems','FontWeight','bold');
%set(gca,'fontsize', 20)

%%
figure(2)
%edges = [0:2:100,200:80000/30:80000];
edges = [0:(150000/30):150000];
histogram(death_cleanClean,edges);
xlabel('Generations, t','FontWeight','bold')
ylabel('Collapsed Ecosystems','FontWeight','bold');
set(gca,'fontsize', 15)

nbins = 30;
[N,edges] = histcounts(death_cleanClean,nbins);



%% code for hist fit plot
xlabel('Generations, t','FontWeight','bold')
ylabel('Collapse Density','FontWeight','bold');
set(gca,'fontsize', 15)

