% Number soup - A study of a simple ecosystem model
%close all
clear
%clc
seed = ceil(rand*10000);
rng(seed);
%rng(1)
tic
numIter = 500000000; %number of iterations
a = 5; % constant affecting reproduction rate
dim = 2; %number of metabolites in system
u = 1; % influx metabolite
p = 0.001; % dictates the mutation probability
my = 500; %constant dictating the increase of metabolite
r = zeros(1,dim); % metabolites in system
rData = zeros(dim,numIter+1); % metabolites in system

tData = zeros(1,numIter+1);
dtData =zeros(1,numIter+1);
t = 0;
dt = 0;

if(mod(dim,2) == 0)
    species = ones(1,(dim+1)*(dim/2)); % number of species
else
    species = ones(1, dim + (dim)*((dim-1)/2));
end

numSpecies = length(species);
speciesData = zeros(numSpecies,numIter +1);


% the matrix sij stores which metabolites each species consumes
sij = zeros(2,numSpecies);
index = 1;
for i=1:dim
    for j=i:dim
        sij(1,index)=i;
        sij(2,index)=j;
        index = index +1;
    end
end

% simulation
for i=1:numIter
    %i
    %disp("timestep: "+i);
    % first  a random time interval dt is generated
    %dt = exprnd(1/sum(species));
    sumSpecies = sum(species);
    mu = 1/sumSpecies;
    dt = -mu .* log(rand);
    
    % saving dt data
    dtData(i) = dt;
    
    % metabolites are added
    r(u) = r(u) +my*dt;
    
    % a species is randomly selected to reproduce or die
    % the probability that species Sij is selected is Sij(t)/N(t)
    % where N(t) is the total population.
    randN = rand*sumSpecies;
    %disp("randN: " + randN);
    temp = 0;
    for j =1:numSpecies
        temp = temp + species(j);
        if(temp >= randN)
            selectedS = j;
            break
        end
    end
    
    si = sij(1,selectedS);
    sj = sij(2,selectedS);
    %disp("si = " + si);
    %disp("sj = " + sj);
    %disp("dt = " + dt);
    %disp("selectedS: " + selectedS)
    %disp(" species: ")
    % disp(species)
    % disp("resources at i: ")
    % disp(r);
    %disp("added resources dt*my: " + dt*my);
    
    % Species selectedS reproduces with probability q(Rj(t),Ri(t))
    % or dies with probabiliy 1-q(Rj(t),Ri(t))
    if(species(selectedS)> 0)
        if(si ~= sj)
            if(r(si) > 0 && r(sj) > 0)
                
                qRjRi = r(si) * r(sj)/((a + r(si)) * (a + r(sj)));
                if(rand < qRjRi)
                    %disp("1.1: (si != sj) reproduction");
                    % the individual consumes metabolites
                    r(si) = r(si)-1;
                    r(sj) = r(sj)-1;
                    
                    % the individual genereates new metabolites
                    if(si + sj < dim +1)
                        r(si+sj) = r(si+sj)+1;
                    elseif (si + sj == dim +1)
                        r(1) = r(1)+1;
                    else
                        r(1) = r(1)+1;
                        x= mod((si+sj),dim+1);
                        r(x) = r(x)+1;
                    end
                    % the individual reproduces, mutating to another species
                    % with probability p
                    if(rand > p)
                        species(selectedS) = species(selectedS)+1;
                    else
                        % the individual mutates
                        tempVec = 1:numSpecies;
                        tempVec(selectedS) = [];
                        mutant = datasample(tempVec,1);
                        species(mutant) = species(mutant) +1;
                    end
                    
                else
                    %disp("1.2: (si != sj) rand > qRjRi the individual dies");
                    % since rand > qRjRi the individual dies
                    species(selectedS) = species(selectedS)-1;
                end
            else
                %disp("2: (si != sj) not enough metabolties, death");
                % there aren't enough metabolties, the individual dies
                species(selectedS) = species(selectedS)-1;
            end
        elseif(r(si) > 1)
            % si == sj
            %disp("3: si = sj");
            qRjRi = r(si) * (r(sj)-1)/(  (a -1 + r(si) ) * (a + r(sj)));
            if(rand < qRjRi)
                %disp("3.1: (si = sj) reproduction");
                % the individual consumes metabolites
                r(si) = r(si)-2;
                
                % the individual genereates new metabolites
                if(si + sj < dim +1)
                    r(si+sj) = r(si+sj)+1;
                elseif (si + sj == dim +1)
                    r(1) = r(1)+1;
                else
                    r(1) = r(1)+1;
                    x= mod((si+sj),dim+1);
                    r(x) = r(x)+1;
                end
                % the individual reproduces, mutating to another species
                % with probability p
                if(rand > p)
                    species(selectedS) = species(selectedS)+1;
                else
                    % the individual mutates
                    tempVec = 1:numSpecies;
                    tempVec(selectedS) = [];
                    mutant = datasample(tempVec,1);
                    species(mutant) = species(mutant) +1;
                end
            else
                %disp("3.2 (si = sj) death");
                % since rand > qRjRi the individual dies
                species(selectedS) = species(selectedS)-1;
            end
        else
            %disp("4 (si = sj) death");
            % less than two metabolites, death
            species(selectedS) = species(selectedS)-1;
        end
    end
    %disp("resources after timestep " + i +": ")
    %disp(r);
    %disp(" species after timestep " + i +": ")
    %disp(species);
    
    %disp("--------------------------")
    
    % time is updated and data is stored
    tData(i+1) = t +dt;
    t = t+dt;
    speciesData(:,i) = species;
    rData(:,i)=r;
    
    %w = input('press Enter');
end
disp("lapsed runtime: " + toc);


%% plots
figure(2)
specTransp = speciesData';
area(tData(1:200:end),specTransp(1:200:end,:))
xlabel('time');
ylabel('nr individuals');
legend('species 11','species 21','species 22')

figure(3)
rTranspose = rData';
plot(tData(1:50000),rTranspose(1:50000,:))
xlabel('time');
ylabel('resources');
legend('metabolite 1','mteabolite 2')