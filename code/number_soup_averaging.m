% Number soup - A study of a simple ecosystem model
close all
clear
clc
% measuring running time
tic
rng(1677)
%rng(ceil(rand*10000)); %randomizing initial seed
numIter = 2500000; %number of iterations
a = 5; % constant affecting reproduction rate
dim = 6; %number of metabolites in system
u = 2; % influx metabolite
p = 0.001; % dictates the mutation probability
my = 500; %constant dictating the increase of metabolite
num_simulations = 1; % number of simulations, results are an average of these
dataFreq = 40; % storing frequency of data
rData = zeros(num_simulations,dim,ceil(numIter/dataFreq)); % metabolites in system
tData = zeros(num_simulations,1,ceil(numIter/dataFreq)); % time data

% species: vector storing number of individuals in each species
if(mod(dim,2) == 0)
    species = ones(1,(dim+1)*(dim/2));
else
    species = ones(1, dim + (dim)*((dim-1)/2));
end

speciesData = zeros(num_simulations,length(species),ceil(numIter/dataFreq)); % number of individuals in each species over time
numSpecies = length(species);
specSum = zeros(1,ceil(numIter/dataFreq)); % total species count data

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

for y = 1:num_simulations
    y
    % setting random seed
    %rng(y)
    % resetting resources, time,dt and #individuals
    r = zeros(1,dim); % metabolites in system
    t = 0;
    dt = 0;
    species(:) = 1;
    dataIndex = 1; % index for storing data
    
    % simulation
    for i=1:numIter
        %disp("timestep: "+i);
        % first  a random time interval dt is generated
        sumSpecies = sum(species);
        mu = 1/sumSpecies;
        dt = -mu .* log(rand);
        
        % metabolites are added
        r(u) = r(u) +my*dt;
        
        % a species is randomly selected to reproduce or die
        % the probability that species Sij is selected is Sij(t)/N(t)
        % where N(t) is the total population.
        randN = rand*sumSpecies;
        %disp("randN: " + randN);
        
        tempSum = 0;
        for j =1:numSpecies
            tempSum = tempSum + species(j);
            %disp("sum(species(1:j): " + sum(species(1:j)));
            if(tempSum >= randN)
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
        %disp(species)
        %disp("resources at i: ")
        %disp(r);
        %disp("added resources dt*my: " + dt*my);
        %disp(r);
        
        % The species selectedS reproduces with probability q(Rj(t),Ri(t))
        % or dies with probabiliy 1-q(Rj(t),Ri(t))
        if(si ~= sj)
            if(r(si) > 0 && r(sj) > 0)
                
                %disp("1 si != sj");
                qRjRi = r(si) * r(sj)/((a + r(si)) * (a + r(sj)));
                if(rand < qRjRi)
                    %disp("1.1: (si != sj) reproduction");
                    % the individual consumes metabolites
                    r(si) = r(si) -1;
                    r(sj) = r(sj) -1;
                    
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
                        % new individual mutates into another species.
                        speciesVec = 1:numSpecies;
                        speciesVec(selectedS) = [];
                        mutant = datasample(speciesVec,1);
                        species(mutant) = species(mutant) +1;
                        
                    end
                    
                else
                    %disp("1.2: (si != sj) the individual dies");
                    % since rand > qRjRi the individual dies
                    species(selectedS) = species(selectedS)-1;
                end
            else
                %disp("2: (si != sj) no metabolties, death");
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
                    speciesVec = 1:numSpecies;
                    speciesVec(selectedS) = [];
                    mutant = datasample(speciesVec,1);
                    species(mutant) = species(mutant) +1;
                    
                    % cruder but actually faster...
                    %while(true)
                    %   mutant = randi(numSpecies);
                    %    if(mutant ~= selectedS)
                    %        species(mutant) = species(mutant) +1;
                    %        %disp("mutation. new species:" + mutant);
                    %        break
                    %    end
                    %end
                    
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
        
        %disp("resources after timestep " + i +": ")
        %disp(r);
        %disp(" species after timestep " + i +": ")
        %disp(species);
        
        %disp("--------------------------")
        
        % time is updated and data is stored
        
        % time is updated and data is stored
        if(mod(i,dataFreq) == 0)
            tData(dataIndex) = t +dt;
            speciesData(y,:,dataIndex) = species;
            rData(y,:,dataIndex)=r;
            specSum(dataIndex) = sum(species);
            dataIndex = dataIndex +1;
        end
        t = t+dt;
        
        
        
        %w = input('press Enter');
    end
end

disp("lapsed runtime: " + toc);
%% plots

dilp1 = zeros(3,numIter+1);
dilp2 = zeros(3,numIter+1);
dilp3 = zeros(3,numIter+1);
dilp4 = zeros(3,numIter+1);
dilp1(:,:) = speciesData(1,:,:);
%dilp2(:,:) = speciesData(2,:,:);
%dilp3(:,:) = speciesData(3,:,:);
%dilp4(:,:) = speciesData(4,:,:);

%
figure(1)
%subplot(1,2,1)
area(tData(1,:),dilp1')
xlabel('time');
ylabel('nr individuals');
legend('species 11','species 21','species 22')

%subplot(1,2,2)
%area(tData(1,:),dilp2')
%xlabel('time');
%ylabel('nr individuals');
%legend('species 11','species 21','species 22')

% subplot(2,2,3)
% area(tData(1,:),dilp3')
% xlabel('time');
% ylabel('nr individuals');
% legend('species 11','species 21','species 22')

% subplot(2,2,4)
% area(tData(1,:),dilp4')
% xlabel('time');
% ylabel('nr individuals');
% legend('species 11','species 21','species 22')

%figure(2)
%plot(tData(1,:),rData(1,:,:))
%xlabel('time');
%ylabel('resources');
%legend('metabolite 1','mteabolite 2')