% Number soup - A study of a simple ecosystem model
close all
clear
clc

numGen = 100000; %number of (maximal) generations, can be lower if ecosystem dies out
a = 5; % constant affecting reproduction rate
dim =4; %number of metabolites in system
u = 4; % influx metabolite
p = 0.001; % dictates the mutation probability
my = 500; %constant dictating the increase of metabolite
numRuns = 1000; % number of runs
dataFreq = 1; % storing frequency of data
death = ones(1,numRuns); % list for saving survival time length, 0 if no extinction
death = -1.*death; % all simulations that dont die out in numGen generations are -1 by default

tic
parfor run = 1:numRuns
%for run = 1:numRuns
    selectedS = 0;
    rng(run,'twister');
    
    r = zeros(1,dim); % metabolites in system
    t = 0;
    dt = 0;
    
    if(mod(dim,2) == 0)
        species = ones(1,(dim+1)*(dim/2)); % number of species
    else
        species = ones(1, dim + (dim)*((dim-1)/2));
    end
    
    numSpecies = length(species); % number of species
%     speciesData = zeros(numSpecies,numGen); %species data
%     rData = zeros(dim,numGen); % metabolites in system
%     tData = zeros(1,numGen); % time data
%     specSum = zeros(1,numGen); % total species count data
%     dataIndex = 2; % index for storing data, index 1 is the initial conditions
    
    
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
    i = 1;
    while (t < numGen)
        %i
        %disp("timestep: "+i);
        % first  a random time interval dt is generated
        %dt = exprnd(1/sum(species));
        %if(t> (numGen-0.001))
        %    death(run) = -1*t;
        %end
        
        sumSpecies = sum(species);
        if(sumSpecies == 0)
            death(run) = t;
            break
        end
        
        mu = 1/sumSpecies;
        dt = -mu .* log(rand);
        
        if(sumSpecies ~= 0)
            mu = 1/sumSpecies;
            dt = -mu .* log(rand);
        else
            break
        end
        
        % metabolites are added
        r(u) = r(u) +my*dt;
        
        % a species is randomly selected to reproduce or die.
        % The probability that species Sij is selected is Sij(t)/N(t)
        % where N(t) is the total population.
        randN = rand*sumSpecies;
        %disp("randN for species selection: " + randN);
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
        %disp(species)
        %disp("resources at i: ")
        %disp(r);
        %disp("added resources dt*my: " + dt*my);
        
        % Species selectedS reproduces with probability q(Rj(t),Ri(t))
        % or dies with probabiliy 1-q(Rj(t),Ri(t))
        if(species(selectedS)> 0)
            if(si ~= sj)
                if(r(si) > 0 && r(sj) > 0)
                    
                    qRjRi = r(si) * r(sj)/((a + r(si)) * (a + r(sj)));
                    randn = rand;
                    %disp("qRjRi: " + qRjRi);
                    %disp("randn: " + randn);
                    if(randn < qRjRi)
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
                        randn = rand;
                        %disp("mutation p: " + p);
                        %disp("randn: " + randn);
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
                
                randn = rand;
                %disp("qRjRi: " + qRjRi);
                %disp("randn: " + randn);
                
                if(randn < qRjRi)
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
                    randn = rand;
                    %disp("mutation p: " + p);
                    %disp("randn: " + randn);
                    if(randn > p)
                        species(selectedS) = species(selectedS)+1;
                    else
                        % the individual mutates
                        tempVec = 1:numSpecies;
                        tempVec(selectedS) = [];
                        mutant = datasample(tempVec,1);
                        species(mutant) = species(mutant) +1;
                    end
                else
                    %disp("3.2 (si = sj) since rand > qRjRi the individual dies");
                    % since rand > qRjRi the individual dies
                    species(selectedS) = species(selectedS)-1;
                end
            else
                %disp("4 (si = sj) less than two metabolites, death");
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
        
%         if(mod(i,dataFreq) == 0)
%             tData(dataIndex) = t +dt;
%             speciesData(:,dataIndex) = species;
%             rData(:,dataIndex)=r;
%             specSum(dataIndex) = sum(species);
%             dataIndex = dataIndex +1;
%         end
        
        t = t+dt;
        i = i+1;
        %w = input('press Enter');
        
    end
end

disp("lapsed runtime: " + toc);
%% postprocessing

% for x=1:ceil
%     if(death(x) ~= -1)
%         death_postp(x) = death(x);
%     end
% end
%% plots
% figure(1)
% specTransp = speciesData';
% area(tData(1,1:1:end),speciesData(:,1:1:end)')
% xlabel('time');
% ylabel('nr individuals');
% legend('species 11','species 21','species 22')
% 
% figure(2)
% plot(tData,specSum)
% xlabel('time');
% ylabel('# individuals');
% 
% figure(3)
% plot(tData,rData)
% xlabel('time');
% ylabel('resources');
% legend('metabolite 1','mteabolite 2')