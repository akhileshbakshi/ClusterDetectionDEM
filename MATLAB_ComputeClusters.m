%% ----------------------------------------------------------------------
% this code detects clusters based on particle size dp and location data 
% code output: clusterdistribution which is Ncx4 matrix 
% where number of rows Nc = number of clusters detected 
% and each row has format [# particles with cluster, x-centroid, y-centroid, z-centroid]
% NOTE: scaling of code has not been tested 
clear all; clc;

%% user variables 
dp = 0.5;                                                   % distance between particle centroids for cluster detection. default = [particle diameter] 
particlelocation = load("currentparticlelocation.txt");     % particle location = [x, y, z]

%% initializing variables
numparticles = size(particlelocation,1);
neighbormatrix = [linspace(1,numparticles,numparticles)', zeros(numparticles, 1), zeros(numparticles, 1)];
% neighbormatrix = [particle-ID, cluster-ID, considered in first pass?]
currentlist = [];
pendinglist = linspace(1,numparticles,numparticles)';
clusterindex = 0;

%% first pass: link every particle with its neighbors in pendinglist using clusterindex 
while ~isempty(pendinglist)
    iparticle = pendinglist(1); 
    neighbormatrix(iparticle,3) = 1;
    currentparticlelocation = particlelocation(iparticle,:);
    particlesincontact = find(sqrt(sum((currentparticlelocation-particlelocation).^2,2))<dp);
   % particlesincontact = particles (absolute indices) in contact with iparticle, including iparticle itself 
    if length(particlesincontact)==1
             pendinglist = setdiff(pendinglist,iparticle);
    else  
        clusterindex = clusterindex + 1;    
        neighbormatrix(intersect(pendinglist,particlesincontact),2) = clusterindex;
        currentlist = [currentlist; intersect(pendinglist,particlesincontact)];
        pendinglist = setdiff(pendinglist,intersect(pendinglist,particlesincontact));  
    end
end


%% second pass: check clusterindices of neighbors of all particles skipped in first pass 
particleIDs = linspace(1,numparticles,numparticles); 
for iparticle = particleIDs(neighbormatrix(:,3)==0)   
    currentparticlelocation = particlelocation(iparticle,:);
    particlesincontact = find(sqrt(sum((currentparticlelocation-particlelocation).^2,2))<dp);
    equivalentclusters = unique(neighbormatrix(particlesincontact,2));
    if length(equivalentclusters)>1
        neighbormatrix(ismember(neighbormatrix(:,2),equivalentclusters),2) = min(neighbormatrix(particlesincontact,2)); 
    end
end


%% compute distribution of clusters in format: [num-particles, x, y, z] 
clusterdistribution = [];
uniqueclusters = unique(neighbormatrix(:,2));
uniqueclusters = uniqueclusters(2:end);
for iclusters = 1:length(uniqueclusters)
    TF = neighbormatrix(:,2) == uniqueclusters(iclusters);
    centroid = mean(particlelocation(neighbormatrix(TF,1),:));
    clusterdistribution = [clusterdistribution; sum(TF), centroid];
end



