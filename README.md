# ClusterDetectionDEM

ClusterDetectionDEM  detects clusters in gas-solid flows based on particle diameter dp and their coordinates (x,y,z). If you use this tool, please reference the following publication: 

Bakshi, A., Altantzis, C., Bates, R.B., and Ghoniem, A. F.,"Multiphase-flow Statistics using 3D Detection and Tracking Algorithm (MS3DATA): Methodology and application to large-scale fluidized beds.", Chemical Engineering Journal 293 (2016): 355-364


The following algorithm is implemented: 

1. currentlist = [ ], pendinglist = list of all particles, clusterindex = 0  

2. "first pass": for iparticle in pendinglist
      - particlesincontact = list of all particles (absolute indiced) in contact with iparticle
      - clusterindex += 1
      - associate iparticle and intersection(particlesincontact, pendinglist) with clusterindex 
      - append currentlist with iparticle, particlesincontact
      - remove iparticle, particlesincontact from pendinglist   
3. "second pass": for iparticle in list of all particles skipped in "first pass" 
      - particlesincontact = list of all particles (absolute indiced) in contact with iparticle
      - equivalentclusters = list of cluster indices associated with particlesincontact 
      - set cluster indices of all particles associated with equivalentclusters to min(equivalentclusters) 
4. compute properties of clusters by aggregating particles: [# particles, x-centroid, y-centroid, z-centroid]

Notes: 
1. Current implementation assumes all particles of uniform diameter
2. To relax criteria for cluster detection (esp. if solids loading is dilute), use dp = (1+tolerance) x particle-diameter 
3. Python implementation is ~2x slower, probably because of nested if within for loop 
4. Scalability of code has not been tested. Suggestions welcome! 

