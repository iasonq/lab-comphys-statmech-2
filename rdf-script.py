# %%
import numpy as np

# %%
#we need to read the file
#../input_data/trajNaCl.pdb
#500 timeframes

#bash command to filter
#grep ATOM trajNaCl.pdb | tr -s ' ' | cut -d' ' -f6-8 > filter_trajNaCl.txt

# %%
data = np.loadtxt("../input_data/filter_trajNaCl.txt") #500 blocks of 40 lines
#Distances in Angstrom
cell = np.zeros(3)+31.47
#Since the box is square, cell 1-dimensional based on the check needed to do leater?

# %%
#Initialize data
bins = ? #total number of bins
ngr = 0 # what is ngr for ?
delg = cell/(2*bins) #delg is bin size... ?
g = np.zeros(bins) 
#First chunk over
npart = 40 #amount of particles

# %%
#I have to figure out how the periodic conditions are given, is box a vector or not?
def rdf(g, ngr, cell, data, delg):
    ngr = ngr + 1 #should be the right place for ngr, it seems to count timeframes?
    for k in range(0, 500, npart): #k for the different timeframes
        #the actual pairing
        #do it seperately for Na-Na, Na-Cl, Cl-Cl? Seems all at once is the right way
        for i in range(k+0, k+npart-1):
            for j in range(k+1, k+npart):
                #distance between a particle pair
                r = data[i]-data[j]
                #implement periodic boundary conditions
                r = r - cell*np.round(r/cell)
                r = np.linalg.norm(r)
                #half-box check & contribution
                if r <= cell/2: #this needs to be true in all fields to proceed right?
                    l = int(r/delg)
                    g(l) = g(l) + 2
    return g, ngr

# %%
def grdf(g, ngr):
    rho = ? #Density of water+NaCl or what?
    for i in range(0, ngr):
        r = delg*(i+0.5) #what is this supposed to be?
        vol = ((i+1)**3  - i**3 )*delg**3 #why is this a volume and then its used in the 4/3 pi r^3 rho formula, weird
        nid = (4/3)*np.pi*vol*rho  #according to the description a density
        g(i) = g(i)/(ngr*npart*nid) #ngr is the pair number? npart?  
    return g


