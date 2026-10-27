import numpy as np 

# In Python 3, input() returns a string, so it must be explicitly cast to float
length = float(input('Enter Profile length:  '))
reso = float(input('Enter resolution along length:  '))
data = np.loadtxt('post_processing_output.dat', usecols=(0, 1, 2, 3, 4, 5, 6))
vp_vs = 1.73

### Data format:
#  X    Z   T   P   Vp  Vs  density     Material mk??

## Calculating Vp in the crust according to Brocher (2005)
f = 0
for i in range(len(data)):
    if data[i, 5] == 0:
        dens = data[i, 6] / 1000.0
        
        # Brocher (2005) empirical relation (Eq. 1) for Vp based on density
        data[i, 4] = 39.128 * dens - 63.064 * dens**2 + 37.083 * dens**3 - 9.1819 * dens**4 + 0.8215 * dens**5
        
        # Calculate Vs using the constant Vp/Vs ratio
        data[i, 5] = data[i, 4] / vp_vs
        
        ### According to equation 6 Brocher BSSA 2005
        # data[i,5] = 0.7858 - 1.2344 * data[i,4] + 0.7949 * data[i,4]**2 - 0.1238 * data[i,4]**3 + 0.0064 * data[i,4]**4
        
        ### According to equation 7 Brocher BSSA 2005
        # data[i,5] = (data[i,4] - 1.36) / 1.16 # does not work

        ### According to equation 8 Brocher BSSA 2005
        # data[i,5] = 2.88 + 0.52 * (data[i,4] - 5.25)  # works but not perfect

    else:
        pass

# Safely manage the file handle using the 'with open' context manager
with open("velocities_crust_mantle.dat", "w") as f_out:
    for i in range(len(data)):
        if i == len(data) - 1:
            f_out.write(" %f  %f  %f  %f  %f  " % (data[i, 0], -data[i, 1], data[i, 6] / 1000.0, data[i, 4], data[i, 5]))
        else:
            f_out.write(" %f  %f  %f  %f  %f  \n" % (data[i, 0], -data[i, 1], data[i, 6] / 1000.0, data[i, 4], data[i, 5])) 

# Reload the newly formatted dataset (which now has 5 columns instead of 7)
data = np.loadtxt("velocities_crust_mantle.dat")

# In Python 3, ensure the argument used for range() is cast to an integer
x_nodes = int(length / reso)

for i in range(x_nodes):
    s = str(i * reso) + "_vel.dat"
    with open(s, "w") as f_out:
        f_out.write("MODEL\n\
TEST MODEL \n\
ISOTROPIC \n\
KGS \n\
FLAT EARTH \n\
1-D \n\
CONSTANT VELOCITY \n\
LINE08 \n\
LINE09 \n\
LINE10 \n\
LINE11 \n\
HR        VP       VS        RHO    QP    QS    ETAP  ETAS FREFP FREFS \n")
        
        # Iterating through blocks of 94 nodes per vertical profile
        for j in range(i + (94 * i), i + (94 * i) + 94):
            f_out.write("%f    \t  %f    \t  %f   \t  %f   \t   %s   \t  %s  \t %s \t  %s \t  %s  \t %s \n"
             % (data[j+1, 1] - data[j, 1],  data[j, 3],   data[j, 4],  data[j, 2], "200.0",  "100", "0.0",  "0.0",  "1.0", "1.0"))
             
        # Write the data for the final boundary node
        f_out.write("%f    \t  %f    \t  %f   \t  %f   \t   %s   \t  %s  \t %s \t  %s \t  %s  \t %s \n"
             % (data[-1, 1] - data[-2, 1],  data[-1, 3],   data[-1, 4],  data[-1, 2], "200.0",  "100", "0.0",  "0.0",  "1.0", "1.0"))