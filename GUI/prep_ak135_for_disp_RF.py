import numpy as np 

# Read the standard ak135f Earth model, skipping the first 3 header lines
# In Python 3, np.loadtxt handles this natively and perfectly.
data = np.loadtxt("ak135f.txt", skiprows=3)

s = "modl.d"

# Safely manage the file handle using the 'with open' context manager
with open(s, "w") as f:
    # Write the strict header required by Computer Programs in Seismology (surf96)
    f.write("MODEL\n\
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

    for j in range(len(data) - 1):
        # This script acts as a "code generator". 
        # It formats the AK135 data into literal Python 'f.write("...")' statements,
        # which can then be copy-pasted into the mantle section of your other scripts.
        # [FIX] The string formatting has been corrected to generate perfectly valid 
        # Python code, complete with the closing quotation marks and newlines.
        f.write("f.write(\"%-10.6f \\t %-10.6f \\t %-10.6f \\t %-10.6f \\t %-10.6f \\t %-10.6f \\t 0.0 \\t 0.0 \\t 1.0 \\t 1.0 \\n\")\n"
             % (data[j+1, 0] - data[j, 0], data[j, 2], data[j, 3], data[j, 1], data[j, 4], data[j, 5]))