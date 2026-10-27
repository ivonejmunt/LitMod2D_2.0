import numpy as np

# In Python 3, input() returns a string, so it must be explicitly cast to float for calculations
length = float(input('Enter Profile length:  '))
reso = float(input('Enter resolution along length: '))
data = np.loadtxt("data_phase_vel")

# Ensure the number of nodes used for the loop is an integer
x_nodes = int(length / reso)
prog = 'SURF96'

# Process string input directly, accommodating the user pressing Enter (default)
control_wave = input("What waves do you want: Love ( 0 ) or Rayleigh ( 1 - default value ) : ").strip()
if control_wave == '0':
    surf_wave = 'L'
else:
    surf_wave = 'R'

control_vel = input("Do you want phase ( 0 ) or group ( 1 - default ) velocities : ").strip()
if control_vel == '0':
    vel_type = 'C'
else:
    vel_type = 'U'

mod = 0
xx = 'T'
error = 0.0
# Standard time periods expected by the inversion
time_period = [4, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 45, 50, 59, 67, 77, 87, 100, 111, 125, 143, 167]

for i in range(x_nodes):
    # Safety check: If x_nodes is greater than the number of data rows, terminate early to prevent IndexError
    if i >= len(data):
        print(f"Reached end of data array at node {i}. Stopping file generation.")
        break

    # The correct way to check for NaN in Python 3 using NumPy
    if np.isnan(data[i, 0]):
        pass 
    elif data[i, 0] > 0:
        s = str(i) + "_SURF96.inp"
        
        # Safely manage file writing using the 'with open' context manager
        with open(s, "w") as f:
            for j in range(len(time_period)):
                # Extra protection: Ensure the time period index j does not exceed the number of columns in data
                if j < len(data[i]):
                    f.write("%s \t %s \t %s \t %s \t %d \t %8.4f \t %f \t %f \n" % 
                            (prog, surf_wave, vel_type, xx, mod, time_period[j], data[i, j], error))
    else:
        continue