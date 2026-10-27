#### ploting liberary
import numpy as np
import os
import copy
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.widgets import MultiCursor
from matplotlib import gridspec

matplotlib.rcParams.update({'errorbar.capsize': 2})
origin = 'lower'
label_size = 4
marker_size = 2
matplotlib.rcParams['xtick.labelsize'] = label_size 

def smoothListGaussian(lst, degree=5):  
     window = degree * 2 - 1  
     weight = np.array([1.0] * window)  
     weightGauss = []  
     for i in range(window):  
         i = i - degree + 1  
         frac = i / float(window)  
         gauss = 1 / (np.exp((4 * (frac))**2))  
         weightGauss.append(gauss)  
     weight = np.array(weightGauss) * weight  
     smoothed = [0.0] * (len(lst) - window)  
     for i in range(len(smoothed)):  
         smoothed[i] = sum(np.array(lst[i:i+window]) * weight) / sum(weight)  
     return smoothed  

def plot_lib(x_nodes, l, b, topo_in_file, bouger_in_file, geoid_in_file, FA_in_file, heatflux_in_file, anomaly, c, l_nodes, b_nodes, mat, max_depth,
  label_size, marker_size, anomaly_x, anomaly_y, anomaly_amount, anomaly_type, anomaly_compo, anomaly_color):
  
  try:
    topo = np.loadtxt(topo_in_file)
    max_depth = -max_depth
    min_depth = max(topo[:,1]/1000)
  except Exception:
    max_depth = -400
    min_depth = 3
  
  # Initializing the figure
  fig  = plt.figure()
  # Set the figure window title using the Matplotlib 3 canvas manager.
  fig.canvas.manager.set_window_title('LitMod - Profile Plot')
  gs = gridspec.GridSpec(6, 1, height_ratios=[1,1,1,1,1,3])
  
  # ploting the bodies geometry
  ax1 = plt.subplot(gs[5]) 
  for i in range(len(l_nodes)):
                ax1.plot(l_nodes[i], b_nodes[i], color='black', lw=1, marker=None)
                poly = Polygon(list(zip(l_nodes[i], b_nodes[i])), facecolor=c[i], label=mat[i], lw=1)
                ax1.add_patch(poly)
  for i in range(len(anomaly_x)):
    ax1.plot(anomaly_x[i], anomaly_y[i], color="green", lw=1, marker='o', markersize=4)
    poly = Polygon(list(zip(anomaly_x[i], anomaly_y[i])), facecolor=anomaly_color[i], label="anomaly", lw=1)
    l_text = min(anomaly_x[i]) + (max(anomaly_x[i]) - min(anomaly_x[i])) / 2.0
    b_text = min(anomaly_y[i]) + (max(anomaly_y[i]) - min(anomaly_y[i])) / 2.0
    ax1.add_patch(poly)
    ax1.text(l_text, b_text, str(int(i+1)) + " " + str(anomaly_type[i]) + " Material"
      + str(anomaly_compo[i]) + "" + str(anomaly_amount[i]), fontsize=8, fontstyle='italic', color='black', bbox=dict(facecolor='white', alpha=0.9))
  
  ax1.grid(True)
  ax1.tick_params(labelsize=10)
  ax1.set_ylim((max_depth, min_depth))
  ax1.set_title('Profile Geometry')
  ax1.set_xlabel('Distance ($km$)')
  ax1.set_ylabel('Depth ($km$)')
  
  # plotting the topography
  ax2 = plt.subplot(gs[4])
  # topography data
  f_out = open("topo_out.dat", "r")
  data_out = f_out.readlines()
  try:
    f_in = open(topo_in_file, "r") 
    data_in = f_in.readlines()
    f_in.close()
  except Exception:
    data_in = data_out
    
  profile_out = []
  profile_in = []
  
  if anomaly == 1:
    topo_in = []
    topo_error = []
    topo_out_coupled = []
    topo_out_uncoupled = []
    f_out.close()
    for i in range(len(data_in)):
      try:
        data_temp = data_out[i].split()
      except IndexError:
        data_temp = []
      try:
        if data_temp:
            profile_out.append(float(str(data_temp[0])))
            topo_out_coupled.append(float(str(data_temp[1])))
            topo_out_uncoupled.append(float(str(data_temp[2])))
      except ValueError:
        data_temp = []
      try:
        data_temp = data_in[i].split()
      except IndexError:
        data_temp = []
      try:
        if data_temp:
            profile_in.append(float(str(data_temp[0])))
            topo_in.append(float(str(data_temp[1])))
            try:
              topo_error.append(float(str(data_temp[2])))
            except:
              topo_error.append(0)
      except ValueError:
        data_temp = []
        
    ax2.errorbar(profile_in, topo_in, yerr=topo_error, fmt='o-', markersize=marker_size, ecolor='b', label='Observed')
    ax2.plot(profile_out, topo_out_coupled, 'k-', label='Coupled')
    ax2.plot(profile_out, topo_out_uncoupled, 'r-', label='Uncoupled')
    ax2.grid(True)
    ax2.set_title('Elevation ')
    ax2.xaxis.set_ticklabels([])
    ax2.tick_params(labelsize=10)
    ax2.set_ylabel('Elevation ($m$)')
    ax2.legend(fancybox=True, framealpha=0.5, fontsize=8, loc='best')
  else:
    topo_in = []
    topo_error = []
    topo_out = []
    f_out.close()
    for i in range(len(data_in)):
      try:
        data_temp = data_out[i].split()
      except IndexError:
        pass
      try:
        profile_out.append(float(str(data_temp[0])))
        topo_out.append(float(str(data_temp[1])))
      except ValueError:
        pass
      try:
        data_temp = data_in[i].split()
      except IndexError:
        pass
      try:
        profile_in.append(float(str(data_temp[0])))
        topo_in.append(float(str(data_temp[1])))
        try:
          topo_error.append(float(str(data_temp[2])))
        except:
          topo_error.append(0)
      except ValueError:
        pass    
       
    ax2.errorbar(profile_in, topo_in, yerr=np.abs(topo_error), fmt='o-', markersize=marker_size, ecolor='b', label='Observed')    
    ax2.plot(profile_out, topo_out, 'r-', label='Calculated')
    ax2.grid(True)
    ax2.set_title('Elevation')
    ax2.xaxis.set_ticklabels([])
    ax2.tick_params(labelsize=10)
    ax2.set_ylabel('Elevation ($m$)')
    ax2.legend(fancybox=True, framealpha=0.5, fontsize=8, loc='best')

  # bouguer data
  f_out = open("bouguer_out.dat", "r")
  data_out = f_out.readlines()
  try:
    f_in = open(bouger_in_file, "r") 
    data_in = f_in.readlines()
    f_in.close()
  except Exception:
    data_in = data_out
  profile_out = []
  profile_in = []
  boug_in = []
  boug_out = []
  boug_error = []
  for i in range(len(data_in)):
    try:
      data_temp = data_out[i].split()
    except IndexError:
      pass
    try:
      profile_out.append(float(str(data_temp[0])))
      boug_out.append(float(str(data_temp[1])))
    except ValueError:
      pass
    try:
      data_temp = data_in[i].split()
    except IndexError:
      pass
    try:
      profile_in.append(float(str(data_temp[0])))
      boug_in.append(float(str(data_temp[1])))
      try:
          boug_error.append(float(str(data_temp[2])))
      except:
          boug_error.append(0.0)
    except ValueError:
      pass  

  ax3 = plt.subplot(gs[3])
  ax3.errorbar(profile_in, boug_in, yerr=np.abs(boug_error), fmt='o-', markersize=marker_size, ecolor='b', label='Observed')
  ax3.plot(profile_out, boug_out, 'r-', label='Calculated')
  ax3.grid(True)
  ax3.set_title('Bouguer')
  ax3.xaxis.set_ticklabels([])
  ax3.tick_params(labelsize=10)
  ax3.legend(fancybox=True, framealpha=0.5, fontsize=8, loc='best')
  ax3.set_ylabel('Bouguer ($mGal$)')
  f_out.close()
  
  # geoid data
  f_out = open("geoid_out.dat", "r")
  data_out = f_out.readlines()
  try:
    f_in = open(geoid_in_file, "r")
    data_in = f_in.readlines()
    f_in.close()
  except Exception:
    data_in = data_out
  profile_out = []
  profile_in = []
  geoid_in = []
  geoid_out = []
  geoid_error = []
  for i in range(len(data_in)):
    try:
      data_temp = data_out[i].split()
    except IndexError:
      pass
    try:
      profile_out.append(float(str(data_temp[0])))
      geoid_out.append(float(str(data_temp[1])))
    except ValueError:
      pass
    try:
      data_temp = data_in[i].split()
    except IndexError:
      pass
    try:
      profile_in.append(float(str(data_temp[0])))
      geoid_in.append(float(str(data_temp[1])))
      try:
        geoid_error.append(float(str(data_temp[2])))
      except IndexError:
        geoid_error.append(0.0)        
    except ValueError:
      pass  

  ax4 = plt.subplot(gs[2])
  # ax4.plot(profile_in, geoid_in, 'b-')
  ax4.errorbar(profile_in, geoid_in, yerr=np.abs(geoid_error), fmt='o-', markersize=marker_size, ecolor='b', label='Observed')
  ax4.plot(profile_out, geoid_out, 'r-', label='Calculated')
  ax4.grid(True)
  ax4.set_title('Geoid')
  ax4.xaxis.set_ticklabels([])
  ax4.tick_params(labelsize=10)
  ax4.legend(fancybox=True, framealpha=0.5, fontsize=8, loc='best')
  ax4.set_ylabel('Geoid ($m$)')
  geoid_limits = []
  if geoid_in:
    geoid_limits.extend((np.array(geoid_in) - np.abs(geoid_error)).tolist())
    geoid_limits.extend((np.array(geoid_in) + np.abs(geoid_error)).tolist())
  if geoid_out:
    geoid_limits.extend(geoid_out)
  geoid_limits = [value for value in geoid_limits if np.isfinite(value)]
  if geoid_limits:
    geoid_min = min(geoid_limits)
    geoid_max = max(geoid_limits)
    if geoid_min == geoid_max:
      geoid_pad = max(abs(geoid_min) * 0.05, 1.0)
    else:
      geoid_pad = (geoid_max - geoid_min) * 0.08
    ax4.set_ylim((geoid_min - geoid_pad, geoid_max + geoid_pad))
  
  f_out.close()

  # heat flux data
  f_out = open("SHF_out.dat", "r")
  data_out = f_out.readlines()
  f_out.close()
  try:
    f_in = open(heatflux_in_file, "r")
    data_in = f_in.readlines()
    f_in.close()
  except Exception:
    data_in = data_out
  profile_out = []
  profile_in = []
  flux_in = []
  flux_out = []
  flux_error = []
  for i in range(len(data_out)):
    try:
      data_temp = data_out[i].split()
    except IndexError:
      pass
    try:
      profile_out.append(float(str(data_temp[0])))
      flux_out.append(float(str(data_temp[1])))
    except ValueError:
      pass
  for i in range(len(data_in)):
    try:
      data_temp = data_in[i].split()
    except IndexError:
      pass
    try:
      if data_temp:
          profile_in.append(float(str(data_temp[0])))
          flux_in.append(float(str(data_temp[1])))
          try:
            flux_error.append(float(str(data_temp[2])))
          except IndexError:
            flux_error.append(0.0)
    except ValueError:
      pass

  ax5 = plt.subplot(gs[0])
  ax5.errorbar(profile_in, flux_in, yerr=np.abs(flux_error), fmt='o-', markersize=marker_size, ecolor='b', label='Observed')
  ax5.plot(profile_out, flux_out, 'r-', label='Calculated')
  ax5.grid(True)
  ax5.set_title('Heat Flux')
  ax5.xaxis.set_ticklabels([])
  ax5.tick_params(labelsize=10)
  ax5.legend(fancybox=True, framealpha=0.5, fontsize=8, loc='best')
  ax5.set_ylabel('Heat-flux ($mW/m^2$)')

  # FA data
  f_out = open("FA_out.dat", "r")
  data_out = f_out.readlines()
  f_out.close()
  try:
    f_in = open(FA_in_file, "r")
    data_in = f_in.readlines()
    f_in.close()
  except Exception:
    data_in = data_out
  profile_out = []
  profile_in = []
  FA_in = []
  FA_out = []
  FA_error = []
  for i in range(len(data_in)):
    try:
      data_temp = data_out[i].split()
    except IndexError:
      pass
    try:
      profile_out.append(float(str(data_temp[0])))
      FA_out.append(float(str(data_temp[1])))
    except ValueError:
      pass
    try:
      data_temp = data_in[i].split()
    except IndexError:
      pass
    try:
      if data_temp:
          profile_in.append(float(str(data_temp[0])))
          FA_in.append(float(str(data_temp[1])))
          try:
            FA_error.append(float(str(data_temp[2])))
          except IndexError:
            FA_error.append(0.0)
    except ValueError:
      pass  
      
  ax6 = plt.subplot(gs[1])
  ax6.errorbar(profile_in, FA_in, yerr=np.abs(FA_error), fmt='o-', markersize=marker_size, ecolor='b', label='Observed')
  ax6.plot(profile_out, FA_out, 'r-', label='Calculated')
  ax6.grid(True)
  ax6.set_title('Free-Air')
  ax6.xaxis.set_ticklabels([])
  ax6.tick_params(labelsize=10)
  ax6.legend(fancybox=True, framealpha=0.5, fontsize=8, loc='best')
  ax6.set_ylabel('FA (mGal)')

  multi = MultiCursor(fig.canvas, (ax1, ax2, ax3, ax4, ax5, ax6), color='r', lw=1)

  ##############33 temperature profile plotting
  gs_ = gridspec.GridSpec(4, 1)
  fig_out  = plt.figure()
  # Set the result window title using the Matplotlib 3 canvas manager.
  fig_out.canvas.manager.set_window_title('LitMod - Temp & Dens & Velocities')
  ax_temp = plt.subplot(gs_[0]) 
  f = open("tempout.dat", "r")
  data_out = f.readlines()
  temp_x = []
  temp_y = []
  temp_z = []
  for i in range(2, len(data_out)):
    data_temp = data_out[i].split()
    temp_x.append(float(str(data_temp[0])))
    temp_y.append(float(str(data_temp[1])))
    temp_z.append(float(str(data_temp[2])))
  xlist = []
  ylist = np.array(temp_y[0:96])
  z = temp_z 
  Z = []
  temp = 0
  for i in range(0, x_nodes):
    tmp = []
    xlist.append(temp_x[temp])
    for j in range(0, 96):
      tmp.append(z[temp])
      temp = temp + 1
    Z.append(tmp)
  X, Y = np.meshgrid(ylist, xlist)   
  
  for i in range(len(l)):
    ax_temp.plot(l[i][:], b[i][:], color='grey', linewidth=2)
  levels = [300, 500, 600, 1000, 1200, 1300, 1450, 1550, 1650]
  
  # -------------------------------------------------------------------------
  # Prefer with_extremes to avoid Matplotlib deprecation warnings.
  # -------------------------------------------------------------------------
  try:
      cmap_coolwarm = plt.get_cmap('coolwarm').resampled(20).with_extremes(under='white', over='black')
  except AttributeError:
      # Fall back for older Matplotlib versions.
      cmap_coolwarm = copy.copy(plt.get_cmap('coolwarm', 20))
      cmap_coolwarm.set_under('white')
      cmap_coolwarm.set_over('black')
  
  CS3 = ax_temp.contourf(Y, X, Z, levels, extend='both', cmap=cmap_coolwarm)
  CS4 = ax_temp.contour(Y, X, Z, levels, linewidths=0.2, linestyles='dashed', colors='black')
  ax_temp.clabel(CS4, inline=True, inline_spacing=10, fontsize=10, fmt='%1.0f', colors='black')
  plt.colorbar(CS3, ax=ax_temp, label='$^oC$')
  ax_temp.set_title('Temperature', fontsize=10, fontweight='bold')
  ax_temp.tick_params(labelsize=10)
  ax_temp.set_ylabel('Depth ($km$)')

  ##############33 density profile plotting
  f = open("dens_node2.dat", "r")
  data_out = f.readlines()
  temp_x = []
  temp_y = []
  temp_z = []
  for i in range(len(data_out)):
    data_temp = data_out[i].split()
    temp_x.append(float(str(data_temp[0])))
    temp_y.append(float(str(data_temp[1])))
    temp_z.append(float(str(data_temp[2])))
  xlist = []
  ylist = np.array(temp_y[0:95])
  z = temp_z
  Z = []
  temp = 0
  for i in range(0, x_nodes):
    tmp = []
    xlist.append(temp_x[temp])
    for j in range(0, 95):
      tmp.append(z[temp])
      temp = temp + 1
    Z.append(tmp)
  X, Y = np.meshgrid(ylist, xlist)   
  
  ax_dens = plt.subplot(gs_[1])
  for i in range(len(l)):
    ax_dens.plot(l[i][:], b[i][:], color='grey', linewidth=2)
  levels_contourf = np.arange(2800, 3800, 50)
  
  # Matplotlib 3 compatibility handling.
  try:
      cmap_rainbow = plt.get_cmap('rainbow').resampled(20).with_extremes(under='white', over='black')
  except AttributeError:
      cmap_rainbow = copy.copy(plt.get_cmap('rainbow', 20))
      cmap_rainbow.set_under('white')
      cmap_rainbow.set_over('black')
  
  CS3 = ax_dens.contourf(Y, X, Z, levels_contourf, extend='both', cmap=cmap_rainbow, origin=origin)
  levels_contour = np.arange(3000, 3800, 50)
  CS4 = ax_dens.contour(Y, X, Z, levels_contour, linewidths=0.2, linestyles='dashed', colors='grey')
  ax_dens.clabel(CS4, inline=True, inline_spacing=2, fontsize=10, fmt='%1.0f', colors='black')
  plt.colorbar(CS3, ax=ax_dens, label='$kg/m^3$')
  ax_dens.set_title('Density', fontsize=10, fontweight='bold')
  ax_dens.tick_params(labelsize=10)
  ax_dens.set_ylabel('Depth ($km$)')

  ##############33 Vp profile plotting
  f = open("velocities.dat", "r")
  data_out = f.readlines()
  temp_x = []
  temp_y = []
  temp_z_vp = []
  temp_z_vs = []
  for i in range(len(data_out)):
    data_temp = data_out[i].split()
    temp_x.append(float(str(data_temp[0])))
    temp_y.append(float(str(data_temp[1])))
    temp_z_vp.append(float(str(data_temp[2])))
    temp_z_vs.append(float(str(data_temp[3])))
  xlist = []
  ylist = np.array(temp_y[0:95])
  z_vp = temp_z_vp 
  z_vs = temp_z_vs
  Z_vp = []
  Z_vs = []
  temp = 0
  for i in range(0, x_nodes):
    tmp1 = []
    tmp2 = []
    xlist.append(temp_x[temp+1])
    for j in range(0, 95):
      tmp1.append(z_vp[temp])
      tmp2.append(z_vs[temp])
      temp = temp + 1
    Z_vp.append(tmp1)
    Z_vs.append(tmp2)
  X, Y = np.meshgrid(ylist, xlist)   
  
  ax_vp = plt.subplot(gs_[2])
  for i in range(len(l)):
    ax_vp.plot(l[i][:], b[i][:], color='grey', linewidth=2)
  levels_contourf_vp = np.arange(7.2, 8.8, 0.01)
  
  # Matplotlib 3 compatibility handling.
  try:
      cmap_rdbu_vp = plt.get_cmap('RdYlBu').resampled(20).with_extremes(under='white', over='black')
  except AttributeError:
      cmap_rdbu_vp = copy.copy(plt.get_cmap('RdYlBu', 20))
      cmap_rdbu_vp.set_under('white')
      cmap_rdbu_vp.set_over('black')
  
  CS3 = ax_vp.contourf(Y, X, Z_vp, levels_contourf_vp, extend='both', cmap=cmap_rdbu_vp, origin=origin)
  levels_contour_vp = np.arange(7.7, 8.8, 0.05)
  CS4 = ax_vp.contour(Y, X, Z_vp, levels_contour_vp, linewidths=0.2, linestyles='dashed', colors='black')
  ax_vp.clabel(CS4, inline=True, inline_spacing=10, fontsize=10, fmt='%1.1f', colors='black')
  plt.colorbar(CS3, ax=ax_vp, label='$km/s$')
  ax_vp.tick_params(labelsize=10)
  ax_vp.set_title('P wave velocity', fontsize=10, fontweight='bold')
  ax_vp.set_ylabel('Depth ($km$)') 

  ax_vs = plt.subplot(gs_[3])
  for i in range(len(l)):
    ax_vs.plot(l[i][:], b[i][:], color='grey', linewidth=2)
  levels_contourf_vs = np.arange(4.2, 4.8, 0.01)
  
  # Matplotlib 3 compatibility handling.
  try:
      cmap_rdbu_vs = plt.get_cmap('RdYlBu').resampled(20).with_extremes(under='white', over='black')
  except AttributeError:
      cmap_rdbu_vs = copy.copy(plt.get_cmap('RdYlBu', 20))
      cmap_rdbu_vs.set_under('white')
      cmap_rdbu_vs.set_over('black')
  
  CS3 = ax_vs.contourf(Y, X, Z_vs, levels_contourf_vs, extend='both', cmap=cmap_rdbu_vs, origin=origin)
  levels_contour_vs = np.arange(4.3, 5.1, 0.05)
  CS4 = ax_vs.contour(Y, X, Z_vs, levels_contour_vs, linewidths=0.2, linestyles='dashed', colors='black')
  ax_vs.clabel(CS4, inline=True, inline_spacing=10, fontsize=10, fmt='%1.1f', colors='black')
  plt.colorbar(CS3, ax=ax_vs, label='$km/s$')
  ax_vs.set_title('S wave velocity', fontsize=10, fontweight='bold')
  ax_vs.set_xlabel('Distance ($km$)')
  ax_vs.set_ylabel('Depth ($km$)')
  ax_vs.tick_params(labelsize=10)
  
  plt.savefig('Temp_Dens_velocities.png', dpi=300) 
  plt.show()
