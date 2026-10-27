# -*- coding: utf-8 -*-
#####################################################################################################
# This is the library which contains functions used in the main.py - Updated for Python 3.12 & PyQt6
#####################################################################################################
import os
import sys, string, platform
from shutil import copyfile
import math
from datetime import datetime

# ==================================================================================
# Core Modification: Force TraitsUI and Matplotlib to use PyQt5 backend, completely abandoning Tkinter!
# ==================================================================================
os.environ['ETS_TOOLKIT'] = 'qt'
os.environ['QT_API'] = 'pyqt5'

import matplotlib
matplotlib.use('Qt5Agg') # Force matplotlib to use Qt engine for rendering
# ==================================================================================

import numpy as np
from shapely.geometry.polygon import  Polygon as Shape
import shapely as shapely
from traits.api import HasTraits, Str, Int, Float, Enum, Array, File, Directory
from traitsui.api import View, Item, Group, HSplit, Handler, Color
from traitsui.tabular_adapter import TabularAdapter
from traitsui.menu import OKButton, CancelButton, ApplyButton, UndoButton
from traits.api import *

import matplotlib.image as image
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Cursor
from matplotlib.figure import Figure
from matplotlib.patches import Polygon

import webbrowser # for web links
plt.ioff()
sys.setrecursionlimit(1500)
######################################################
### importing custom made libraries for the GUI
import plot_lib
import edit_body as edit
import shape_sum as shape_sum
import shape_substract as shape_substract
import shape_split as shape_split
import envelop_maker as envelop_maker
import time

#################################################################################
#### This a class where matplotlib is embedded with plotting area, functions relate to
#### to the GUI. 
class Canvas(HasTraits, object):
    #####################3##########################################################
    #### model setup
    X_S = 0; X_E = 1000; Y_length = 400; X_resolution = 10; dir = Directory
    #### body properties
    name = "mantle"
    color = Enum('green','blue','yellow','red','pink','orange','purple','cyan','magenta','black','white','gold','grey','violet','maroon','aqua','lime','tan','rosybrown')
    color_edi= Color('blue')
    Body_number= "0"
    ref_density = Str("3.245E+03", desc="Reference density (kg/m3)")
    depth_density_factor= Str("0.000E+00", desc="Depth dependence of density (m)")
    Temp_density_factor= Str("0", desc="Temperature dependence of density")
    Thermal_expansion= Str("0.000E+00", desc="Coefficient of thermnal expansion(only for crust)")
    material = Int(desc="0:no thermodynamic table, 1-10: crust with thermodynamic table, >20: lithospheric mantle with Thermo-table")
    ### thermal conductivity parameters
    Kx =Str("2.100E+00", desc="Conductivity Tensor component 1 (W/mK)")
    Ky =Str("2.100E+00", desc="Conductivity Tensor component 2 (W/mK)")
    Kz =Str("0.000E+00", desc="Anisotropy for Conductivity") #angle of anisotropy
    Temp_conduc_factor=Str("0.000E+00", desc="Temperature dependence of thermal conductivity factor")
    Depth_conduc_factor=Str("0.000E+00", desc="Depth dependence of thermal conductivity factor (m)")
    ### heat production parameters
    radiogenic_heat_producion=Str("2.000E-08", desc="Ao of Heat Production (W/m3)")
    depth_heat_factor=Str("0.000E+00", desc="Ref Thickness for exponential Term (m)")
    x_s = "0.0"; y_s = "-100"; x_e = "1000"; y_e = "-100"  
    ## temperature boundary conditions
    T_surf="0"; T_lab="1320"
    #############################################################################3
    #### observed data and some parameteres for input file
    topo    = File
    digitized =File(" ", desc="This is file in which you might have nodes for the bodies you wanat to draw. May be Moho or Lab depth which will be plotted in the background and you can click on them to add you points.")
    bouger  = File(" ", desc="Make sure that sampliing interval and length is exactly same as that of your profile") 
    geoid   = File(" ", desc="Make sure that sampliing interval and length is exactly same as that of your profile") 
    FA      = File(" ", desc="Make sure that sampliing interval and length is exactly same as that of your profile") 
    SHF     = File(" ", desc="Make sure that sampliing interval and length is exactly same as that of your profile") 
    Anomaly_file =File( ) 
    thermo_data = 5; grain_size = 5; oscillation = 75; anomaly = 0; print_option = 1; topo_response= Enum("no","yes"); digitized_response= Enum("no","yes")
    ## for anomaly
    body_type=Enum("normal","anomaly","split")
    z_levels=Enum(1,2,3,4,5)
    anomaly_type=Enum("Vp(%)","Vs(%)","Vp(>6.5)/Vs(<6.5)","Thermal","Compositional","File")
    anomaly_value=Float
    anomaly_mat=Int(99)
    anomaly_loaded_alert= "Hit ok and look at the command line for more info."
    ## for seismic anomaly from a file
    anom_levels=10
    anom_nodes=10
    anomaly_file = File("", desc="Select file with seismic anomalies in a specified format. If no file is selected then LitMod does not look for a file.")
    ## warning message
    body_alert_message="Body is not valid. Please delete this one"
    split_body_alert_message="Could not split this body!! Make sure you have nodes at both end"
    split_body_no=Int(0)
    merge_body_no_1=Int(0)
    merge_body_no_2=Int(0)
    delete_body_no=Int(0)
    shape_alert_message=" Please finish drwaing all bodies first!! Model is finished by clicking the Close Model button at mid upper panel."
    exit_alert=Enum("no","yes")
    body_delete_alert=Enum("no","yes")
    ans=Enum("no","yes")
    close_model_alert="Model is not closed yet! Please close the model first otherwise you will loose the progress :( ."
    plot_output_alert=" "
    load_model_message="You do not have bodies_GUI.out file. It can be generated by running LitMod in the command with a input file. Close GUI and open when you have bodies_GUI.out. "
    ## Plotting part user interactive variables
    max_depth=Float(400, desc='Depth upto which profile will be plotted')
    label_size=Float(3, desc='Size of the data labels')
    marker_size=Float(2, desc='Data points marker size')

    # Views definitions
    model_set_up_view = View(Group(Item(name = 'X_S',label="Start of the profile (km)"),
                 Item(name = 'X_E',label="End of the profile (km)"),
                 Item(name = 'Y_length',label="Depth of the profile (must not be changed)"),
                 Item(name = 'X_resolution',label="Resolution along length of the profile (5,10 km)"),
                 Item(name = 'dir'),
                 show_border = True),
                 title="Model Setup",
                 buttons = [OKButton])
                 
    save_view = View(Group(Item(name = 'topo',label="Topography Input File"),
                 Item(name = 'bouger',label="Bouguer Input File"),
                 Item(name = 'geoid',label="Geoid Input File"),
                 Item(name = 'FA',label="Free-Air Input File"),
                 Item(name = 'SHF',label="Heat-Flow Input File"),
                 Item(name = 'thermo_data',label="Thermodynamic Table(default:5)"),
                 Item(name = 'anomaly_file',label="File with seismic anomalies."),
                 Item(name = 'anomaly',label="Switch for Anomalous Body"),
                 Item(name = 'print_option'),
                 Item(name = 'T_surf',label="Temperature at surface"),
                 Item(name = 'T_lab',label="Temperature at LAB"),
                 show_border = True),
                 title="Input Files and some control inputs",
                 buttons = [OKButton])   
                 
    add_view = View(Group(Item(name = 'name'),
                 Item(name = 'Body_number',label = "Body Number"),
                 Item(name = 'color_edi',label="Select color"),
                 Item(name = 'material'),
                 Item(name = 'body_type'),
                 show_border = True),
                 title="Property Editor",resizable=True,
                 buttons = [OKButton])
                 
    add_view_no_mat = View(Group(Item(name = 'name'),
                 Item(name = 'Body_number',label = "Body Number"),
                 Item(name = 'color_edi',label="Select color"),
                 Item(name = 'material'),
                 Group(Item(name = 'ref_density',label = "Reference Density(kg/m3)"),Item(name = 'depth_density_factor',label="P dependence factor"),\
                    Item(name = 'Temp_density_factor',label="T dependence factor"),Item(name = 'Thermal_expansion',label="Thermal expansion  factor"),\
                    orientation='horizontal', layout='tabbed', springy=True,show_labels=True, label="Density Parameters"),
                 Group(Item(name = 'Kx',label="Conductivity in X"),Item(name = 'Ky',label="Conductivity in Y"),\
                    Item(name = 'Kz',label="Anisotropy angle"),Item(name = 'Temp_conduc_factor',label="T Conductivity Factor"),\
                    Item(name = 'Depth_conduc_factor',label="Depth Conductivity Factor")
                    ,orientation='horizontal', layout='tabbed', springy=True,show_labels=True,label="Thermal Conductivity"),
                 Group(Item(name = 'radiogenic_heat_producion',label="Heat production A0(W/m3)"),Item(name = 'depth_heat_factor',label="Reference thickness of exponential production(m)"),\
                    orientation='horizontal', layout='tabbed', springy=True,show_labels=True, label="Heat Production  Parameters"),
                 Item(name = 'body_type'),
                 show_border = True),
                 title="Property Editor",resizable=True,
                 buttons = [OKButton])
                 
    add_view_mat = View(Group(Item(name = 'name'),
                 Item(name = 'Body_number',label = "Body Number"),
                 Item(name = 'color_edi',label="Select color"),
                 Item(name = 'material'),
                 Group(Item(name = 'ref_density',label = "Reference Density(kg/m3)"),Item(name = 'depth_density_factor',label="P dependence factor"),\
                    Item(name = 'Temp_density_factor',label="T dependence factor"),Item(name = 'Thermal_expansion',label="Thermal expansion  factor"),\
                    orientation='horizontal', layout='tabbed', springy=True,show_labels=True, label="Density Parameters"),
                 Group(Item(name = 'Kx',label="Conductivity in X"),Item(name = 'Ky',label="Conductivity in Y"),\
                    Item(name = 'Kz',label="Anisotropy angle"),Item(name = 'Temp_conduc_factor',label="T Conductivity Factor"),\
                    Item(name = 'Depth_conduc_factor',label="Depth Conductivity Factor")
                    ,orientation='horizontal', layout='tabbed', springy=True,show_labels=True,label="Thermal Conductivity"),
                 Group(Item(name = 'radiogenic_heat_producion',label="Heat production A0(W/m3)"),Item(name = 'depth_heat_factor',label="Reference thickness of exponential production(m)"),\
                    orientation='horizontal', layout='tabbed', springy=True,show_labels=True, label="Heat Production  Parameters"),
                 Item(name = 'body_type'),
                 show_border = True),
                 title="Property Editor",resizable=True,
                 buttons = [OKButton])
                 
    edit_view = View(Group(Item(name = 'name'),
                 Item(name = 'Body_number'),
                 Item(name = 'color_edi',label="Select Color"),
                 Item(name = 'material'),
                 Group(Item(name = 'ref_density',label = "Reference Density(kg/m3)"),Item(name = 'depth_density_factor',label="P dependence factor"),\
                    Item(name = 'Temp_density_factor',label="T dependence factor"),Item(name = 'Thermal_expansion',label="Thermal expansion  factor"),\
                    orientation='horizontal', layout='tabbed', springy=True,show_labels=True, label="Density Parameters"),
                 Group(Item(name = 'Kx',label="Conductivity in X"),Item(name = 'Ky',label="Conductivity in Y"),\
                    Item(name = 'Kz',label="Anisotropy angle"),Item(name = 'Temp_conduc_factor',label="T Conductivity Factor"),\
                    Item(name = 'Depth_conduc_factor',label="Depth Conductivity Factor")
                    ,orientation='horizontal', layout='tabbed', springy=True,show_labels=True,label="Thermal Conductivity"),
                 Group(Item(name = 'radiogenic_heat_producion',label="Heat production A0(W/m3)"),Item(name = 'depth_heat_factor',label="Reference thickness of exponential production(m)"),\
                    orientation='horizontal', layout='tabbed', springy=True,show_labels=True, label="Heat Production  Parameters"),
                 Item(name = 'body_type'),
                 show_border = True),
                 title="Property Editor",resizable=True,
                 buttons = [OKButton])
                 
    edit_view_mantle = View(Group(Item(name = 'name'),
                 Item(name = 'color_edi',label="Select Color"),
                 Item(name = 'material'),
                 show_border = True),
                 title="Property Editor",resizable=True,
                 buttons = [OKButton])
                 
    topo_view = View(Group(Item(name = 'topo_response', label = "Do you have topography file"),
                 Item(name = 'topo'),
                 show_border = True),
                 title="Topo file",resizable=True,
                 buttons = [OKButton])
                 
    digitized_view = View(Group(Item(name = 'digitized_response', label = "Do you have digitized file for bodies?"),
                 Item(name = 'digitized'),
                 show_border = True),
                 title="Digitized file",resizable=True,
                 buttons = [OKButton])
                 
    body_no_view = View(Group(Item(name = 'body_type'),
                 Item(name = 'Body_number',label="Enter number of the Body to edit"),
                 show_border = True),
                 title="Property Editor",resizable=True,
                 buttons = [OKButton])
                 
    anomaly_select = View(Group(Item(name = 'anomaly_type',label="Choose the type of anomaly"),
                 show_border = True),
                 title="Anomaly",resizable=True,
                 buttons = [OKButton])
                 
    anomaly_view = View(Group(Item(name = 'anomaly_type',label="Choose the type of anomaly"),
                 Item(name = 'anomaly_value',label = "Enter anomaly value"),
                 Item(name = 'anomaly_mat',label = "Enter material file name for anomaly"),   
                 show_border = True),
                 title="Anomaly",resizable=True,
                 buttons = [OKButton])

    anomaly_file_setup_view =  View(Group(Item(name = 'anom_nodes',label = "Enter anomaly value"),
                 Item(name = 'anom_levels',label = "Enter anomaly value")),
                 title='Anomaly file setu', width=0.3, height=0.8, resizable=True, buttons = [OKButton])
                 
    anomaly_file_add_view =  View(Item(name ='anomaly_file',label='File with anomlies:'),
                 title='Anomaly file ', resizable=True, buttons = [OKButton])
    
    anomaly_edit_view=View(Group(Item(name = 'name',label='Anomaly Name'),
                 Item(name = 'anomaly_type'),
                 Item(name = 'color_edi',label="Select Color"),
                 Item(name = 'anomaly_value',label = "Enter anomaly value"),
                 Item(name = 'anomaly_mat',label = "Enter material file name for anomaly"),   
                 show_border = True),
                 title="Anomaly Edit",resizable=True,
                 buttons = [OKButton])
                 
    anomaly_loaded_view = View(Group(Item(name = 'anomaly_loaded_alert',label="You have anomalies in this model."),   
                 show_border = True),
                 title="Anomalous body loaded!",style="readonly",resizable=True,
                 buttons = [OKButton])
                 
    body_alert = View(Group(Item(name = 'body_alert_message',label=""),   
                 show_border = True),
                 title="Alert: Invalid Body!!",style="readonly",resizable=True,
                 buttons = [OKButton])
                 
    split_body_alert = View(Group(Item(name = 'split_body_alert_message',label=""),   
                 show_border = True),
                 title="Alert!!! ",style="readonly",resizable=True,
                 buttons = [OKButton])
    
    shape_alert = View(Group(Item(name = 'shape_alert_message',label=""),   
                 show_border = True),
                 title="Model is not final yet!",style="readonly",resizable=True,
                 buttons = [OKButton])
                 
    load_model_alert = View(Group(Item(name = 'load_model_message',label=""),   
                 show_border = True),
                 title="Can not load model!",style="readonly",resizable=True,
                 buttons = [OKButton])
                 
    body_split = View(Group(Item(name = 'name'),
                 Item(name = 'Body_number'),
                 Item(name = 'color_edi',label="Select color"),
                 Item(name = 'material'),
                 Group(Item(name = 'ref_density',label = "Reference Density(kg/m3)"),Item(name = 'depth_density_factor',label="P dependence factor"),\
                    Item(name = 'Temp_density_factor',label="T dependence factor"),Item(name = 'Thermal_expansion',label="Thermal expansion  factor"),\
                    orientation='horizontal', layout='tabbed', springy=True,show_labels=True, label="Density Parameters"),
                 Group(Item(name = 'Kx',label="Conductivity in X"),Item(name = 'Ky',label="Conductivity in Y"),\
                    Item(name = 'Kz',label="Anisotropy angle"),Item(name = 'Temp_conduc_factor',label="T Conductivity Factor"),\
                    Item(name = 'Depth_conduc_factor',label="Depth Conductivity Factor")
                    ,orientation='horizontal', layout='tabbed', springy=True,show_labels=True,label="Thermal Conductivity"),
                 Group(Item(name = 'radiogenic_heat_producion',label="Heat production A0(W/m3)"),Item(name = 'depth_heat_factor',label="Reference thickness of exponential production(m)"),\
                    orientation='horizontal', layout='tabbed', springy=True,show_labels=True, label="Heat Production  Parameters"),
                 Group(Item(name = 'radiogenic_heat_producion')),
                 Item(name = 'body_type'),
                 show_border = True),
                 title="Property Editor",resizable=True,
                 buttons = [OKButton])
                 
    split_view = View(Group(Item(name = 'ans',label="Are you sure you want to Split:"),
                 Item(name = 'split_body_no',label="Enter number of the body to split:"),   
                 show_border = True),
                 title="Split Body",resizable=True,
                 buttons = [OKButton])
                 
    merge_view = View(Group(Item(name = 'ans',label="Are you sure you want to Merge:"),
                 Item(name = 'merge_body_no_1',label="Enter number of the body to Merge:"),
                 Item(name = 'merge_body_no_2',label="Enter number of the body to Merge (Properties of this body will be retained):"),   
                 show_border = True),
                 title="Merge Body",resizable=True,
                 buttons = [OKButton])   
                 
    exit_view = View(Group(Item(name = 'exit_alert',label="Are you sure you want to exit:"),   
                 show_border = True),
                 title="Exit Model",resizable=True,
                 buttons = [OKButton])
                 
    anomaly_delete_view = View(Group(Item(name = 'delete_body_no',label="Enter No of the anomaly to delete:"),
                 Item(name = 'body_delete_alert',label="Are you sure you want to delete this anomaly:"),                
                 show_border = True),
                 title="Delete Anomaly",resizable=True,
                 buttons = [OKButton])
                 
    body_delete_view = View(Group(Item(name = 'body_delete_alert',label="Are you sure you want to delete current body:"),
                 Item(name = 'body_type'),
                 show_border = True),
                 title="Delete Body",resizable=True,
                 buttons = [OKButton])
                 
    close_model_view = View(Group(Item(name = 'close_model_alert',label="you will loose the progress :("),   
                 show_border = True),
                 title="Alert:  Model not closed yet.",style="readonly",resizable=True,
                 buttons = [OKButton])
                 
    plot_output_view = View(Group(Item(name = 'plot_output_alert',label="ALERT!!!! Something went wronge. Look in the command line for more details."),   
                 show_border = True),
                 title="Alert: Cannot perform the action!!",style="readonly",resizable=True,
                 buttons = [OKButton])
                 
    plot_output_cosmectics = View(Group(Item(name = 'max_depth',label="Maximum Depth of profile to plot"),
                 Item(name = 'label_size',label="Size of the data labels"),
                 Item(name = 'marker_size',label="Data points size"),
                 show_border = True),
                 title="Output plot cosmectics",resizable=True,
                 buttons = [OKButton])

    def __init__(self,litmod_path,fig,ax,save_option,add_option,edit_option,plot_option,run_option,edit_shape_option,delete_option,split_option,load_option,quit_option,close_model,merge_option,plot_results,X_S,X_E,Y_length,X_resolution):
        self.litmod_path=litmod_path
        ## buttons likning
        self.save_option = save_option;self.add_option = add_option;self.edit_option = edit_option;self.plot_option = plot_option
        self.run_option = run_option; self.edit_shape_option= edit_shape_option; self.delete_option=delete_option;self.split_option=split_option; self.quit_option=quit_option
        self.load_option= load_option; self.close_model=close_model;self.merge_option=merge_option;self.plot_results=plot_results 
        ## profile lenth and resolution in x direction
        self.X_S=float(str(X_S))
        self.X_E=float(str(X_E))
        self.Y=float(str(Y_length))
        self.reso=float(str(X_resolution))
        self.Z=[8.5,8.2,8.0,7.5,7.2,7.0,6.5,6.2,6.0,5.5,5.2,5.0,4.5,4.2,4.0,3.8,3.5,3.2,3.0,2.8,2.6,2.5,2.4,2.3,2.0,1.8,1.6,1.5,1.4,1.2,1.1,1.0,0.8,0.6,0.5,0.4,0.3,0.2,0.1,0.0,
          -0.5,  -1.5,  -2,  -2.5,  -3,  -3.5,  -4, \
        -4.5,  -5.5,  -6,  -6.5,  -7,  -8,  -9, -10, \
        -11, -12, -13, -14, -15, -17, -19, -21, \
        -23, -25, -27, -29, -31, -33, -35, -37, \
        -40, -43, -46, -49, -52, -55, -58, -61, \
        -64, -67, -70, -75, -80, -85, -90, -95, \
        -100,-105,-110,-115,-120,-125,-130,-135, \
        -140,-145,-150,-155,-160,-165,-170,-175, \
        -180,-185,-190,-195,-200,-205,-210,-215, \
        -220,-225,-230,-240,-250,-260,-270,-280, \
        -290,-300,-305,-310,-315,-320,-330,-340, \
        -350,-360,-370,-380,-390,-400,-410,-420]
        self.Z_for_shape=[0.0,  -0.5,  -1.5,  -2,  -2.5,  -3,  -3.5,  -4, \
        -4.5,  -5.5,  -6,  -6.5,  -7,  -8,  -9, -10, \
        -11, -12, -13, -14, -15, -17, -19, -21, \
        -23, -25, -27, -29, -31, -33, -35, -37, \
        -40, -43, -46, -49, -52, -55, -58, -61, \
        -64, -67, -70, -75, -80, -85, -90, -95, \
        -100,-105,-110,-115,-120,-125,-130,-135, \
        -140,-145,-150,-155,-160,-165,-170,-175, \
        -180,-185,-190,-195,-200,-205,-210,-215, \
        -220,-225,-230,-240,-250,-260,-270,-280, \
        -290,-300,-305,-310,-315,-320,-330,-340, \
        -350,-360,-370,-380,-390,-400]
        self.ax = ax;self.n = [];self.c = [];self.n_bod=0;self.mat = []; self.fig= fig
        #### density parameters
        self.dens = [];self.dens_depth = [];self.dens_temp = [];self.therm_expan = []
        ### thermal parameters
        self.conducx=[];self.conducy=[];self.conducz=[];self.depth_conduc_factor=[];self.temp_conduc_factor=[]
        self.heat_produ = [];self.heat_produ_depth = [];self.x_st= [];self.y_st= []
        ## variable related to coordinated of bodies and stuff
        self.x = [];self.y = [];self.x_end = [];self.y_end = [];self.l = [];self.b = []; self.x_st_index=[];self.y_st_index=[];self.l_text=[];self.b_text=[];self.x_envelope=[];self.y_envelope=[];
        ## variable for split body
        self.split_st_index=[];self.split_end_index=[];self.l_split=[];self.b_split=[]

        self.shape=Shape([(self.X_S,0),(self.X_E,0),(self.X_E,-400),(self.X_S,-400)])
        ## varibale related to anomaly coordinates
        self.anomaly_conut=0;self.anomaly_x=[];self.anomaly_y=[];self.anomaly_amount=[];self.body_Type=[];self.anomaly_compo=[];self.anomaly_Type=[];self.anomaly=0;self.anomaly_color=[]
        ## control variables
        self.body_status=1;
        
        self.ax.set_xlim((self.X_S, self.X_E))
        self.ax.set_ylim((-self.Y,8))
        self.ax.tick_params(labelsize=10)

        self.ax.grid(True)
        self.edit_traits(kind='modal', view='digitized_view')
        if self.digitized_response=="yes":
            f=np.loadtxt(self.digitized,comments=">")
            self.ax.plot(f[:,0],f[:,1],'*',lw=1,color='grey')
        else:
            pass
            
        self.canvas=self.fig.canvas
        self.canvas.mpl_connect('button_press_event', self._add_point)
        self.canvas.mpl_connect('button_press_event', self._close_polygon)
        self.canvas.mpl_connect('button_press_event', self._delete_point)

    def quit(self,event):
        self.edit_traits(kind='modal', view='exit_view')
        if self.exit_alert=="yes":
            exit()
        else:
            pass

    def set_location(self,event):
        if event.inaxes:
            self.x = round(event.xdata/float(self.reso))*float(self.reso)
            if self.body_type =="anomaly":
                self.y=round(event.ydata)
            else:
                self.y=event.ydata

    def _add_point(self,event):
        if event.button==2:
            self.vert.append((self.x,self.y))
            try:
                x = [self.vert[k][0] for k in range(len(self.vert))]
                y = [self.vert[k][1] for k in range(len(self.vert))]
                self.path.set_data(x,y)
            except Exception:
                pass
            self.path, = self.ax.plot([self.x],[self.y],'black',lw=1,marker='o')
            plt.draw()

    def _delete_point(self,event):
        if event.dblclick:
            if event.button==1:
                if len(self.vert)>0:
                    self.vert.pop()
                    x = [self.vert[k][0] for k in range(len(self.vert))]
                    y = [self.vert[k][1] for k in range(len(self.vert))]
                    self.path.set_data(x,y)
                    plt.draw()

    def read_model(self,event):
        if event.inaxes==self.load_option:
            if os.path.exists("litmod.inp")==True:
                if os.path.exists("bodies_GUI.out")==True:
                            self.X_S=[]
                            self.X_E=[]
                            self.Y=[]
                            self.reso=[]
                            self.n = [];self.c = [];self.n_bod=0;self.mat = []
                            self.dens = [];self.dens_depth = [];self.dens_temp = [];self.therm_expan = []
                            self.conducx=[];self.conducy=[];self.conducz=[];self.depth_conduc_factor=[];self.temp_conduc_factor=[]
                            self.heat_produ = [];self.heat_produ_depth = [];self.x_st= [];self.y_st= []
                            self.x = [];self.y = [];self.x_end = [];self.y_end = [];self.l = [];self.b = []; self.x_st_index=[];self.y_st_index=[];self.l_text=[];self.b_text=[]
                            self.x_envelope=[];self.y_envelope=[]
                            self.split_st_index=[];self.split_end_index=[];self.l_split=[];self.b_split=[]
                            self.anomaly=[];self.anomaly_conut=0;self.anomaly_x=[];self.anomaly_y=[];self.anomaly_amount=[];self.body_Type=[];self.anomaly_Type=[];self.anomaly_color=[]
                            fopen = open("litmod.inp")
                            content=fopen.readlines()
                            fopen.close()
                            
                            print(str(content[0].split()))
                            print(len((content[0].split())))
                            if len((content[0].split())) >=2:
                                self.topo=str(content[0].split()[0])
                            if len((content[1].split()))>=2:
                                self.bouger=str(content[1].split()[0])
                            if len((content[2].split()))>=2:
                                self.FA=str(content[2].split()[0])
                            if len((content[3].split()))>=2:
                                self.geoid=str(content[3].split()[0])
                            if len((content[11].split()))>=3:
                                self.anomaly_file=str(content[11].split()[0])
                                
                            self.n_mat=int(str(content[23].split()[0]))
                            self.n_bod=int(str(content[23].split()[1]))
                            self.anomaly=int(str(content[23].split()[2]))
                            self.print_option=int(str(content[23].split()[3]))
                            
                            j=1;temp=[]
                            for i in range(self.n_bod):
                                temp=content[24+(3*i)].split()
                                self.mat.append(int(str(temp[1])))
                                self.n.append(str(temp[2]))
                                temp=content[24+(3*i)+1].split()
                                self.heat_produ.append((str(temp[0])))
                                self.heat_produ_depth.append((str(temp[1])))
                                self.conducx.append((str(temp[2])))
                                self.conducy.append((str(temp[3])))
                                self.conducz.append((str(temp[4])))
                                self.temp_conduc_factor.append((str(temp[5])))
                                self.depth_conduc_factor.append((str(temp[6])))
                                temp=content[24+(3*i)+2].split()
                                self.dens.append((str(temp[0])))
                                self.dens_depth.append((str(temp[1])))
                                self.dens_temp.append((str(temp[2])))
                                self.therm_expan.append((str(temp[3])))
                                
                                if self.dens_temp[-1]==str(int(-1)) and float(str(self.dens[-1]))<0.0:
                                    self.anomaly_Type.append("Thermal")
                                    self.body_Type.append("anomaly")
                                    self.anomaly_amount.append(self.dens_depth[-1])
                                    self.anomaly_compo.append(self.mat[-1])
                                    self.n.pop(-1); self.mat.pop(-1); self.dens.pop(-1); self.dens_depth.pop(-1); self.dens_temp.pop(-1); self.therm_expan.pop(-1)
                                    self.conducx.pop(-1); self.conducy.pop(-1); self.conducz.pop(-1); self.depth_conduc_factor.pop(-1); self.temp_conduc_factor.pop(-1)
                                    self.heat_produ.pop(-1); self.heat_produ_depth.pop(-1)
                                elif self.dens_temp[-1]==str(int(1)) and float(str(self.dens[-1]))<0.0 and float(str(self.therm_expan[-1])) != 0.0 :
                                    self.anomaly_Type.append("Vp(%)")
                                    self.body_Type.append("anomaly")
                                    self.anomaly_amount.append(self.therm_expan[-1])
                                    self.anomaly_compo.append(self.mat[-1])
                                    self.n.pop(-1); self.mat.pop(-1); self.dens.pop(-1); self.dens_depth.pop(-1); self.dens_temp.pop(-1); self.therm_expan.pop(-1)
                                    self.conducx.pop(-1); self.conducy.pop(-1); self.conducz.pop(-1); self.depth_conduc_factor.pop(-1); self.temp_conduc_factor.pop(-1)
                                    self.heat_produ.pop(-1); self.heat_produ_depth.pop(-1)
                                elif self.dens_temp[-1]==str(int(2)):
                                    self.anomaly_Type.append("Vs(%)")
                                    self.body_Type.append("anomaly")
                                    self.anomaly_amount.append(self.therm_expan[-1])
                                    self.anomaly_compo.append(self.mat[-1])
                                    self.n.pop(-1); self.mat.pop(-1); self.dens.pop(-1); self.dens_depth.pop(-1); self.dens_temp.pop(-1); self.therm_expan.pop(-1)
                                    self.conducx.pop(-1); self.conducy.pop(-1); self.conducz.pop(-1); self.depth_conduc_factor.pop(-1); self.temp_conduc_factor.pop(-1)
                                    self.heat_produ.pop(-1); self.heat_produ_depth.pop(-1)
                                elif self.dens_temp[-1]==str(int(3)):
                                    self.anomaly_Type.append("Vp(>6.5)/Vs(<6.5)")
                                    self.body_Type.append("anomaly")
                                    self.anomaly_amount.append(self.therm_expan[-1])
                                    self.anomaly_compo.append(self.mat[-1])
                                    self.n.pop(-1); self.mat.pop(-1); self.dens.pop(-1); self.dens_depth.pop(-1); self.dens_temp.pop(-1); self.therm_expan.pop(-1)
                                    self.conducx.pop(-1); self.conducy.pop(-1); self.conducz.pop(-1); self.depth_conduc_factor.pop(-1); self.temp_conduc_factor.pop(-1)
                                    self.heat_produ.pop(-1); self.heat_produ_depth.pop(-1)
                                elif self.dens_temp[-1] == str(int(1)) and float(str(self.dens[-1])) < 0.0 and float(str(self.therm_expan[-1])) == 0.0 :
                                    self.anomaly_Type.append("Compositional")
                                    self.body_Type.append("anomaly")
                                    self.anomaly_amount.append("0.000E+00")
                                    self.anomaly_compo.append(self.mat[-1])
                                    self.n.pop(-1); self.mat.pop(-1); self.dens.pop(-1); self.dens_depth.pop(-1); self.dens_temp.pop(-1); self.therm_expan.pop(-1)
                                    self.conducx.pop(-1); self.conducy.pop(-1); self.conducz.pop(-1); self.depth_conduc_factor.pop(-1); self.temp_conduc_factor.pop(-1)
                                    self.heat_produ.pop(-1); self.heat_produ_depth.pop(-1)
                                else:
                                    self.body_Type.append("normal")

                            print("Types of bodies read:")
                            print(self.body_Type)
                            f=open("bodies_GUI.out",'r')
                            coords=f.readlines()
                            fopen.close()
                            
                            j=-1 
                            norm_body_count=-1
                            for i in range(1,len(coords)):
                                temp=coords[i].split()
                                if temp[0]!='>':
                                    if (self.body_Type[j] == "normal"):
                                        self.l[norm_body_count].append(float(str(temp[0])))
                                        self.b[norm_body_count].append(float(str(temp[1])))
                                        self.x_st[norm_body_count]=self.l[norm_body_count][0]
                                        self.x_st_index[norm_body_count]=self.l[norm_body_count][0]
                                        self.y_st[norm_body_count]=self.b[norm_body_count][0]
                                        self.y_st_index[norm_body_count]=self.b[norm_body_count][0]
                                        self.x_end[norm_body_count]=self.l[norm_body_count][-1]
                                        self.y_end[norm_body_count]=self.b[norm_body_count][-1]
                                    elif (self.body_Type[j] == "anomaly"): 
                                        self.anomaly_x[self.anomaly_conut-1].append(float(str(temp[0])))
                                        self.anomaly_y[self.anomaly_conut-1].append(float(str(temp[1])))
                                    else:
                                        pass
                                elif temp[0]=='>':
                                    j=j+1
                                    if (self.body_Type[j] == "normal"):
                                        norm_body_count=norm_body_count+1
                                        self.l.append([])
                                        self.b.append([])
                                        self.x_st.append([])
                                        self.x_end.append([])
                                        self.y_st.append([])
                                        self.y_end.append([])
                                        self.x_st_index.append([])
                                        self.y_st_index.append([])
                                        self.c.append(temp[2])
                                    else:
                                        self.anomaly_conut=self.anomaly_conut+1
                                        self.anomaly_x.append([])
                                        self.anomaly_y.append([])
                                        self.anomaly_color.append(temp[2])
                                else:
                                    pass
                            self.n_bod=len(self.l)
                            print('Number of bodies read =', self.n_bod) 
                            info=content[-15].split()
                            self.X_S=float(str(info[0]))/1000.0
                            self.reso=float(str(info[1]))/1000.0
                            X_nodes=int(str(info[2]))
                            self.X_E=(self.reso*X_nodes)-self.reso
                            self.X_E=self.X_E+self.X_S
                            print("Start and End of the profile")
                            print(self.X_S, self.X_E)
                            self.Y=400.0
                            self.shape=Shape([(self.X_S,0),(self.X_E,0),(self.X_E,-400),(self.X_S,-400)])
                            self.ax.set_xlim((self.X_S, self.X_E))
                            self.ax.set_ylim((-self.Y,10))
                            self.ax.tick_params(labelsize=10)
                            self.ax.grid(True)
                            
                            f=open("bodies_GUI_envelops.out",'r')
                            coords=f.readlines()
                            j=-1 
                            norm_body_count=-1
                            for i in range(1,len(coords)):
                                temp=coords[i].split()
                                if temp[0]!='>':
                                    self.x_envelope[norm_body_count].append(float(str(temp[0])))
                                    self.y_envelope[norm_body_count].append(float(str(temp[1])))
                                elif temp[0]=='>':
                                    j=j+1
                                    norm_body_count=norm_body_count+1
                                    self.x_envelope.append([])
                                    self.y_envelope.append([])
                            
                            self.x_envelope.append(self.l[len(self.l)-1])
                            self.y_envelope.append(self.b[len(self.l)-1])
                            print('No of envelops read =', len(self.x_envelope))
                            self.vert = []
                            self.path = []
                            self.shape=Shape( [ (self.x_envelope[self.n_bod-1][k],self.y_envelope[self.n_bod-1][k]) 
                                for k in range(len(self.x_envelope[self.n_bod-1]))] )
                            xa=[];ya=[]
                            xa,ya=self.shape.exterior.xy
                            for i in range(len(xa)):
                                self.vert.append((xa[i],ya[i]))
                            x = [self.vert[k][0] for k in range(len(self.vert))]
                            y = [self.vert[k][1] for k in range(len(self.vert))]
                            self.path, = self.ax.plot([self.x_e],[self.y_e],'o-',lw=1,color=self.c[self.n_bod-1])
                            self.path, = self.ax.plot(x,y,'o-',lw=1,color=self.c[self.n_bod-1])
                            self.path, = self.ax.plot(x,y,color=self.c[self.n_bod-1],lw=1,marker='o')
                            if self.anomaly_conut != 0:
                                self.edit_traits(kind='modal', view='anomaly_loaded_view')
                                print("Anomalous bodies information")
                                print(self.anomaly_Type)
                                print(self.anomaly_amount)
                                print(self.anomaly_compo)
                                self.anomaly=1
                            else:
                                pass
                else:
                    print("No bodies_GUI.out file")
                    self.edit_traits(kind='modal', view='plot_output_view') 
            else:
                print("No litmod.inp input file")
                self.edit_traits(kind='modal', view='plot_output_view') 
                
    def Plot_results(self,event):
        if event.inaxes == self.plot_results:
            self.edit_traits(kind='modal', view='plot_output_cosmectics')
            plot_lib.plot_lib(len(range(int(float(self.X_S)),int(float(self.X_E)),int(float(self.reso))))+1,
                self.l,self.b,self.topo,self.bouger,self.geoid,self.FA,self.SHF,self.anomaly,self.c,self.l,self.b,self.mat,self.max_depth,self.label_size,
                self.marker_size,self.anomaly_x,self.anomaly_y,self.anomaly_amount,self.anomaly_Type,self.anomaly_compo,self.anomaly_color)

    def run_model(self,event):
            if self.body_status==1 and self.mat[self.n_bod-1]==99:
                if os.path.exists("litmod.inp")==True:
                    print("running the model")
                    
                    # --- LINUX ---
                    if platform.system()=="Linux":
                        os.system(self.litmod_path+"/LITMOD_V4_LINUX")
                        try:
                            self.edit_traits(kind='modal', view='plot_output_cosmectics')
                            plot_lib.plot_lib(len(range(int(float(self.X_S)),int(float(self.X_E)),int(float(self.reso))))+1,
                                self.l,self.b,self.topo,self.bouger,self.geoid,self.FA,self.SHF,self.anomaly,self.c,self.l,self.b,self.mat,self.max_depth,self.label_size,
                                self.marker_size,self.anomaly_x,self.anomaly_y,self.anomaly_amount,self.anomaly_Type,self.anomaly_compo,self.anomaly_color)
                        except Exception as e:
                            import traceback
                            traceback.print_exc()
                            print("Something is wrong. Look above if LitMod ran without errors. If yes than something is wrong with plotting part.")
                            self.edit_traits(kind='modal', view='plot_output_view')
                            
                    # --- WINDOWS ---
                    elif platform.system()=="Windows":
                        os.system(f'"{self.litmod_path}\\LITMOD_V4_Windows.exe"')
                        try:
                            self.edit_traits(kind='modal', view='plot_output_cosmectics')
                            plot_lib.plot_lib(len(range(int(float(self.X_S)),int(float(self.X_E)),int(float(self.reso))))+1,
                                self.l,self.b,self.topo,self.bouger,self.geoid,self.FA,self.SHF,self.anomaly,self.c,self.l,self.b,self.mat,self.max_depth,self.label_size,
                                self.marker_size,self.anomaly_x,self.anomaly_y,self.anomaly_amount,self.anomaly_Type,self.anomaly_compo,self.anomaly_color)
                        except Exception as e:
                            import traceback
                            traceback.print_exc()
                            print("Something is wrong. Look above if LitMod ran without errors. If yes than something is wrong with plotting part.")
                            self.edit_traits(kind='modal', view='plot_output_view')
                            
                    # --- MACOS (Darwin) ---
                    elif platform.system()=="Darwin":
                        # Note: Ensure the macOS executable is named 'LITMOD_V4_MAC' and has execute permissions (chmod +x).
                        os.system(self.litmod_path+"/LITMOD_V4_MAC")
                        try:
                            self.edit_traits(kind='modal', view='plot_output_cosmectics')
                            plot_lib.plot_lib(len(range(int(float(self.X_S)),int(float(self.X_E)),int(float(self.reso))))+1,
                                self.l,self.b,self.topo,self.bouger,self.geoid,self.FA,self.SHF,self.anomaly,self.c,self.l,self.b,self.mat,self.max_depth,self.label_size,
                                self.marker_size,self.anomaly_x,self.anomaly_y,self.anomaly_amount,self.anomaly_Type,self.anomaly_compo,self.anomaly_color)
                        except Exception as e:
                            import traceback
                            traceback.print_exc()
                            print("Something is wrong. Look above if LitMod ran without errors. If yes than something is wrong with plotting part.")
                            self.edit_traits(kind='modal', view='plot_output_view')
                            
                else:
                    print("LitMod input file (litmod.inp) is not there")
                    self.edit_traits(kind='modal', view='plot_output_view')          
            else:
                print("Cannot run LitMod. Either model is not closed or you have a body added which is not closed. Close the body and close the model.")
                self.edit_traits(kind='modal', view='plot_output_view')

    def plot_model(self,event):
        if self.body_status == 1:
            self.ax.cla()
            self.ax.set_xlim((self.X_S, self.X_E))
            self.ax.set_ylim((-self.Y,10))
            self.ax.tick_params(labelsize=10)
            self.ax.set_xlabel('Distance (Km)')
            self.ax.set_ylabel('Depth (Km)')
            self.ax.grid(True)            
            if os.path.exists("normal_body_shape.dat")==True:
                self.x_envelope=[]
                self.y_envelope=[]
                f=open("normal_body_shape.dat")
                coords=f.readlines()
                f.close()
                j=-1
                for i in range(1,len(coords)):
                    temp=coords[i].split()
                    if(temp[0]!='>'):
                         self.x_envelope[j].append(float(str(temp[0])))
                         self.y_envelope[j].append(float(str(temp[1])))
                    if(temp[0]=='>'):
                        j=j+1
                        self.x_envelope.append([])
                        self.y_envelope.append([])
                print("Total number of envelops recovered")
                print(len(self.x_envelope))
                os.remove("normal_body_shape.dat")
                self.l=[]
                self.b=[]
                ### getting envelopes
                for i in range(len(self.x_envelope)-1):
                    poly_1=[]
                    poly_2=[]
                    a=Shape() 
                    self.l.append([])
                    self.b.append([])
                    poly_1=Shape( [(self.x_envelope[i][j],self.y_envelope[i][j]) for j in range(len(self.x_envelope[i])) ])
                    poly_2=Shape( [(self.x_envelope[i+1][j],self.y_envelope[i+1][j]) for j in range(len(self.x_envelope[i+1])) ])
                    print(str(i) + ' body extracted')
                    a=poly_1.difference(poly_2)
                    self.l[i],self.b[i]=a.exterior.coords.xy
                self.l.append([])
                self.b.append([])
                self.y_envelope.append(self.y_envelope[-1])
                self.x_envelope.append(self.x_envelope[-1])
                self.l[-1],self.b[-1]=poly_2.exterior.coords.xy
                print("Total number of bodies recovered")
                print(len(self.l))
            
            if os.path.exists("anomaly_body_shape.dat")==True:
                self.anomaly_x=[]
                self.anomaly_y=[]
                f=open("anomaly_body_shape.dat")
                coords=f.readlines()
                f.close()
                j=-1
                for i in range(1,len(coords)):
                    temp=coords[i].split()
                    if(temp[0]!='>'):
                         self.anomaly_x[j].append(float(str(temp[0])))
                         self.anomaly_y[j].append(float(str(temp[1])))
                    if(temp[0]=='>'):
                        j=j+1
                        self.anomaly_x.append([])
                        self.anomaly_y.append([])
                print("Total number of envelops recovered")
                print(len(self.x_envelope))
                os.remove("anomaly_body_shape.dat")

            if self.n_bod>0:       
                for i in range(len(self.l)):
                    self.ax.plot(self.l[i],self.b[i],color='black',lw=1,marker='o',markersize=4)
                    poly = Polygon(list(zip(self.l[i],self.b[i])),facecolor=self.c[i],label=self.mat[i], lw=1)
                    self.ax.add_patch(poly)
                    self.l_text=min(self.l[i][:]) + (max(self.l[i][:])-min(self.l[i][:]))/2.0
                    self.b_text=min(self.b[i][:]) + (max(self.b[i][:])-min(self.b[i][:]))/2.0
                    self.ax.text(self.l_text, self.b_text, str(int(i+1))+" "+str(self.n[i])+" Material "
                        +str(self.mat[i]), fontsize=8,fontstyle='italic',color='black',bbox=dict(facecolor='white', alpha=0.9))
                self.ax.plot(self.x_envelope[-1],self.y_envelope[-1])
                
                for i in range(len(self.anomaly_x)):
                   self.ax.plot(self.anomaly_x[i],self.anomaly_y[i],color="black",lw=1,marker='o',markersize=4)
                   poly = Polygon(list(zip(self.anomaly_x[i],self.anomaly_y[i])),facecolor=self.anomaly_color[i],label="anomaly", lw=1)
                   self.l_text=min(self.anomaly_x[i])+ (max(self.anomaly_x[i])-min(self.anomaly_x[i]))/2.0
                   self.b_text=min(self.anomaly_y[i])+ (max(self.anomaly_y[i])-min(self.anomaly_y[i]))/2.0
                   self.ax.add_patch(poly)
                   self.ax.text(self.l_text, self.b_text, str(int(i+1))+" "+str(self.anomaly_Type[i])+" Material "
                    +str(self.anomaly_compo[i])+"Value "+str(self.anomaly_amount[i]), fontsize=8,fontstyle='italic',color='black',bbox=dict(facecolor='white', alpha=0.9))
                try:
                    f=np.loadtxt(self.digitized,comments=">")
                    self.ax.plot(f[:,0],f[:,1],'*',lw=1,color='grey')
                except:
                    pass
                plt.draw()
            else:
                try:
                    f=np.loadtxt(self.digitized,comments=">")
                    self.ax.plot(f[:,0],f[:,1],'*',lw=1,color='grey')
                except:
                    pass
                print("There is nothing to plot")
                self.edit_traits(kind='modal', view='plot_output_view') 
                pass
        else:
            print("Cannot plot. You have a body added which is not closed yet.")
            self.edit_traits(kind='modal', view='plot_output_view')

    def save_body(self,event):
        if self.body_status == 1: 
            if self.anomaly_conut==0:
                self.anomaly=0
            else:
                self.anomaly=1
            self.edit_traits(kind='modal', view='save_view')
            
            # ==============================================================
            # 1. Filter duplicate coordinate points (prevents Fortran geometry intersection errors)
            # ==============================================================
            for j in range(len(self.l)):
                clean_x, clean_y = [], []
                for i in range(len(self.l[j])):
                    px, py = self.l[j][i], self.b[j][i]
                    if not clean_x or (px != clean_x[-1] or py != clean_y[-1]):
                        clean_x.append(px); clean_y.append(py)
                self.l[j] = clean_x; self.b[j] = clean_y
                
            for j in range(len(self.anomaly_x)):
                clean_x, clean_y = [], []
                for i in range(len(self.anomaly_x[j])):
                    px, py = self.anomaly_x[j][i], self.anomaly_y[j][i]
                    if not clean_x or (px != clean_x[-1] or py != clean_y[-1]):
                        clean_x.append(px); clean_y.append(py)
                self.anomaly_x[j] = clean_x; self.anomaly_y[j] = clean_y
            
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            backup_targets = {
                "litmod.inp": "litmod_inp",
                "bodies_GUI.out": "bodies_GUI",
                "bodies_GUI_envelops.out": "bodies_GUI_envelops",
            }
            for src, stem in backup_targets.items():
                if os.path.exists(src):
                    try:
                        copyfile(src, f"{stem}_{timestamp}.bak")
                    except Exception as exc:
                        print(f"Could not create backup for {src}: {exc}")
            
            f=open("litmod.inp","w")
            f.writelines('{:74s} {:s}\n'.format(self.topo,"ELEVin"))
            f.writelines('{:74s} {:s}\n'.format(self.bouger,"Bouguer"))
            f.writelines('{:74s} {:s}\n'.format(self.FA,"Free-air"))
            f.writelines('{:74s} {:s}\n'.format(self.geoid,"GEOIDin"))
            f.writelines('{:74s} {:s}\n'.format("topo_out.dat","Topograp"))
            f.writelines('{:74s} {:s}\n'.format("bouguer_out.dat","Bouguer"))
            f.writelines('{:74s} {:s}\n'.format("FA_out.dat","Free air"))
            f.writelines('{:74s} {:s}\n'.format("        ","Mixed"))
            f.writelines('{:74s} {:s}\n'.format("geoid_out.dat","GRAVout"))
            f.writelines('{:74s} {:s}\n'.format("tempout.dat","Temperat"))
            f.writelines('{:74s} {:s}\n'.format("SHF_out.dat","SHF out"))
            if len(self.anomaly_file)>0:
                f.writelines('{:74s} {:s}\n'.format(self.anomaly_file,"CSTR out"))
            else:   
                f.writelines('{:74s} {:s}\n'.format("            ","CSTR out"))
            f.writelines('{:74s} {:s}\n'.format("            ","ESTR out"))
            f.writelines('{:74s} {:s}\n'.format("            ","COL out"))
            f.writelines('{:74s} {:s}\n'.format("            ","BOUN out"))
            f.writelines('{:74s} {:s}\n'.format("bodies.out","BODY out"))
            f.writelines('{:74s} {:s}\n'.format("            ","ELEM out"))
            f.writelines('{:74s} {:s}\n'.format("            ","ISTR out"))
            f.writelines('{:74s} {:s}\n'.format("P_T_out.dat","P/T out"))
            f.writelines('{:74s} {:s}\n'.format("dens.dat","Density"))
            f.writelines('{:74s} {:s}\n'.format("            ","RECA out"))
            f.writelines(["%s                                                                         %s\n" % (self.thermo_data,"Thermo D")])
            f.writelines(["%s" % ("comment\n")])
            
            if self.anomaly_conut==0:
                if len(self.anomaly_file)>0: self.anomaly=1
                else: self.anomaly=0
                f.writelines(["   %s    %s    %s    %s\n" % (len(self.l),len(self.l),self.anomaly,self.print_option)])                
                for i in range(len(self.l)):
                    f.writelines(["      %s      %s      %s\n" % (i+1,self.mat[i],self.n[i],)])
                    f.writelines(["  %s  %s  %s  %s  %s  %s  %s\n" % (self.heat_produ[i],self.heat_produ_depth[i], \
                            self.conducx[i],self.conducy[i],self.conducz[i],self.temp_conduc_factor[i],self.depth_conduc_factor[i])])
                    f.writelines(["  %s  %s  %s  %s\n" % (self.dens[i],self.dens_depth[i],self.dens_temp[i],self.therm_expan[i])])
                
                for j in range(len(self.l)):
                    f.writelines(["%s                  %s\n" % (j+1,self.n[j])])
                    for i in range(len(self.l[j])):
                        # [Maintain Legacy Format]: Must use 7d and %4! The Fortran 8F8.0 format strictly requires this!
                        if i==(len(self.l[j])-1) or (i+1)%4 ==0:
                            f.writelines('{:7d}{:s}{:7d}{:s}'.format(int(float(self.l[j][i]*1000)),".",int(float(self.b[j][i]*1000)),"."))
                            f.writelines("\n")
                        else:
                            f.writelines('{:7d}{:s}{:7d}{:s}'.format(int(float(self.l[j][i]*1000)),".",int(float(self.b[j][i]*1000)),"."))
                f.writelines(["%s %s %s %s %s\n" % ("0","1","1","0","0")])
                f.writelines(["%s%s%s%s%s\n" % ("0.,","0.,","0.,","0.,","0.")])
                f.writelines(["%s%s%s%s%s\n" % ("0.,",self.T_surf+".,",self.T_lab+".,","0.,","0.")])
                f.writelines(["%s " " %s " " %s " " %s\n" % (int(self.X_S*1000),int(self.reso*1000),len(range(int(float(self.X_S)),int(float(self.X_E)),int(float(self.reso))))+1,"96")])
                
                # [Maintain Legacy Format]: Restore compact layout for Z-levels
                f.write("      0.   -500.  -1500.  -2000.  -2500.  -3000.  -3500.  -4000.\n\
  -4500.  -5500.  -6000.  -6500.  -7000.  -8000.  -9000. -10000.\n\
 -11000. -12000. -13000. -14000. -15000. -17000. -19000. -21000.\n\
 -23000. -25000. -27000. -29000. -31000. -33000. -35000. -37000.\n\
 -40000. -43000. -46000. -49000. -52000. -55000. -58000. -61000.\n\
 -64000. -67000. -70000. -75000. -80000. -85000. -90000. -95000.\n\
-100000.-105000.-110000.-115000.-120000.-125000.-130000.-135000.\n\
-140000.-145000.-150000.-155000.-160000.-165000.-170000.-175000.\n\
-180000.-185000.-190000.-195000.-200000.-205000.-210000.-215000.\n\
-220000.-225000.-230000.-240000.-250000.-260000.-270000.-280000.\n\
-290000.-300000.-305000.-310000.-315000.-320000.-330000.-340000.\n\
-350000.-360000.-370000.-380000.-390000.-400000.-410000.-420000.\n\
0,4,5,3,4\n\
0,0,0,0,0")
                f.close()
                
                f=open("bodies_GUI.out","w")
                f.writelines(">  \t  %d\n " %(len(self.l)))
                for j in range(len(self.l)):
                    for i in range(len(self.l[j])):
                        if i==(len(self.l[j])-len(self.l[j])):
                            f.writelines(">  \t  %d  %s %s %s\n " % (len(self.l[j]),self.c[j],self.dens[j],self.n[j]))
                        f.writelines(["%f \t %f \n" % (self.l[j][i] , self.b[j][i])])
                f.close()
                
                # ==============================================================
                # [Critical Sync]: Force synchronization of node data with the Fortran engine!
                # Fixes crashes caused by Fortran misreading line counts!
                # ==============================================================
                import shutil
                shutil.copyfile("bodies_GUI.out", "bodies.out")
                
                f=open("bodies_GUI_envelops.out","w")
                f.writelines(">  \t  %d\n " %(len(self.l)))
                for j in range(len(self.x_envelope)-1):
                    for i in range(len(self.x_envelope[j])):
                        if i==(len(self.x_envelope[j])-len(self.x_envelope[j])):
                            f.writelines(">  \t  %d  %s\n " % (len(self.x_envelope[j]),self.c[j]))
                        f.writelines(["%f \t %f \n" % (self.x_envelope[j][i] , self.y_envelope[j][i])])
                f.close()

            else:
                self.anomaly=1
                for anom_i in range(self.anomaly_conut):
                    anom_mat = self.anomaly_compo[anom_i]
                    try:
                        if int(anom_mat) <= 0:
                            anom_mat = 99
                    except (TypeError, ValueError):
                        anom_mat = 99
                    self.c.insert(len(self.l)-1,self.anomaly_color[anom_i])
                    self.n.insert(len(self.l)-1,self.anomaly_Type[anom_i])
                    self.conducx.insert(len(self.l)-1,"3.300E+00")
                    self.conducy.insert(len(self.l)-1,"3.300E+00")
                    self.conducz.insert(len(self.l)-1,"0.000E+00")
                    self.temp_conduc_factor.insert(len(self.l)-1,"1.500E-03")
                    self.depth_conduc_factor.insert(len(self.l)-1,"1.300E+02")
                    self.heat_produ.insert(len(self.l)-1,"0.000E+00")
                    if self.anomaly_Type[anom_i]=="Thermal":
                        self.mat.insert(len(self.l)-1,anom_mat)
                        self.heat_produ_depth.insert(len(self.l)-1,"0.000E+00")
                        self.dens.insert(len(self.l)-1,"-3.100E+03")
                        self.dens_depth.insert(len(self.l)-1,self.anomaly_amount[anom_i])
                        self.dens_temp.insert(len(self.l)-1,"-1")
                        self.therm_expan.insert(len(self.l)-1,"0.000E+00")
                        self.l.insert(len(self.l)-1,self.anomaly_x[anom_i][:])
                        self.b.insert(len(self.b)-1,self.anomaly_y[anom_i][:])
                    elif self.anomaly_Type[anom_i]=="Vp(%)":
                        self.mat.insert(len(self.l)-1,anom_mat)
                        self.heat_produ_depth.insert(len(self.l)-1,"0.000E+00")
                        self.dens.insert(len(self.l)-1,"-3.100E+03")
                        self.dens_depth.insert(len(self.l)-1,"0.000E+00")
                        self.dens_temp.insert(len(self.l)-1,"1")
                        self.therm_expan.insert(len(self.l)-1,self.anomaly_amount[anom_i])
                        self.l.insert(len(self.l)-1,self.anomaly_x[anom_i][:])
                        self.b.insert(len(self.b)-1,self.anomaly_y[anom_i][:])
                    elif self.anomaly_Type[anom_i]=="Vs(%)":
                        self.mat.insert(len(self.l)-1,anom_mat)
                        self.heat_produ_depth.insert(len(self.l)-1,"0.000E+00")
                        self.dens.insert(len(self.l)-1,"-3.100E+03")
                        self.dens_depth.insert(len(self.l)-1,"0.000E+00")
                        self.dens_temp.insert(len(self.l)-1,"2")
                        self.therm_expan.insert(len(self.l)-1,self.anomaly_amount[anom_i])
                        self.l.insert(len(self.l)-1,self.anomaly_x[anom_i][:])
                        self.b.insert(len(self.b)-1,self.anomaly_y[anom_i][:])
                    elif self.anomaly_Type[anom_i]=="Vp(>6.5)/Vs(<6.5)":
                        self.mat.insert(len(self.l)-1,anom_mat)
                        self.heat_produ_depth.insert(len(self.l)-1,"0.000E+00")
                        self.dens.insert(len(self.l)-1,"-3.100E+03")
                        self.dens_depth.insert(len(self.l)-1,"0.000E+00")
                        self.dens_temp.insert(len(self.l)-1,"3")
                        self.therm_expan.insert(len(self.l)-1,self.anomaly_amount[anom_i])
                        self.l.insert(len(self.l)-1,self.anomaly_x[anom_i][:])
                        self.b.insert(len(self.b)-1,self.anomaly_y[anom_i][:])
                    elif self.anomaly_Type[anom_i]=="Compositional":
                        self.mat.insert(len(self.l)-1,anom_mat)
                        self.heat_produ_depth.insert(len(self.l)-1,"0.000E+00")
                        self.dens.insert(len(self.l)-1,"-3.100E+03")
                        self.dens_depth.insert(len(self.l)-1,"0.000E+00")
                        self.dens_temp.insert(len(self.l)-1,"1")
                        self.therm_expan.insert(len(self.l)-1,"0.000E+00")
                        self.l.insert(len(self.l)-1,self.anomaly_x[anom_i][:])
                        self.b.insert(len(self.b)-1,self.anomaly_y[anom_i][:])
                f.writelines(["   %s    %s    %s    %s\n" % (len(self.l),len(self.l),self.anomaly,self.print_option)])                
                for i in range(len(self.l)):
                    f.writelines(["      %s      %s      %s\n" % (i+1,self.mat[i],self.n[i],)])
                    f.writelines(["  %s  %s  %s  %s  %s  %s  %s\n" % (self.heat_produ[i],self.heat_produ_depth[i], \
                            self.conducx[i],self.conducy[i],self.conducz[i],self.temp_conduc_factor[i],self.depth_conduc_factor[i])])
                    f.writelines(["  %s  %s  %s  %s\n" % (self.dens[i],self.dens_depth[i],self.dens_temp[i],self.therm_expan[i])])
                
                for j in range(len(self.l)):
                    f.writelines(["%s                  %s\n" % (j+1,self.n[j])])
                    for i in range(len(self.l[j])):
                        if i==(len(self.l[j])-1) or (i+1)%4 ==0:
                            f.writelines('{:7d}{:s}{:7d}{:s}'.format(int(float(self.l[j][i]*1000)),".",int(float(self.b[j][i]*1000)),"."))
                            f.writelines("\n")
                        else:
                            f.writelines('{:7d}{:s}{:7d}{:s}'.format(int(float(self.l[j][i]*1000)),".",int(float(self.b[j][i]*1000)),"."))
                f.writelines(["%s %s %s %s %s\n" % ("0","1","1","0","0")])
                f.writelines(["%s%s%s%s%s\n" % ("0.,","0.,","0.,","0.,","0.")])
                f.writelines(["%s%s%s%s%s\n" % ("0.,",self.T_surf+".,",self.T_lab+".,","0.,","0.")])
                f.writelines(["%s " " %s " " %s " " %s\n" % (int(self.X_S*1000),int(self.reso*1000),len(range(int(float(self.X_S)),int(float(self.X_E)),int(float(self.reso))))+1,"96")])
                
                f.write("      0.   -500.  -1500.  -2000.  -2500.  -3000.  -3500.  -4000.\n\
  -4500.  -5500.  -6000.  -6500.  -7000.  -8000.  -9000. -10000.\n\
 -11000. -12000. -13000. -14000. -15000. -17000. -19000. -21000.\n\
 -23000. -25000. -27000. -29000. -31000. -33000. -35000. -37000.\n\
 -40000. -43000. -46000. -49000. -52000. -55000. -58000. -61000.\n\
 -64000. -67000. -70000. -75000. -80000. -85000. -90000. -95000.\n\
-100000.-105000.-110000.-115000.-120000.-125000.-130000.-135000.\n\
-140000.-145000.-150000.-155000.-160000.-165000.-170000.-175000.\n\
-180000.-185000.-190000.-195000.-200000.-205000.-210000.-215000.\n\
-220000.-225000.-230000.-240000.-250000.-260000.-270000.-280000.\n\
-290000.-300000.-305000.-310000.-315000.-320000.-330000.-340000.\n\
-350000.-360000.-370000.-380000.-390000.-400000.-410000.-420000.\n\
0,4,5,3,4\n\
0,0,0,0,0")
                f.close()
                
                f=open("bodies_GUI.out","w")
                f.writelines(">  \t  %d\n " %(len(self.l)))
                for j in range(len(self.l)):
                    for i in range(len(self.l[j])):
                        if i==(len(self.l[j])-len(self.l[j])):
                            f.writelines(">  \t  %d  %s %s %s\n " % (len(self.l[j]),self.c[j],self.dens[j],self.n[j]))
                        f.writelines(["%f \t %f \n" % (self.l[j][i] , self.b[j][i])])
                f.close()
                
                # ==============================================================
                # Synchronize copies of the selected geometry.
                # ==============================================================
                import shutil
                shutil.copyfile("bodies_GUI.out", "bodies.out")

                f=open("bodies_GUI_envelops.out","w")
                f.writelines(">  \t  %d\n " %(len(self.l)))
                for j in range(len(self.x_envelope)-1):
                    for i in range(len(self.x_envelope[j])):
                        if i==(len(self.x_envelope[j])-len(self.x_envelope[j])):
                            f.writelines(">  \t  %d  %s\n " % (len(self.x_envelope[j]),self.c[j]))
                        f.writelines(["%f \t %f \n" % (self.x_envelope[j][i] , self.y_envelope[j][i])])
                f.close()
                for anom_i in range(self.anomaly_conut):
                    self.mat.pop(len(self.l)-2)
                    self.c.pop(len(self.l)-2)
                    self.n.pop(len(self.l)-2)
                    self.conducx.pop(len(self.l)-2)
                    self.conducy.pop(len(self.l)-2)
                    self.conducz.pop(len(self.l)-2)
                    self.temp_conduc_factor.pop(len(self.l)-2)
                    self.depth_conduc_factor.pop(len(self.l)-2)
                    self.heat_produ.pop(len(self.l)-2)
                    self.heat_produ_depth.pop(len(self.l)-2)
                    self.dens.pop(len(self.l)-2)
                    self.dens_depth.pop(len(self.l)-2)
                    self.dens_temp.pop(len(self.l)-2)
                    self.therm_expan.pop(len(self.l)-2)
                    self.l.pop(len(self.l)-2)
                    self.b.pop(len(self.b)-2)

        else:
            print("Cannot save model yet. Either you have body added which is not closed yet or you have not closed the model by clicking Close Model button at top mid of the window.")
            self.edit_traits(kind='modal', view='plot_output_view')
            
    #################################################################################
#TO split a body into two
#################################################################################
    def split_body(self,event):
        if self.n_bod>0 and self.body_status == 1:
            self.edit_traits(kind='modal', view='split_view')
            if self.split_body_no > 0 and self.ans=='yes':
                self.body_type='split'
                
                # Automatically inherit original name and append _split for identification
                self.name = self.n[self.split_body_no-1] + "_split" 
                
                self.edit_traits(kind='modal', view='add_view')
                if self.material >0:
                    self.Temp_density_factor='1'
                    self.ref_density = "3.245E+03"
                    self.depth_density_factor= "0.000E+00"
                    self.Thermal_expansion= "0.000E+00"
                    self.Kx ="2.100E+00"
                    self.Ky ="2.100E+00"
                    self.Kz ="0.000E+00"
                    self.Depth_conduc_factor="1.300E+02"
                    self.Temp_conduc_factor="1.500E-03" 
                    pass
                else:
                    self.Depth_conduc_factor="0.000E+00"
                    self.Temp_conduc_factor="0.000E+00" 
                    self.Temp_density_factor='0'
                    self.body_type='split'
                    self.edit_traits(kind='modal', view='body_split')
                
                self.body_status=0;
                self.l_split=self.l[self.split_body_no-1]
                self.b_split=self.b[self.split_body_no-1]
                self.path=[]
                self.vert = []
            elif self.ans=='no':
                pass
            else:
                self.body_type='normal'
                print("Entered body number is not valid.")
                self.edit_traits(kind='modal', view='plot_output_view') 
        else: 
            print("Cannot split. You have a body added which is not closed yet or you do not have any body in the model.")
            self.edit_traits(kind='modal', view='plot_output_view')

    def merge_bodies(self,event):
        if self.n_bod>0 and self.body_status == 1:
            self.edit_traits(kind='modal', view='merge_view')
            if self.merge_body_no_1 <= self.n_bod and self.merge_body_no_2 <= self.n_bod and self.ans=='yes':
                
                idx1 = self.merge_body_no_1 - 1
                idx2 = self.merge_body_no_2 - 1
                
                if idx1 == idx2:
                    print("Cannot merge the same body.")
                    return
                
                # Determine the indices to retain and remove.
                keep_idx = min(idx1, idx2)
                del_idx = max(idx1, idx2)
                prop_idx = idx2  # Retain the physical properties of the second selected body.
                
                # ==============================================================
                # Merge the selected body polygons with Shapely.
                # Use the polygon union to obtain the merged boundary.
                # ==============================================================
                p1 = Shape(list(zip(self.l[idx1], self.b[idx1]))).buffer(0)
                p2 = Shape(list(zip(self.l[idx2], self.b[idx2]))).buffer(0)
                
                merged_poly = p1.union(p2).buffer(0)
                
                if merged_poly.geom_type == 'MultiPolygon':
                    merged_poly = max(merged_poly.geoms, key=lambda p: p.area)
                    
                m_x, m_y = merged_poly.exterior.coords.xy
                
                self.l[keep_idx] = list(m_x)
                self.b[keep_idx] = list(m_y)
                
                # ==============================================================
                # [Fix 2: Perfect fusion of envelopes (No more "spider web" artifacts!)]
                # Core Logic: Maintain strictly increasing X coordinates, set Y to the "shallowest depth (maximum value)"
                # ==============================================================
                env1_x = self.x_envelope[idx1]
                env1_y = self.y_envelope[idx1]
                env2_x = self.x_envelope[idx2]
                env2_y = self.y_envelope[idx2]
                
                if len(env1_x) == len(env2_x):
                    # Merge corresponding envelope samples using the shallower boundary.
                    new_env_y = [max(y1, y2) for y1, y2 in zip(env1_y, env2_y)]
                    self.x_envelope[keep_idx] = list(env1_x)
                    self.y_envelope[keep_idx] = new_env_y
                else:
                    # For different sample counts, choose an envelope using the average depth.
                    avg1 = sum(env1_y) / len(env1_y)
                    avg2 = sum(env2_y) / len(env2_y)
                    if avg2 > avg1:
                        self.x_envelope[keep_idx] = list(env2_x)
                        self.y_envelope[keep_idx] = list(env2_y)
                
                # ==============================================================
                # [Fix 3: Property inheritance and memory cleanup]
                # ==============================================================
                self.n[keep_idx] = self.n[prop_idx]
                self.c[keep_idx] = self.c[prop_idx]
                self.conducx[keep_idx] = self.conducx[prop_idx]
                self.conducy[keep_idx] = self.conducy[prop_idx]
                self.conducz[keep_idx] = self.conducz[prop_idx]
                self.depth_conduc_factor[keep_idx] = self.depth_conduc_factor[prop_idx]
                self.temp_conduc_factor[keep_idx] = self.temp_conduc_factor[prop_idx]
                self.mat[keep_idx] = self.mat[prop_idx]
                self.heat_produ[keep_idx] = self.heat_produ[prop_idx]
                self.heat_produ_depth[keep_idx] = self.heat_produ_depth[prop_idx]
                self.dens[keep_idx] = self.dens[prop_idx]
                self.dens_depth[keep_idx] = self.dens_depth[prop_idx]
                self.dens_temp[keep_idx] = self.dens_temp[prop_idx]
                self.therm_expan[keep_idx] = self.therm_expan[prop_idx]
                self.body_Type[keep_idx] = self.body_Type[prop_idx]
                
                self.x_st_index[keep_idx] = list(m_x)[0]
                self.y_st_index[keep_idx] = list(m_y)[0]
                self.x_end[keep_idx] = list(m_x)[-1]
                self.y_end[keep_idx] = list(m_y)[-1]
                self.x_st[keep_idx] = list(m_x)[0]
                self.y_st[keep_idx] = list(m_y)[0]
                
                # Remove the merged-away body from all parallel property arrays.
                del self.l[del_idx]; del self.b[del_idx]
                del self.n[del_idx]; del self.c[del_idx]
                del self.conducx[del_idx]; del self.conducy[del_idx]; del self.conducz[del_idx]
                del self.depth_conduc_factor[del_idx]; del self.temp_conduc_factor[del_idx]
                del self.mat[del_idx]; del self.heat_produ[del_idx]; del self.heat_produ_depth[del_idx]
                del self.dens[del_idx]; del self.dens_depth[del_idx]; del self.dens_temp[del_idx]
                del self.therm_expan[del_idx]; del self.body_Type[del_idx]
                del self.x_st_index[del_idx]; del self.y_st_index[del_idx]
                del self.x_end[del_idx]; del self.y_end[del_idx]
                del self.x_st[del_idx]; del self.y_st[del_idx]
                
                # Ensure redundant envelopes are deleted!
                del self.x_envelope[del_idx]; del self.y_envelope[del_idx]
                
                self.n_bod -= 1
                
                # Refresh drawing canvas
                self.ax.cla()
                class DummyEvent: pass
                self.plot_model(DummyEvent())
                print(">>> Selected bodies merged successfully.")
                
            elif self.ans=='no':
                pass
            else:
                print("Entered body number not is valid.")
                self.edit_traits(kind='modal', view='plot_output_view') 
        else: 
            print("Cannot Merge. You have a body added which is not closed yet or you do not have any body in the model.")
            self.edit_traits(kind='modal', view='plot_output_view') 

    def add_body(self,event):
        if self.body_status==1:
            if self.n_bod == 0:
                self.edit_traits(kind='modal', view='topo_view')
                self.edit_traits(kind='modal', view='add_view') 
                if self.material>0 and self.body_type=='normal':
                    self.Temp_density_factor='1'
                    self.ref_density = "3.245E+03"
                    self.depth_density_factor= "0.000E+00"
                    self.Thermal_expansion= "0.000E+00"
                    self.Kx ="2.100E+00"
                    self.Ky ="2.100E+00"
                    self.Kz ="0.000E+00"
                    self.Depth_conduc_factor="1.300E+02"
                    self.Temp_conduc_factor="1.500E-03" 
                    self.radiogenic_heat_producion="2.000E-08"
                    self.depth_heat_factor="0.000E+00"
                else:
                    self.Temp_density_factor='0'
                    self.edit_traits(kind='modal', view='add_view_no_mat')                            
            else:
                self.topo_response="no"
                self.edit_traits(kind='modal', view='add_view') 
                if self.material>0:
                    self.Temp_density_factor='1'
                    self.ref_density = "3.245E+03"
                    self.depth_density_factor= "0.000E+00"
                    self.Thermal_expansion= "0.000E+00"
                    self.Kx ="2.100E+00"
                    self.Ky ="2.100E+00"
                    self.Kz ="0.000E+00"
                    self.Depth_conduc_factor="1.300E+02"
                    self.Temp_conduc_factor="1.500E-03" 
                    self.radiogenic_heat_producion="2.000E-08"
                    self.depth_heat_factor="0.000E+00"
                else:
                    self.Temp_density_factor='0'
                    self.edit_traits(kind='modal', view='add_view_no_mat')
            self.body_status=0;
            
            if self.body_type=="normal":
                self.x_st.append([])
                self.y_st.append([])
                self.x_end.append([])
                self.y_end.append([])
                self.l.append([])
                self.b.append([])
                self.n.append([])
                self.c.append([])
                self.conducx.append([])
                self.conducy.append([])
                self.conducz.append([])
                self.temp_conduc_factor.append([])
                self.depth_conduc_factor.append([])
                self.mat.append([])
                self.heat_produ.append([])
                self.heat_produ_depth.append([])
                self.dens.append([])
                self.dens_depth.append([])
                self.dens_temp.append([])
                self.therm_expan.append([])
                self.n_bod=self.n_bod+1
                
                temp_d=[abs(float(str(self.y_s))-self.Z[i]) for i in range(len(self.Z))]
                temp_index=temp_d.index(min(temp_d))
                
                self.n[self.n_bod-1]=self.name
                self.c[self.n_bod-1]=self.color_edi.name()
                self.conducx[self.n_bod-1]=self.Kx
                self.conducy[self.n_bod-1]=self.Ky
                self.conducz[self.n_bod-1]=self.Kz
                self.depth_conduc_factor[self.n_bod-1]=self.Depth_conduc_factor
                self.temp_conduc_factor[self.n_bod-1]=self.Temp_conduc_factor
                self.mat[self.n_bod-1]=self.material
                self.heat_produ[self.n_bod-1]=self.radiogenic_heat_producion
                self.heat_produ_depth[self.n_bod-1]=self.depth_heat_factor
                self.dens[self.n_bod-1]=self.ref_density
                self.dens_depth[self.n_bod-1]=self.depth_density_factor
                self.dens_temp[self.n_bod-1]=self.Temp_density_factor
                self.therm_expan[self.n_bod-1]=self.Thermal_expansion
                self.body_Type.append(self.body_type)
                
                if self.n_bod==1 and self.topo_response=="no":
                    self.x_envelope=[]
                    self.y_envelope=[]
                    shape_=[]
                    
                    for line in range(int(self.X_S),int(self.X_E)+int(self.reso),int(self.reso)):
                        self.ax.plot(line,0.0,color="red",lw=0.1,marker='o',markersize=5)
                        shape_.append((line,0.0))
                    shape_.append((self.X_E,-400))
                    shape_.append((self.X_S,-400))
                    self.shape=Shape(shape_)
                    self.x_envelope.append([])
                    self.y_envelope.append([])
                    self.x_envelope[0],self.y_envelope[0]=self.shape.exterior.xy
                    
            elif self.body_type =="anomaly": 
                self.edit_traits(kind='modal', view='anomaly_select')
                if self.anomaly_type != "File":
                    self.edit_traits(kind='modal', view='anomaly_view')
                    self.body_status=0;
                    self.x=float(str(self.x_s))
                    temp_d=[abs(float(str(self.y_s))-self.Z[i]) for i in range(len(self.Z))]
                    temp_index=temp_d.index(min(temp_d))
                    self.anomaly=1
                    self.y=self.Z[temp_index]
                    self.anomaly_conut=self.anomaly_conut+1
                    self.anomaly_x.append([])
                    self.anomaly_y.append([])
                    self.anomaly_Type.append([])
                    self.anomaly_amount.append([])
                    self.anomaly_compo.append([])
                    self.anomaly_color.append([])
                    self.body_Type.append(self.body_type)
                    self.anomaly_amount[self.anomaly_conut-1]=self.anomaly_value
                    self.anomaly_Type[self.anomaly_conut-1]=self.anomaly_type
                    self.anomaly_color[self.anomaly_conut-1]=self.color_edi.name()
                    if self.anomaly_Type[self.anomaly_conut-1] not in ("Thermal", "Compositional"):
                        self.anomaly_compo[self.anomaly_conut-1]=self.material
                    else:
                        self.anomaly_compo[self.anomaly_conut-1]=self.anomaly_mat
                else :
                    self.edit_traits(kind='modal', view='anomaly_file_add_view')
                    self.body_status=1
                    
            if self.topo_response== "yes":
                shape_=[]
                f=open(self.topo,"r")
                data=f.readlines()
                for line in data:
                    t=line.split()
                    self.ax.plot(float(str(t[0])),float(str(t[1]))/1000.0,color="red",lw=0.1,marker='o',markersize=5)
                    shape_.append((float(str(t[0])),float(str(t[1]))/1000.0))
                shape_.append((self.X_E,-400))
                shape_.append((self.X_S,-400))
                self.shape=[]
                self.shape=Shape(shape_)
                self.x_envelope=[]
                self.y_envelope=[]
                self.x_envelope.append([])
                self.y_envelope.append([])
                self.x_envelope[0],self.y_envelope[0]=self.shape.exterior.xy
                
                self.ax.set_xlim((self.X_S, self.X_E))
                self.ax.set_ylim((-self.Y,10))
                self.ax.tick_params(labelsize=10)
                self.ax.grid(True)
                plt.draw()
            else:
                pass
            print("You have selected following color")
            print(self.color_edi.name())
            print("This is the coverted color")
            print(str(self.c))
            self.vert = []
            print("addig this body")
            print("Total number of bodies")
            print(len(self.l))
            print("body number count")
            print(self.n_bod)
            print("Total number of envelpos")
            print(len(self.x_envelope))
        else:
            print("You cannot add next body untill you close the current added body.")        
            self.edit_traits(kind='modal', view='plot_output_view') 

    def delete_body(self,event):
        if self.n_bod > 0:
            # Use TraitsUI to create a body-selection dialog.
            from traits.api import HasTraits, Int
            from traitsui.api import View, Item
            
            class DeleteTarget(HasTraits):
                body_no = Int()
                view = View(
                    Item('body_no', label=f'Enter Body Number to Delete (1 to {self.n_bod}):'),
                    title='Targeted Delete',
                    buttons=['OK', 'Cancel']
                )
                
            # Select the last body by default.
            dlg = DeleteTarget(body_no=self.n_bod)
            res = dlg.edit_traits(kind='modal')
            
            if res.result:  # Proceed only if the dialog is accepted.
                idx = dlg.body_no - 1
                if 0 <= idx < self.n_bod:
                    # ==============================================================
                    # Remove all arrays associated with the selected body.
                    # Subsequent body indices shift down after deletion.
                    # ==============================================================
                    del self.l[idx]; del self.b[idx]
                    del self.n[idx]; del self.c[idx]
                    del self.conducx[idx]; del self.conducy[idx]; del self.conducz[idx]
                    del self.depth_conduc_factor[idx]; del self.temp_conduc_factor[idx]
                    del self.mat[idx]; del self.heat_produ[idx]; del self.heat_produ_depth[idx]
                    del self.dens[idx]; del self.dens_depth[idx]; del self.dens_temp[idx]
                    del self.therm_expan[idx]; del self.body_Type[idx]
                    
                    try:
                        del self.x_st_index[idx]; del self.y_st_index[idx]
                        del self.x_end[idx]; del self.y_end[idx]
                        del self.x_st[idx]; del self.y_st[idx]
                    except Exception:
                        pass
                        
                    del self.x_envelope[idx]; del self.y_envelope[idx]
                    
                    self.n_bod -= 1
                    
                    # Clear the axes and redraw the updated model.
                    self.ax.cla()
                    class DummyEvent: pass
                    self.plot_model(DummyEvent())
                    print(f">>> Deleted body {dlg.body_no}; lower body indices were compacted.")
                else:
                    print(">>> Invalid body number; delete cancelled.")
        else:
            print(">>> No bodies available to delete.")

    def edit_body(self,event):
        if self.body_status==1:
            if self.n_bod > 0 :
                print("in edit")
                self.edit_traits(kind='modal', view='body_no_view')
                n_b=int(self.Body_number)
                if n_b <= self.n_bod and self.body_type=='normal' :
                    if self.mat[n_b-1]==0:
                        self.name = self.n[n_b-1]
                        self.color_edi = self.c[n_b-1]
                        self.body_type=self.body_Type[n_b-1]
                        self.ref_density = self.dens[n_b-1]
                        self.depth_density_factor= self.dens_depth[n_b-1]
                        self.Temp_density_factor=self.dens_temp[n_b-1]
                        self.Thermal_expansion=self.therm_expan[n_b-1]
                        self.material = self.mat[n_b-1]
                        self.Kx =self.conducx[n_b-1]
                        self.Ky =self.conducy[n_b-1]
                        self.Kz =self.conducz[n_b-1]
                        self.Depth_conduc_factor=self.depth_conduc_factor[n_b-1]
                        self.Temp_conduc_factor=self.temp_conduc_factor[n_b-1]
                        self.radiogenic_heat_producion =self.heat_produ[n_b-1]
                        self.depth_heat_factor=self.heat_produ_depth[n_b-1]
                        self.edit_traits(kind='modal', view='edit_view')
                        print(self.l[n_b-1])
                        print(self.b[n_b-1])
                        self.n[n_b-1]=self.name
                        self.c[n_b-1]=self.color_edi.name()
                        self.conducx[n_b-1]=self.Kx
                        self.conducy[n_b-1]=self.Ky
                        self.conducz[n_b-1]=self.Kz
                        self.depth_conduc_factor[n_b-1]=self.Depth_conduc_factor
                        self.temp_conduc_factor[n_b-1]=self.Temp_conduc_factor
                        self.mat[n_b-1]=self.material
                        self.heat_produ[n_b-1]=self.radiogenic_heat_producion
                        self.heat_produ_depth[n_b-1]=self.depth_heat_factor
                        self.dens[n_b-1]=self.ref_density
                        self.dens_depth[n_b-1]=self.depth_density_factor
                        self.dens_temp[n_b-1]=self.Temp_density_factor
                        self.therm_expan[n_b-1]=self.Thermal_expansion
                        self.body_Type[n_b-1]=self.body_type
                        self.ax.plot(self.l[n_b-1],self.b[n_b-1],color="green",lw=1,marker='o')
                        self.ax.tick_params(labelsize=10)
                        poly = Polygon(list(zip(self.l[n_b-1],self.b[n_b-1])),facecolor=self.c[n_b-1],label=self.mat[n_b-1], lw=1)
                        self.ax.add_patch(poly)
                        poly=[]
                        plt.draw()
                    elif self.mat[n_b-1] > 0 :
                        self.name = self.n[n_b-1]
                        self.color_edi = self.c[n_b-1]
                        self.material = self.mat[n_b-1]
                        self.edit_traits(kind='modal', view='edit_view_mantle')
                        print(self.l[n_b-1])
                        print(self.b[n_b-1])
                        self.n[n_b-1]=self.name
                        self.c[n_b-1]=self.color_edi.name()
                        self.mat[n_b-1]=self.material
                        self.ax.plot(self.l[n_b-1],self.b[n_b-1],color="green",lw=1,marker='o')
                        poly = Polygon(list(zip(self.l[n_b-1],self.b[n_b-1])),facecolor=self.c[n_b-1],label=self.mat[n_b-1], lw=1)
                        self.ax.add_patch(poly)
                        self.ax.tick_params(labelsize=10)
                        poly=[]
                        plt.draw()
                elif n_b <= len(self.anomaly_compo) and self.body_type =="anomaly":
                    print(n_b)
                    self.anomaly_mat=int(self.anomaly_compo[n_b-1])
                    self.anomaly_value=float(self.anomaly_amount[n_b-1])
                    self.anomaly_type=self.anomaly_Type[n_b-1]
                    self.color_edi = self.anomaly_color[n_b-1]
                    self.edit_traits(kind='modal', view='anomaly_edit_view')
                    self.anomaly_compo[n_b-1]=self.anomaly_mat
                    self.anomaly_amount[n_b-1]=self.anomaly_value
                    self.anomaly_Type[n_b-1]=self.anomaly_type
                    self.anomaly_color[n_b-1]=self.color_edi.name()
                else:
                    print("Desired body number exceeds the total number of bodies")
                    self.edit_traits(kind='modal', view='plot_output_view')
            else:
                print("You have no body to edit. First add some.")
                self.edit_traits(kind='modal', view='plot_output_view') 
        else:
            print("You cannot edit bodies yet. You have a body added which is not closed yet.")
            self.edit_traits(kind='modal', view='plot_output_view') 

    def edit_shape(self,event):
        if self.body_status==1 and self.n_bod>=1:
            if (self.mat[self.n_bod-1]!=99 ): 
                self.edit_traits(kind='modal', view='shape_alert')
            else:
                self.edit_traits(kind='modal', view='plot_output_cosmectics')
                if self.anomaly_conut== 0:
                    print("in shape edit")
                    fig_edit = plt.figure()
                    fig_edit.canvas.manager.set_window_title('LitMod Normal Body shape edit')
                    ax_edit = fig_edit.add_subplot(111)
                    ax_edit.set_xlim((self.X_S, self.X_E))
                    ax_edit.set_ylim((-self.max_depth,10))
                    ax_edit.tick_params(labelsize=10)
                    try:
                        f=np.loadtxt(self.digitized,comments=">")
                        ax_edit.plot(f[:,0],f[:,1],'*',lw=1,color='grey')
                    except:
                        pass
                    try:
                        f=np.loadtxt(self.topo)
                        ax_edit.plot(f[:,0],f[:,1]/1000.0,color="black",lw=0.1,marker='d',markersize=self.marker_size)
                    except:
                        pass
                    ax_edit.grid(True)
                    for i in range(len(self.l)):
                        self.l_text=min(self.l[i][:])+ (max(self.l[i][:])-min(self.l[i][:]))/2.0
                        self.b_text=min(self.b[i][:])+ (max(self.b[i][:])-min(self.b[i][:]))/2.0
                        ax_edit.text(self.l_text, self.b_text, str(int(i+1))+" "+str(self.n[i])+" Material "
                            +str(self.mat[i]), fontsize=self.label_size,fontstyle='italic',color='black',bbox=dict(facecolor='white', alpha=0.9))
                    l_edit=[];b_edit=[];no_nodes=[]
                    cursor = Cursor(ax_edit, useblit=True, color='black', linewidth=1)
                    for i in range(len(self.x_envelope)-1):
                        no_nodes.append([])
                        no_nodes[i]=len(self.x_envelope[i])
                    print("printing no of nodes")
                    print(no_nodes)
                    for i in range(len(self.x_envelope)-1):
                        l_edit=np.append(l_edit,self.x_envelope[i])
                        b_edit=np.append(b_edit,self.y_envelope[i])
                    edit_poly = Polygon(list(zip(l_edit, b_edit)), fill=False, edgecolor='none', linewidth=1.5)# change: blue to none 
                    ax_edit.add_patch(edit_poly)
                    
                    self.shape_cursor = cursor
                    self.poly_interactor = edit.PolygonInteractor(ax_edit,edit_poly,no_nodes,self.digitized,self.reso,self.Z,self.topo,"normal_body_shape.dat")
                    
                    no_nodes=[]
                    edit_poly=[]
                    plt.show()
                else:
                    print("in shape edit")
                    fig_edit = plt.figure()
                    fig_edit.canvas.manager.set_window_title('LitMod Normal Body shape edit')
                    ax_edit = fig_edit.add_subplot(111)
                    ax_edit.set_xlim((self.X_S, self.X_E))
                    ax_edit.set_ylim((-self.max_depth,10))
                    ax_edit.tick_params(labelsize=10)
                    try:
                        f=np.loadtxt(self.digitized,comments=">")
                        ax_edit.plot(f[:,0],f[:,1],'*',lw=1,color='grey')
                    except:
                        pass
                    try:
                        f=np.loadtxt(self.topo)
                        ax_edit.plot(f[:,0],f[:,1]/1000.0,color="black",lw=0.1,marker='d',markersize=5)
                    except:
                        pass
                    ax_edit.grid(True)
                    for i in range(len(self.l)):
                        self.l_text=min(self.l[i][:])+ (max(self.l[i][:])-min(self.l[i][:]))/2.0
                        self.b_text=min(self.b[i][:])+ (max(self.b[i][:])-min(self.b[i][:]))/2.0
                        ax_edit.text(self.l_text, self.b_text, str(int(i+1))+" "+str(self.n[i])+" Material "
                            +str(self.mat[i]), fontsize=8,fontstyle='italic',color='black',bbox=dict(facecolor='white', alpha=0.9))
                    l_edit=[];b_edit=[];no_nodes=[]
                    cursor = Cursor(ax_edit, useblit=False, color='black', linewidth=1)
                    for i in range(len(self.x_envelope)-1):
                        no_nodes.append([])
                        no_nodes[i]=len(self.x_envelope[i])
                    print("printing no of nodes")
                    print(no_nodes)
                    for i in range(len(self.x_envelope)-1):
                        l_edit=np.append(l_edit,self.x_envelope[i])
                        b_edit=np.append(b_edit,self.y_envelope[i])
                    edit_poly = Polygon(list(zip(l_edit, b_edit)), fill=False, edgecolor='none', linewidth=1.5) # change: blue to none 
                    
                    print("in shape edit (Anomaly)")
                    fig_edit_anom = plt.figure()
                    fig_edit_anom.canvas.manager.set_window_title('LitMod Anomaly Body shape edit')
                    ax_edit_anom = fig_edit_anom.add_subplot(111)
                    ax_edit_anom.set_xlim((self.X_S, self.X_E))
                    ax_edit_anom.set_ylim((-self.max_depth,10))
                    ax_edit_anom.tick_params(labelsize=10)
                    try:
                        f=np.loadtxt(self.digitized,comments=">")
                        ax_edit_anom.plot(f[:,0],f[:,1],'*',lw=1,color='grey')
                    except:
                        pass
                    try:
                        f=np.loadtxt(self.topo)
                        ax_edit_anom.plot(f[:,0],f[:,1]/1000.0,color="black",lw=0.1,marker='d',markersize=5)
                    except:
                        pass
                    ax_edit_anom.grid(True)
                    for i in range(self.anomaly_conut):
                        self.l_text_anom=min(self.anomaly_x[i][:])+ (max(self.anomaly_x[i][:])-min(self.anomaly_x[i][:]))/2.0
                        self.b_text_anom=min(self.anomaly_y[i][:])+ (max(self.anomaly_y[i][:])-min(self.anomaly_y[i][:]))/2.0
                        ax_edit_anom.text(self.l_text_anom, self.b_text_anom, str(int(i+1))+" "+str(self.anomaly_type[i])+" Material "
                            +str(self.anomaly_amount[i]), fontsize=8,fontstyle='italic',color='black',bbox=dict(facecolor='white', alpha=0.9))
                    l_edit_anom=[];b_edit_anom=[];no_nodes_anom=[]
                    cursor_ = Cursor(ax_edit_anom, useblit=True, color='black', linewidth=1)
                    for i in range(self.anomaly_conut):
                        no_nodes_anom.append([])
                        no_nodes_anom[i]=len(self.anomaly_x[i])
                    print("printing no of nodes")
                    print(no_nodes_anom)
                    for i in range(self.anomaly_conut):
                        l_edit_anom=np.append(l_edit_anom,self.anomaly_x[i])
                        b_edit_anom=np.append(b_edit_anom,self.anomaly_y[i])
                    edit_poly_anom = Polygon(list(zip(l_edit_anom, b_edit_anom)), fill=False, edgecolor='none', linewidth=1.5) #change:red to none
                    ax_edit_anom.add_patch(edit_poly_anom)
                    ax_edit_anom.plot(l_edit,b_edit,marker='o',markersize=3, markerfacecolor='black', linestyle='none')
                    ax_edit.add_patch(edit_poly)
                    ax_edit.plot(l_edit_anom,b_edit_anom,marker='o',markersize=3, markerfacecolor='black', linestyle='none')
                    ax_edit.set_ylim((-self.max_depth,10))
                    
                    self.shape_cursor = cursor
                    self.shape_cursor_anom = cursor_
                    self.poly_interactor = edit.PolygonInteractor(ax_edit,edit_poly,no_nodes,self.digitized,self.reso,self.Z,self.topo,"normal_body_shape.dat")
                    self.poly_interactor_anom = edit.PolygonInteractor(ax_edit_anom,edit_poly_anom,no_nodes_anom,self.digitized,self.reso,self.Z,self.topo,"anomaly_body_shape.dat")
                    
                    no_nodes_anom=[]
                    edit_poly_anom=[]
                    plt.show()
        else:
            print("You cannot change the nodes yet. Model is not closed yet")
            self.edit_traits(kind='modal', view='plot_output_view') 

    #################################################################################
#This function takes care of the closing the body. 
#################################################################################
    def _close_polygon(self,event):
        if event.button==3:
            x = [self.vert[k][0] for k in range(len(self.vert))]
            y = [self.vert[k][1] for k in range(len(self.vert))]
            if self.body_type=="normal":
                self.l[self.n_bod-1]=x
                self.b[self.n_bod-1]=y
                self.x_st[self.n_bod-1]=self.l[self.n_bod-1][0]
                self.y_st[self.n_bod-1]=self.b[self.n_bod-1][0]
                self.x_end[self.n_bod-1]=self.l[self.n_bod-1][len(self.l[self.n_bod-1])-1]
                self.y_end[self.n_bod-1]=self.b[self.n_bod-1][len(self.b[self.n_bod-1])-1]
            elif self.body_type=="thermal_anomaly" or self.body_type=="seismic_anomaly" or self.body_type=="anomaly":
                self.anomaly_x[self.anomaly_conut-1]=x
                self.anomaly_y[self.anomaly_conut-1]=y              
            plt.draw()
            
            # ===================================================================
            # Split operation logic.
            # ===================================================================
            if self.body_type=="split":
                try:
                        xa=[];ya=[]
                        a=Shape()
                        xa,ya,a,err=shape_split.shape_split(self.vert[:][:],self.l_split,self.b_split,self.l[self.split_body_no-1],self.b[self.split_body_no-1])
                        
                        self.l.insert(self.split_body_no,xa)
                        self.b.insert(self.split_body_no,ya)
                        
                        self.n_bod=self.n_bod+1
                        
                        # Insert the new properties and envelope at the selected indices.
                        self.n.insert(self.split_body_no-1,self.name)
                        self.c.insert(self.split_body_no-1,self.color_edi.name())
                        self.conducx.insert(self.split_body_no-1,self.Kx)
                        self.conducy.insert(self.split_body_no-1,self.Ky)
                        self.conducz.insert(self.split_body_no-1,self.Kz)
                        self.depth_conduc_factor.insert(self.split_body_no-1,self.Depth_conduc_factor)
                        self.temp_conduc_factor.insert(self.split_body_no-1,self.Temp_conduc_factor)
                        self.mat.insert(self.split_body_no-1,self.material)
                        self.heat_produ.insert(self.split_body_no-1,self.radiogenic_heat_producion)
                        self.heat_produ_depth.insert(self.split_body_no-1,self.depth_heat_factor)
                        self.dens.insert(self.split_body_no-1,self.ref_density)
                        self.dens_depth.insert(self.split_body_no-1,self.depth_density_factor)
                        self.dens_temp.insert(self.split_body_no-1,self.Temp_density_factor)
                        self.therm_expan.insert(self.split_body_no-1,self.Thermal_expansion)
                        self.body_Type.insert(self.split_body_no-1,'normal')
                        
                        self.x_st_index.append([])
                        self.x_st_index[self.split_body_no-1]=self.l[self.split_body_no-1][0]
                        self.y_st_index.append([])
                        self.y_st_index[self.split_body_no-1]=self.b[self.split_body_no-1][0]
                        
                        # Interpolate the new envelope from the selected body.
                        temp_2_shape=Shape([(self.x_envelope[self.split_body_no-1][k],self.y_envelope[self.split_body_no-1][k]) for k in range(len(self.x_envelope[self.split_body_no-1]))])
                        
                        # Use buffer(0) to repair minor geometry issues.
                        a = a.buffer(0)
                        temp_2_shape = temp_2_shape.buffer(0)
                        b = temp_2_shape.difference(a)
                        if b.geom_type == 'MultiPolygon':
                            b = max(b.geoms, key=lambda p: p.area)
                            
                        env_xa, env_ya = b.exterior.coords.xy 
                        self.x_envelope.insert(self.split_body_no, list(env_xa))
                        self.y_envelope.insert(self.split_body_no, list(env_ya))
                        
                        self.body_type='normal'
                        self.body_status=1
                        self.vert = []
                        self.path = []
                        
                        # Clear and redraw the model.
                        self.ax.cla()
                        class DummyEvent: pass
                        self.plot_model(DummyEvent())
                        
                        # Return after handling the split; do not fall through to new-body drawing.
                        return 
                        
                except Exception:
                    try:
                        # Handle the split in the reverse direction.
                        xa=[];ya=[]
                        a=Shape()
                        xa,ya,a,err=shape_split.shape_split(self.vert[::-1][:],self.l_split,self.b_split,self.l[self.split_body_no-1],self.b[self.split_body_no-1])
                        
                        self.l.insert(self.split_body_no,xa)
                        self.b.insert(self.split_body_no,ya)
                        self.n_bod=self.n_bod+1
                        
                        self.n.insert(self.split_body_no-1,self.name)
                        self.c.insert(self.split_body_no-1,self.color_edi.name())
                        self.conducx.insert(self.split_body_no-1,self.Kx)
                        self.conducy.insert(self.split_body_no-1,self.Ky)
                        self.conducz.insert(self.split_body_no-1,self.Kz)
                        self.depth_conduc_factor.insert(self.split_body_no-1,self.Depth_conduc_factor)
                        self.temp_conduc_factor.insert(self.split_body_no-1,self.Temp_conduc_factor)
                        self.mat.insert(self.split_body_no-1,self.material)
                        self.heat_produ.insert(self.split_body_no-1,self.radiogenic_heat_producion)
                        self.heat_produ_depth.insert(self.split_body_no-1,self.depth_heat_factor)
                        self.dens.insert(self.split_body_no-1,self.ref_density)
                        self.dens_depth.insert(self.split_body_no-1,self.depth_density_factor)
                        self.dens_temp.insert(self.split_body_no-1,self.Temp_density_factor)
                        self.therm_expan.insert(self.split_body_no-1,self.Thermal_expansion)
                        self.body_Type.insert(self.split_body_no-1,'normal')
                        
                        self.x_st_index.append([])
                        self.x_st_index[self.split_body_no-1]=self.l[self.split_body_no-1][0]
                        self.y_st_index.append([])
                        self.y_st_index[self.split_body_no-1]=self.b[self.split_body_no-1][0]
                        
                        temp_2_shape=Shape([(self.x_envelope[self.split_body_no-1][k],self.y_envelope[self.split_body_no-1][k]) for k in range(len(self.x_envelope[self.split_body_no-1]))])
                        
                        # Use buffer(0) to repair minor geometry issues.
                        a = a.buffer(0)
                        temp_2_shape = temp_2_shape.buffer(0)
                        b = temp_2_shape.difference(a)
                        if b.geom_type == 'MultiPolygon':
                            b = max(b.geoms, key=lambda p: p.area)
                            
                        env_xa, env_ya = b.exterior.coords.xy 
                        self.x_envelope.insert(self.split_body_no, list(env_xa))
                        self.y_envelope.insert(self.split_body_no, list(env_ya))
                        
                        self.body_type='normal'
                        self.body_status=1
                        self.vert = []
                        self.path = []
                        
                        self.ax.cla()
                        class DummyEvent: pass
                        self.plot_model(DummyEvent())
                        
                        # Return after handling the reverse split.
                        return 
                        
                    except Exception:                    
                        self.edit_traits(kind='modal', view='split_body_alert')
                        self.body_type='normal'
                        self.body_status=1;
                        return 
            # ===================================================================

            # ===================================================================
            # Manual new-body drawing logic.
            # ===================================================================
            if self.body_type=="normal":
                self.shape=Shape( [ (self.x_envelope[self.n_bod-1][k],self.y_envelope[self.n_bod-1][k]) for k in range(len(self.x_envelope[self.n_bod-1]))] )
                t_x,t_y=self.shape.exterior.xy
                st_index=[]
                end_index=[]
                pol_ext = shapely.geometry.LinearRing(self.shape.exterior.coords)
                ###################################################################################33
                ############fixing problem start point
                print ("Working on Body number")
                print (self.n_bod)
                print ("start @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
                ### finding closest node in the envelop
                temp_d=[ math.sqrt((t_x[i] - self.x_st[self.n_bod-1])**2 + (t_y[i] - self.y_st[self.n_bod-1])**2) for i in range(len(t_x))]
                st_index=temp_d.index(min(temp_d))
                st_x=t_x[st_index]
                st_y=t_y[st_index]
                print("Start Index:");print(st_index)
                print ("Start node");print(st_x,st_y)
                ###################################################################################33
                ############fixing end point
                print ("end #################################################################################")
                temp_d=[ math.sqrt((t_x[i] - self.x_end[self.n_bod-1])**2 + (t_y[i] - self.y_end[self.n_bod-1])**2) for i in range(len(t_x))]
                end_index=temp_d.index(min(temp_d))
                end_x=t_x[end_index]
                end_y=t_y[end_index]
                print("End Index:");print(end_index)
                print ("End node");print(end_x,end_y)
                ## subbody polygon starts here
                temp=[] ## subbody
                self.vert.append((self.x_end[self.n_bod-1],self.y_end[self.n_bod-1]))
                self.x_end[self.n_bod-1]=end_x
                self.y_end[self.n_bod-1]=end_y
                self.x_st[self.n_bod-1]=st_x
                self.y_st[self.n_bod-1]=st_y
                print("Fixing start")
                if (st_index==0):
                    pass
                else:
                    it=0                    
                    try:    
                        while (it < st_index):
                            temp.append((t_x[it],t_y[it]))
                            print(it)
                            it=it+1
                        temp.append((t_x[it],t_y[it]))
                    except IndexError:
                        pass
                for i in range(len(self.vert)):
                    temp.append(self.vert[i])
                print("Fixing End")
                if (t_x[end_index] == self.X_E):
                    pass
                else: 
                    it=end_index;
                    while (t_x[it] < self.X_E):
                        temp.append((t_x[it],t_y[it]))
                        print(it)
                        it=it+1
                    temp.append((t_x[it],t_y[it]))
                temp.append((self.X_E,-400))
                temp.append((self.X_S,-400))
                ######
                self.x_st_index.append([])
                self.x_st_index[self.n_bod-1]=t_x[st_index]
                self.y_st_index.append([])
                self.y_st_index[self.n_bod-1]=t_y[st_index]
                q_x,q_y=self.shape.exterior.xy
                print("here is the problme _ master polygon")
                body_shape=Shape(temp[:][:])
                print(temp)
                print("here is the problem")
                try:
                    self.x_envelope.append([])
                    self.y_envelope.append([])
                    
                    # Repair minor polygon defects in newly drawn bodies.
                    self.shape = self.shape.buffer(0)
                    body_shape = body_shape.buffer(0)
                    a = self.shape.difference(body_shape)
                    if a.geom_type == 'MultiPolygon':
                        a = max(a.geoms, key=lambda p: p.area)
                        
                    self.vert=[]
                    xa,ya=a.exterior.xy
                    self.shape=body_shape
                    self.x_envelope[self.n_bod],self.y_envelope[self.n_bod]=self.shape.exterior.xy
                    for i in range(len(xa)):
                        self.vert.append((xa[i],ya[i]))
                    x = [self.vert[k][0] for k in range(len(self.vert))]
                    y = [self.vert[k][1] for k in range(len(self.vert))]
                    self.path.set_data(x,y)
                    self.l[self.n_bod-1]=x
                    self.b[self.n_bod-1]=y
                except Exception:
                        self.edit_traits(kind='modal', view='body_alert')    
            elif self.body_type=="anomaly":
                self.vert.append((self.vert[0][0],self.vert[0][1]))
                x = [self.vert[k][0] for k in range(len(self.vert))]
                y = [self.vert[k][1] for k in range(len(self.vert))]
                self.path.set_data(x,y)
                self.anomaly_x[self.anomaly_conut-1]=x
                self.anomaly_y[self.anomaly_conut-1]=y
                print("print anomaly coordinates")
                print(self.anomaly_x)
                print(self.anomaly_y)
                self.body_type='normal'

            self.body_status=1;    
            print ("after closing current body")
            print ("Total number of bodies")
            print (len(self.l))
            print ("Body count number")
            print (self.n_bod)
            print ("Total number of envelpos")
            print(len(self.x_envelope))

    def Close_Model(self,event):
        if self.body_status==1:
            if self.n_bod>0:
                self.x_st.append([])
                self.y_st.append([])
                self.x_end.append([])
                self.y_end.append([])
                self.l.append([])
                self.b.append([])
                self.n.append([])
                self.c.append([])
                self.conducx.append([])
                self.conducy.append([])
                self.conducz.append([])
                self.temp_conduc_factor.append([])
                self.depth_conduc_factor.append([])
                self.mat.append([])
                self.heat_produ.append([])
                self.heat_produ_depth.append([])
                self.dens.append([])
                self.dens_depth.append([])
                self.dens_temp.append([])
                self.therm_expan.append([])
                self.body_Type.append([])
                self.n_bod=self.n_bod+1
                
                self.x_st[self.n_bod-1]=self.X_S
                self.y_st[self.n_bod-1]=float(str("-400.00"))
                self.x_end[self.n_bod-1]=self.X_E
                self.y_end[self.n_bod-1]=float(str("-400.00"))
                self.n[self.n_bod-1]="asthenosphere"
                self.c[self.n_bod-1]='#800080' #"purple"
                self.conducx[self.n_bod-1]="6.600E+08"
                self.conducy[self.n_bod-1]="6.600E+08"
                self.conducz[self.n_bod-1]="0.000E+00"
                self.depth_conduc_factor[self.n_bod-1]="0.000E+00"
                self.temp_conduc_factor[self.n_bod-1]="0.000E+00"
                self.mat[self.n_bod-1]=99
                self.heat_produ[self.n_bod-1]="0.000E-08"
                self.heat_produ_depth[self.n_bod-1]="0.000E+00"
                self.dens[self.n_bod-1]="3.603E+03"
                self.dens_depth[self.n_bod-1]="0.000E+00"
                self.dens_temp[self.n_bod-1]="1"
                self.therm_expan[self.n_bod-1]="0.000E+00"
                self.body_Type[self.n_bod-1]='normal'

                # ==============================================================
                # [FIX] Forcefully clean up residual empty lists at the end.
                # Directly copy the data from the previous envelope and append it.
                # This ensures the lengths of x_envelope/y_envelope and n_bod 
                # remain strictly synchronized at all times
                # ==============================================================
                self.x_envelope = self.x_envelope[:self.n_bod]
                self.y_envelope = self.y_envelope[:self.n_bod]
                
                self.x_envelope.append(list(self.x_envelope[self.n_bod-1]))
                self.y_envelope.append(list(self.y_envelope[self.n_bod-1]))
                # ==============================================================

                #===============================================================
                # get envelope point
                points = [(self.x_envelope[self.n_bod-1][k], self.y_envelope[self.n_bod-1][k]) for k in range(len(self.x_envelope[self.n_bod-1]))]
                
                # delet_point
                cleaned = []
                for p in points:
                    if not cleaned or p != cleaned[-1]:
                        cleaned.append(p)
                        
                # creat_pology
                self.shape = Shape(cleaned)
                xa, ya = self.shape.exterior.xy
                
                # alter:clear exterior 
                self.path.set_data(xa, ya)
                plt.draw()
                self.l[self.n_bod-1] = xa
                self.b[self.n_bod-1] = ya
                #===============================================================

                self.x_st_index.append([])
                self.x_st_index[self.n_bod-1]=0.0
                self.y_st_index.append([])
                self.y_st_index[self.n_bod-1]=float(str("-400.00"))
                
                print("after closing current body")
                print("Total number of bodies")
                print(len(self.l))
                print("Body count number")
                print(self.n_bod)
                print("Total number of envelpos")
                print(len(self.x_envelope))
                print("Last two envelops##########################################")
                
                print(self.x_envelope[-1], self.y_envelope[-1])
                print(self.x_envelope[-2], self.y_envelope[-2])
            else:
                print("Can not close model. First add some bodies.")
                self.edit_traits(kind='modal', view='plot_output_view') 
        else:
            print("Cannot close model yet. You have a body added which is not closed yet.")
            self.edit_traits(kind='modal', view='plot_output_view')

##########################################################################################
############### main function for building a model
##########################################################################################
def build_model(litmod_path):
        class Model_setup(HasTraits):
            X_S = 0
            X_E = 1000
            Y_length = 400
            X_resolution = 10
            dir = Directory
            view1 = View(Group(Item(name = 'X_S',label="Start of the profile (km)"),
                 Item(name = 'X_E',label="End of the profile (km)"),
                 Item(name = 'Y_length',label="Depth of the profile (must not be changed)"),
                 Item(name = 'X_resolution',label="Resolution along length of the profile (2,5,10 km)"),
                 Item(name = 'dir',label="Folder"),
                 show_border = True),
                 title="Model Setup",
                 buttons = [OKButton])
        button_color='lightblue'
        button_font=12
        model_setup=Model_setup()
        model_setup.edit_traits(kind='modal')
        try:
            os.chdir(model_setup.dir)
        except Exception:
            model_setup.dir=litmod_path
            os.chdir(model_setup.dir)
        if os.path.exists("bodies_GUI.out")==True:
            os.remove("bodies_GUI.out")
        else:
            pass
        fig = plt.figure(0)
        fig.canvas.manager.set_window_title('LitMod Build Model')
        ax = fig.add_subplot(111)
        cursor = Cursor(ax, useblit=True, color='black', linewidth=1)
        plt.subplots_adjust(bottom=0.2)
        plt.xlabel('Distance (Km)')
        plt.ylabel('Depth (Km)')
        ax.set_xlim((-2, 2))
        ax.set_ylim((-2, 2))
        
        gdl=image.imread(litmod_path+"/Images/gdl.jpg")
        subitop=image.imread(litmod_path+"/Images/subitop.jpg")
        logo_gdl= plt.axes([0.01, 0.92, 0.08, 0.07])
        logo_gdl.imshow(gdl, aspect='auto', extent=(0.2, 0.3, .25, .35), zorder=-1)
        logo_gdl.set_xticks([])
        logo_gdl.set_yticks([])
        logo_subitop= plt.axes([0.1, 0.92, 0.08, 0.07])
        logo_subitop.imshow(subitop, aspect='auto', extent=(0.2, 0.3, .25, .35), zorder=-1)
        logo_subitop.set_xticks([])
        logo_subitop.set_yticks([])
        
        plot_results = plt.axes([0.21, 0.92, 0.1, 0.075])
        load_option = plt.axes([0.32, 0.92, 0.1, 0.075])
        run_option = plt.axes([0.45, 0.92, 0.1, 0.075])
        close_model=  plt.axes([0.60, 0.92, 0.15, 0.075])
        quit_option = plt.axes([0.80, 0.92, 0.1, 0.075])
        split_option = plt.axes([0.1, 0.05, 0.1, 0.075])
        delete_option = plt.axes([0.21, 0.05, 0.1, 0.075])
        edit_shape_option = plt.axes([0.32, 0.05, 0.1, 0.075])
        plot_option = plt.axes([0.43, 0.05, 0.1, 0.075])
        merge_option = plt.axes([0.54, 0.05, 0.1, 0.075])
        save_option = plt.axes([0.65, 0.05, 0.1, 0.075])
        add_option = plt.axes([0.76, 0.05, 0.1, 0.075])
        edit_option = plt.axes([0.87, 0.05, 0.1, 0.075])
        
        fig.cnv = Canvas(litmod_path,fig,ax,save_option,add_option,edit_option,plot_option,run_option,edit_shape_option,delete_option,\
            split_option,load_option,quit_option,close_model,merge_option,plot_results,model_setup.X_S,model_setup.X_E,model_setup.Y_length,model_setup.X_resolution)
        
        fig.bnext = Button(add_option, 'Add Body', color=button_color)
        fig.bnext.on_clicked(fig.cnv.add_body)
        fig.bnext.label.set_fontsize('large')
        
        fig.bprev = Button(save_option, 'Save Model', color=button_color)
        fig.bprev.on_clicked(fig.cnv.save_body)
        fig.bprev.label.set_fontsize('large')
        
        fig.edit = Button(edit_option, 'Edit Property', color=button_color)
        fig.edit.on_clicked(fig.cnv.edit_body)
        fig.edit.label.set_fontsize('large')
        
        fig.plot = Button(plot_option, 'Refresh', color=button_color)
        fig.plot.on_clicked(fig.cnv.plot_model)
        fig.plot.label.set_fontsize('large')

        fig.run = Button(run_option, 'Run Model', color=button_color)
        fig.run.on_clicked(fig.cnv.run_model)
        fig.run.label.set_fontsize('large')

        fig.merge = Button(merge_option, 'Merge Bodies', color=button_color)
        fig.merge.on_clicked(fig.cnv.merge_bodies)
        fig.merge.label.set_fontsize('large')

        fig.shape_edit = Button(edit_shape_option, 'Shape Change', color=button_color)
        fig.shape_edit.on_clicked(fig.cnv.edit_shape)
        fig.shape_edit.label.set_fontsize('large')

        fig.delete = Button(delete_option, 'Delete Body', color=button_color)
        fig.delete.on_clicked(fig.cnv.delete_body)
        fig.delete.label.set_fontsize('large')
        
        fig.split = Button(split_option, 'Split Body', color=button_color)
        fig.split.on_clicked(fig.cnv.split_body)
        fig.split.label.set_fontsize('large')

        fig.load_btn = Button(load_option, 'Load Model', color=button_color)
        fig.load_btn.on_clicked(fig.cnv.read_model)
        fig.load_btn.label.set_fontsize('large')

        fig.results = Button(plot_results, 'Plot Results', color=button_color)
        fig.results.on_clicked(fig.cnv.Plot_results)
        fig.results.label.set_fontsize('large')
        
        fig.quit_btn = Button(quit_option, 'Quit', color=button_color)
        fig.quit_btn.on_clicked(fig.cnv.quit)
        fig.quit_btn.label.set_fontsize('large')
        
        fig.close_btn = Button(close_model, 'Close Model', color=button_color)
        fig.close_btn.on_clicked(fig.cnv.Close_Model)
        fig.close_btn.label.set_fontsize('large')
        
        plt.connect('motion_notify_event', fig.cnv.set_location)
        plt.show()

#################################################################################
#Main function for loading a model
#################################################################################
def load_model(litmod_path):
        class Model_setup(HasTraits):
            X_S = 0
            X_E = 1000
            Y_length = 400
            X_resolution = 10
            dir = Directory
            view1 = View(Group(Item(name = 'X_S',label="Start of the profile (km)"),
                 Item(name = 'X_E',label="End of the profile (km)"),
                 Item(name = 'Y_length',label="Depth of the profile (must not be changed)"),
                 Item(name = 'X_resolution',label="Resolution along length of the profile (5,10 km)"),
                 Item(name = 'dir',label="Folder"),
                 show_border = True),
                 title="Model Setup",
                 buttons = [OKButton])
        button_color='lightblue'
        button_font=12
        model_setup=Model_setup()
        model_setup.edit_traits(kind='modal')
        try:
            os.chdir(model_setup.dir)
        except Exception:
            model_setup.dir=litmod_path
            os.chdir(model_setup.dir)
        fig = plt.figure(0)
        fig.canvas.manager.set_window_title('LitMod Load Model')
        ax = fig.add_subplot(111)
        cursor = Cursor(ax, useblit=True, color='black', linewidth=1)
        plt.subplots_adjust(bottom=0.2)
        plt.xlabel('Distance (Km)')
        plt.ylabel('Depth (Km)')
        ax.set_xlim((-2, 2))
        ax.set_ylim((-2, 2))
        
        gdl=image.imread(litmod_path+"/Images/gdl.jpg")
        subitop=image.imread(litmod_path+"/Images/subitop.jpg")
        logo_gdl= plt.axes([0.01, 0.92, 0.08, 0.07])
        logo_gdl.imshow(gdl, aspect='auto', extent=(0.2, 0.3, .25, .35), zorder=-1)
        logo_gdl.set_xticks([])
        logo_gdl.set_yticks([])
        logo_subitop= plt.axes([0.1, 0.92, 0.08, 0.07])
        logo_subitop.imshow(subitop, aspect='auto', extent=(0.2, 0.3, .25, .35), zorder=-1)
        logo_subitop.set_xticks([])
        logo_subitop.set_yticks([])
        
        plot_results = plt.axes([0.21, 0.92, 0.1, 0.075])
        load_option = plt.axes([0.32, 0.92, 0.1, 0.075])
        run_option = plt.axes([0.45, 0.92, 0.1, 0.075])
        close_model=  plt.axes([0.60, 0.92, 0.15, 0.075])
        quit_option = plt.axes([0.80, 0.92, 0.1, 0.075])
        split_option = plt.axes([0.1, 0.05, 0.1, 0.075])
        delete_option = plt.axes([0.21, 0.05, 0.1, 0.075])
        edit_shape_option = plt.axes([0.32, 0.05, 0.1, 0.075])
        plot_option = plt.axes([0.43, 0.05, 0.1, 0.075])
        merge_option = plt.axes([0.54, 0.05, 0.1, 0.075])
        save_option = plt.axes([0.65, 0.05, 0.1, 0.075])
        add_option = plt.axes([0.76, 0.05, 0.1, 0.075])
        edit_option = plt.axes([0.87, 0.05, 0.1, 0.075])
        
        fig.cnv = Canvas(litmod_path,fig,ax,save_option,add_option,edit_option,plot_option,run_option,edit_shape_option,delete_option,split_option,load_option,quit_option,close_model,merge_option,plot_results,0,1000,400,10)
        
        fig.bnext = Button(add_option, 'Add Body', color=button_color)
        fig.bnext.on_clicked(fig.cnv.add_body)
        fig.bnext.label.set_fontsize('large')
        
        fig.bprev = Button(save_option, 'Save Model', color=button_color)
        fig.bprev.on_clicked(fig.cnv.save_body)
        fig.bprev.label.set_fontsize('large')
        
        fig.edit = Button(edit_option, 'Edit Property', color=button_color)
        fig.edit.on_clicked(fig.cnv.edit_body)
        fig.edit.label.set_fontsize('large')
        
        fig.plot = Button(plot_option, 'Refresh', color=button_color)
        fig.plot.on_clicked(fig.cnv.plot_model)
        fig.plot.label.set_fontsize('large')

        fig.run = Button(run_option, 'Run Model', color=button_color)
        fig.run.on_clicked(fig.cnv.run_model)
        fig.run.label.set_fontsize('large')

        fig.merge = Button(merge_option, 'Merge Bodies', color=button_color)
        fig.merge.on_clicked(fig.cnv.merge_bodies)
        fig.merge.label.set_fontsize('large')

        fig.shape_edit = Button(edit_shape_option, 'Shape Change', color=button_color)
        fig.shape_edit.on_clicked(fig.cnv.edit_shape)
        fig.shape_edit.label.set_fontsize('large')

        fig.delete = Button(delete_option, 'Delete Body', color=button_color)
        fig.delete.on_clicked(fig.cnv.delete_body)
        fig.delete.label.set_fontsize('large')
        
        fig.split = Button(split_option, 'Split Body', color=button_color)
        fig.split.on_clicked(fig.cnv.split_body)
        fig.split.label.set_fontsize('large')

        fig.load_btn = Button(load_option, 'Load Model', color=button_color)
        fig.load_btn.on_clicked(fig.cnv.read_model)
        fig.load_btn.label.set_fontsize('large')

        fig.results = Button(plot_results, 'Plot Results', color=button_color)
        fig.results.on_clicked(fig.cnv.Plot_results)
        fig.results.label.set_fontsize('large')
        
        fig.quit_btn = Button(quit_option, 'Quit', color=button_color)
        fig.quit_btn.on_clicked(fig.cnv.quit)
        fig.quit_btn.label.set_fontsize('large')
        
        fig.close_btn = Button(close_model, 'Close Model', color=button_color)
        fig.close_btn.on_clicked(fig.cnv.Close_Model)
        fig.close_btn.label.set_fontsize('large')

        plt.connect('motion_notify_event', fig.cnv.set_location)
        plt.show()
