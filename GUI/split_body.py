from matplotlib import pyplot as plt
from matplotlib.widgets import Button, Cursor
from matplotlib.patches import Polygon

class Split_axis:
    def __init__(self, fig, ax, done_option):
        self.ax = ax
        self.fig = fig
        self.done_option = done_option
        self.canvas = self.fig.canvas
        self.tx = [] 
        self.ty = [] 
        self.canvas.mpl_connect('button_press_event', self.select_point)
        # self.mouse_button is not strictly needed in modern matplotlib, keeping logic clean
        
    def select_point(self, event):
        if event.inaxes != self.ax: return
        if event.button == 1: # Left mouse click
            self.tx.append(event.xdata)
            self.ty.append(event.ydata)
            print(f"Point added: {event.xdata}, {event.ydata}")
            
            # Draw the line segment. Using 'r-o' (red line with dots) for better visibility
            self.ax.plot(self.tx, self.ty, 'r-o')
            
            # CRITICAL FIX: Use canvas.draw() instead of plt.show() inside event loops. 
            # plt.show() blocks the GUI and causes severe lag/freezes in PyQt environments.
            self.canvas.draw()

    def done(self, event):
        if event.inaxes == self.done_option:
            print("Saving split points and quitting...")
            # Safely manage the file stream using 'with open'
            with open("split_body_dump.dat", "w") as f:
                for i in range(len(self.tx)):
                    f.write("%f \t %f \n" % (self.tx[i], self.ty[i]))
            
            # Safely close the figure window without killing the main application
            plt.close(self.fig)

def insert_body(l_split, b_split, X, Y, major_xticks, major_yticks):
    # Using plt.figure() without hardcoded number to avoid window conflicts
    fig_split = plt.figure()
    
    # Matplotlib 3.x compatibility: set_window_title moved to manager
    fig_split.canvas.manager.set_window_title('LitMod - Split Body')
    ax_split = fig_split.add_subplot(111)
    
    # CRITICAL FIX: Bind cursor to the figure to prevent Python garbage collection 
    # from destroying it and causing the crosshair to disappear.
    fig_split.cursor = Cursor(ax_split, useblit=True, color='white', linewidth=0.5)
    
    plt.subplots_adjust(bottom=0.2)
    done_option = plt.axes([0.1, 0.05, 0.1, 0.075])
    
    # Fix axis limits based on given parameters
    ax_split.set_xlim((0, X))
    ax_split.set_ylim((-Y, 0))
    ax_split.set_yticks(major_yticks)
    ax_split.set_xticks(major_xticks)  
    ax_split.grid(True)
    
    # Python 3 compatible zip
    # Added facecolor and edgecolor to make the polygon clearly visible during splitting
    poly = Polygon(list(zip(l_split, b_split)), lw=1, facecolor='lightblue', edgecolor='blue', alpha=0.5)
    ax_split.add_patch(poly)
    
    p = Split_axis(fig_split, ax_split, done_option)
    
    # CRITICAL FIX: Bind the Button and Split_axis instance to the figure to prevent 
    # garbage collection, which would otherwise make the 'Done' button unresponsive.
    fig_split.p = p
    fig_split.done_btn = Button(done_option, 'Done', color="lightgrey")
    fig_split.done_btn.on_clicked(p.done)
    
    print("Ready to draw the splitting axis. Click on the plot.")
    plt.show()
    
    # CRITICAL FIX: Removed exit() here. 
    # Calling exit() will terminate the entire main PyQt GUI application!

"""
# Test block kept for standalone testing, translating comments to English
if __name__ == "__main__":
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cursor = Cursor(ax, useblit=True, color='white', linewidth=0.5)
    plt.subplots_adjust(bottom=0.2)
    done_option = plt.axes([0.1, 0.05, 0.1, 0.075])
    ax.set_title('Click to build line segments')
    ax.set_xlim((-1, 2))
    ax.set_ylim((0, 5))

    linebuilder = Split_axis(fig, ax, done_option)
    done = Button(done_option, 'Done', color="lightgrey")
    done.on_clicked(linebuilder.done)
    plt.show()
    
    if plt.get_fignums():
        print("Figure is still open")
    else:
        print("Figure is closed")
        print(linebuilder.tx)
"""