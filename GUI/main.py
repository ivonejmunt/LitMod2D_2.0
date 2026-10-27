# -*- coding: utf-8 -*-
#####################################################################################################
# This is the main python program of the gui - Updated for Python 3.12 and PyQt5
#####################################################################################################
import sys
import os
import webbrowser
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QPushButton, QLabel, QStackedWidget, QFrame)
from PyQt5.QtGui import QFont
from PyQt5.QtCore import Qt

try:
    from traitsui.menu import OKButton, CancelButton, ApplyButton, UndoButton
except ImportError:
    print("Warning: traitsui not found. Ensure it's installed for full functionality.")

import my_lib

litmod_path = os.getcwd()
LARGE_FONT = QFont("Times", 12)

class StartPage(QWidget):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        
        layout = QVBoxLayout(self)
        layout.setContentsMargins(20, 20, 20, 20)
        layout.setSpacing(10)
        
        label = QLabel("Welcome to LitMod")
        label.setFont(LARGE_FONT)
        label.setAlignment(Qt.AlignCenter) 
        layout.addWidget(label)
        
        btn_build = QPushButton("Build Model")
        btn_build.clicked.connect(lambda: my_lib.build_model(litmod_path))
        layout.addWidget(btn_build)
        
        btn_load = QPushButton("Load Model")
        btn_load.clicked.connect(lambda: my_lib.load_model(litmod_path))
        layout.addWidget(btn_load)
        
        btn_help = QPushButton("Help!")
        btn_help.clicked.connect(lambda: self.controller.show_frame(HelpPage))
        layout.addWidget(btn_help)
        
        tips_text = (
            "Tips!\n"
            "1. Build\n"
            "* Keep the default format for body properties.\n"
            "* Prepare a profile sketch before drawing.\n"
            "* Hover cursor over fields for help info.\n\n"
            "2. Load Model\n"
            "* Requires 'litmod.inp' and 'bodies.out'.\n\n"
            "3. Input/Output\n"
            "* Keep input files and thermodynamic tables in one folder.\n"
            "* Output filenames are preset for plotting (e.g., GMT)."
        )
        
        tips_label = QLabel(tips_text)
        tips_label.setWordWrap(True)
        tips_label.setFrameStyle(QFrame.Panel | QFrame.Sunken)
        tips_label.setMargin(10)
        layout.addWidget(tips_label)
        layout.addStretch()

class HelpPage(QWidget):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller
        
        layout = QVBoxLayout(self)
        layout.setContentsMargins(20, 20, 20, 20)
        
        label = QLabel("Help")
        label.setFont(LARGE_FONT)
        label.setAlignment(Qt.AlignCenter)
        layout.addWidget(label)
        
        btn_web = QPushButton("LitMod webpage")
        btn_web.clicked.connect(self.controller.web_link)
        layout.addWidget(btn_web)
        
        btn_back = QPushButton("Back to home")
        btn_back.clicked.connect(lambda: self.controller.show_frame(StartPage))
        layout.addWidget(btn_back)
        layout.addStretch()

class LithMod(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("LitMod")
        self.setMinimumSize(400, 300)
        
        self.central_widget = QStackedWidget()
        self.setCentralWidget(self.central_widget)
        
        self.frames = {} 
        for F in (StartPage, HelpPage):
            frame = F(self.central_widget, self)
            self.frames[F] = frame
            self.central_widget.addWidget(frame)
        
        self.show_frame(StartPage)

    def show_frame(self, cont):
        frame = self.frames[cont]
        self.central_widget.setCurrentWidget(frame)

    def web_link(self):
        webbrowser.open_new("https://sites.google.com/site/ictjagdl/people")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = LithMod()
    window.show()
    sys.exit(app.exec_()) 