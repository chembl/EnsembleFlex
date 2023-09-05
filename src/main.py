
import sys
import os
#import ntpath
# IMPORTING ALL THE NECESSARY PySide6 MODULES FOR THE APPLICATION.
from PySide6 import QtCore, QtGui, QtWidgets
from PySide6.QtCore import (QCoreApplication, QPropertyAnimation, QDate, QDateTime, QMetaObject, QObject, QPoint,
                            QRect, QSize, QTime, QUrl, Qt, QEvent)
from PySide6.QtGui import (QBrush, QColor, QConicalGradient, QCursor, QFont, QFontDatabase, QIcon, QKeySequence,
                           QLinearGradient, QPalette, QPainter, QPixmap, QRadialGradient)
from PySide6.QtWidgets import *
#import subprocess
#subprocess.check_call(['Rscript', 'myscript.R'], shell=False)

# def path_leaf(path):
#     head, tail = ntpath.split(path)
#     return tail or ntpath.basename(head)


# APPLICATION MAIN WINDOW :
# -----> MAIN APPLICATION CLASS
class MainWindow(QMainWindow):

    def __init__(self):
        """MainWindow constructor."""
        super(MainWindow, self).__init__()

        # Main UI code goes here

        # ----> SET WINDOW TITLE AND ICON
        applicationName = "EnsembleFlex - Flexibility Analysis of Structure Ensembles"
        self.setWindowTitle(applicationName)  # SETS THE APPLICATION NAME IN THE WINDOW TOPBAR
        #############################################################

        # Fonts
        self.font1 = QFont("Sans", 24, QFont.Bold)
        self.font2 = QFont("Sans", 18)
        self.font3 = QFont("Sans", 16)
        self.font4 = QFont("Sans", 14)
        self.font5 = QFont("Sans", 12)

        self.initUI()

    # def callProgram(self):
    #     # run the process # `start` takes the exec and a list of arguments
    #     self.process.start('python',['-u','Robot.py'])

    #############################################################
    # Global Functions
    #############################################################

    def dataReady(self):
        cursor = self.logviewer.textCursor()
        cursor.movePosition(cursor.End)
        cursor.insertText(str(self.process.readAll().data().decode('utf-8', 'ignore').strip()))
        self.logviewer.ensureCursorVisible()

    #############################################################
    # ENSEMBLE NAME
    #############################################################
    def confirmName(self):
        self.name = self.input_name_edit.text()
        print('Ensemble Name: ' + self.name)

    def createGroupName(self):
        groupbox = QGroupBox("")
        groupbox.setStyleSheet("background-color:rgb(47,79,79);")  # darkslategrey

        input_name_lab = QLabel()
        input_name_lab.setText("Ensemble Name:")
        input_name_lab.setFont(self.font2)

        self.input_name_edit = QLineEdit(self)
        self.input_name_edit.setMaxLength(8)
        self.input_name_edit.setStyleSheet("background-color:rgb(190, 230, 230);")  # light cyan
        self.input_name_edit.setAlignment(Qt.AlignRight)
        self.input_name_edit.setFont(self.font2)
        self.input_name_edit.editingFinished.connect(self.confirmName)

        self.confirm_name_btn = QPushButton('Confirm')
        self.confirm_name_btn.setEnabled(True)
        self.confirm_name_btn.clicked.connect(self.confirmName)

        layout = QGridLayout()
        # # addWidget(*Widget, row, column, rowspan, colspan)
        layout.addWidget(input_name_lab, 0, 0)
        layout.addWidget(self.input_name_edit, 0, 1)
        layout.addWidget(self.confirm_name_btn, 0, 2)
        groupbox.setLayout(layout)
        return groupbox

    #############################################################
    # INPUT/OUTPUT
    #############################################################

    def launchSelectionDialog(self):
        option = self.options.index(self.combo.currentText())
        if option == 0:
            self.inputdir = self.getDirectory()
            self.structureFiles = [os.path.join(self.inputdir, f) for f in os.listdir(self.inputdir) if
                                   f.endswith(".pdb") or f.endswith(".ent") or f.endswith(".cif") or f.endswith(".cif.txt")]
            self.selection_response_lab.setText(
                str(len(self.structureFiles)) + ' files selected from folder: ' + str(self.inputdir))
        # elif option == 1:
        #     self.structureFiles = self.getFileNames()
        #     self.selection_response_lab.setText(str(len(self.structureFiles)) + ' selected files: ' + str(
        #         [path_leaf(path) for path in self.structureFiles]))
        # elif option == 2:
        #     response = self.getSaveFileName()
        else:
            print('Got Nothing')

    def launchSelectionDialog2(self):
        self.outputdir = self.getDirectory()
        # Checking if the directory exists and is empty or not
        if os.path.exists(self.outputdir) and not os.listdir(self.outputdir):
            self.selection_response_lab2.setText('Output directory set to: '+str(self.outputdir))
        elif os.path.exists(self.outputdir) and os.listdir(self.outputdir):
            self.selection_response_lab2.setText('Output directory set to: '+str(self.outputdir)+
                                                 '\n Attention: The directory may not be empty and existing files may '
                                                 'be overwritten.')

    def getFileNames(self):
        #file_filter = 'Data File (*.xlsx *.csv *.dat);; Excel File (*.xlsx *.xls)'
        #file_filter = 'Data File (*.pdb *.cif *.cif.txt)'
        response = QFileDialog.getOpenFileNames(
            parent=self,
            caption='Select a data file',
            dir=os.getcwd(),
            filter='Data File (*.pdb *.ent *.cif *.cif.txt)'
        )
        return response[0]

    def getDirectory(self):
        response = QFileDialog.getExistingDirectory(
            parent=self,
            caption='Select a folder'
        )
        return response


    def createGroupInputOutput(self):
        groupbox = QGroupBox("1")
        groupbox.setStyleSheet("background-color:rgb(47,79,79);")  # darkslategrey

        # Input
        selection_label = QLabel()
        selection_label.setMinimumSize(QSize(0, 0))
        selection_label.setMaximumSize(QSize(16777215, 55))
        selection_label.setFont(self.font1)
        # selection_label.setStyleSheet("background-color:rgb(255,255,255);")
        # selection_label.setText(QCoreApplication.translate(u"MainWindow", u"<html><head/><body><p><span style=\" color:#ffffff;\">Select Input Structures</span></p></body></html>", None))
        selection_label.setText("Select Input Structures")

        # self.options = ('Select File Names', 'Select Folder Dir', 'Save File Name')
        # self.options = ('Select whole folder', 'Select distinct files')
        self.options = ('Select whole folder',)

        self.combo = QComboBox()
        self.combo.addItems(self.options)

        selection_btn = QPushButton('Browse')
        selection_btn.clicked.connect(self.launchSelectionDialog)

        self.selection_response_lab = QTextEdit()
        self.selection_response_lab.setReadOnly(True)
        self.selection_response_lab.setFont(self.font4)
        self.selection_response_lab.setMaximumHeight(45)
        self.selection_response_lab.setStyleSheet("background-color:rgb(118, 168, 168);")

        # Output
        selection_label2 = QLabel()
        selection_label2.setMinimumSize(QSize(0, 0))
        selection_label2.setMaximumSize(QSize(16777215, 55))
        selection_label2.setFont(self.font2)
        selection_label2.setText("Create/select Output folder")

        self.selection_btn2 = QPushButton('Browse')
        self.selection_btn2.clicked.connect(self.launchSelectionDialog2)

        self.selection_response_lab2 = QTextEdit()
        self.selection_response_lab2.setReadOnly(True)
        self.selection_response_lab2.setFont(self.font4)
        self.selection_response_lab2.setMaximumHeight(45)
        self.selection_response_lab2.setStyleSheet("background-color:rgb(118, 168, 168);")

        layout = QGridLayout()
        # # addWidget(*Widget, row, column, rowspan, colspan)
        layout.addWidget(selection_label, 0, 0)
        layout.addWidget(self.combo, 1, 0)
        layout.addWidget(selection_btn, 1, 1)
        layout.addWidget(self.selection_response_lab, 2, 0, 1, 2)
        layout.addWidget(selection_label2, 3, 0)
        layout.addWidget(self.selection_btn2, 3, 1)
        layout.addWidget(self.selection_response_lab2, 4, 0, 1, 2)
        groupbox.setLayout(layout)
        return groupbox


    #############################################################
    # SUPERIMPOSE (global)
    #############################################################

    # def btngroup(self, btn):
    #     print(btn.text() + " is selected")

    def hide_confirm_btn(self, state):
        if state == QtCore.Qt.Checked.value:
            self.confirm_btn.setEnabled(True)
        else:
            self.confirm_btn.setEnabled(False)

    def confirmSuper(self):
        self.superimposed = self.inputdir
        self.superimpose_response.setText("Superimposed structures will be taken from " + self.superimposed)
        self.runBio3D_btn.setEnabled(True)
        self.runProDy_btn.setEnabled(True)

    def hide_superimpose_btn1(self, state):
        # if checkbutton is checked and input files and output folder are given
        if state == QtCore.Qt.Checked.value and self.structureFiles and len(self.structureFiles) > 1 and \
                os.path.exists(self.outputdir):
            self.superimpose_btn1.setEnabled(True)
        else:
            self.superimpose_btn1.setEnabled(False)

    def hide_superimpose_btn2(self, state):
        # if checkbutton is checked and input files and output folder are given
        if state == QtCore.Qt.Checked.value and self.structureFiles and len(self.structureFiles) > 1 and \
                os.path.exists(self.outputdir):
            self.superimpose_btn2.setEnabled(True)
        else:
            self.superimpose_btn2.setEnabled(False)

    def callSuperimposeBio3D(self):
        # run the process # `start` takes the exec and a list of arguments
        self.process.start('/usr/local/bin/Rscript',
                           ['./superimpose_bio3d.R', '-d', str(self.inputdir), '-o', str(self.outputdir)])
        self.superimposed = self.outputdir + '/superimposed'
        self.superimpose_response.setText("Superimposed structures are saved in " + self.superimposed)
        self.runBio3D_btn.setEnabled(True)
        self.runProDy_btn.setEnabled(True)

    def callSuperimposeProDy(self):
        # run the process # `start` takes the exec and a list of arguments
        self.process.start('python3',
                           ['./superimpose_prody.py', '-i', str(self.inputdir), '-o', str(self.outputdir)])
        self.superimposed = self.outputdir + '/superimposed/split_ensemble'
        self.superimpose_response.setText("Superimposed structures are saved separately in " + self.superimposed)
        self.runBio3D_btn.setEnabled(True)
        self.runProDy_btn.setEnabled(True)


    def createGroupSuperimp(self):
        groupbox = QGroupBox("2")
        groupbox.setStyleSheet("background-color:rgb(65, 122, 122);")  # darker cyan

        superimpose_label = QLabel()
        superimpose_label.setMinimumSize(QSize(0, 40))
        superimpose_label.setMaximumSize(QSize(16777215, 55))
        superimpose_label.setFont(self.font1)
        superimpose_label.setText("Superimpose Structures (global)")

        self.checkbox1 = QCheckBox("My input structures are already superimposed.")
        self.checkbox1.setChecked(True)
        #self.checkbox1.stateChanged.connect(lambda: self.btnstate(self.checkbox1))
        self.checkbox2 = QCheckBox("Superimpose structures using Bio3D.")
        self.checkbox2.setChecked(False)
        #self.checkbox2.stateChanged.connect(lambda: self.btnstate(self.checkbox2))
        self.checkbox3 = QCheckBox("Superimpose structures using ProDy.")
        self.checkbox3.setChecked(False)
        self.bg = QButtonGroup()
        self.bg.addButton(self.checkbox1, 1)
        self.bg.addButton(self.checkbox2, 2)
        self.bg.addButton(self.checkbox3, 3)

        self.confirm_btn = QPushButton('Confirm')
        self.confirm_btn.setEnabled(True)
        self.checkbox1.stateChanged.connect(self.hide_confirm_btn)
        self.confirm_btn.clicked.connect(self.confirmSuper)

        self.superimpose_btn1 = QPushButton('Run')
        # Disable the button as default - enabled when more than one file is selected in launchSelectionDialog
        self.superimpose_btn1.setEnabled(False)
        self.checkbox2.stateChanged.connect(self.hide_superimpose_btn1)
        self.superimpose_btn1.clicked.connect(self.callSuperimposeBio3D)

        self.superimpose_btn2 = QPushButton('Run')
        # Disable the button as default - enabled when more than one file is selected in launchSelectionDialog
        self.superimpose_btn2.setEnabled(False)
        self.checkbox3.stateChanged.connect(self.hide_superimpose_btn2)
        self.superimpose_btn2.clicked.connect(self.callSuperimposeProDy)

        self.logviewer = QTextEdit()
        self.logviewer.setReadOnly(True)
        self.logviewer.setMinimumHeight(100)
        self.logviewer.setStyleSheet("background-color:rgb(115, 148, 148);")  # grey-cyan

        # QProcess object for external app
        self.process = QtCore.QProcess(self)
        # QProcess emits `readyRead` when there is data to be read
        self.process.readyRead.connect(self.dataReady)
        # Just to prevent accidentally running multiple times
        # Disable the button when process starts, and enable it when it finishes
        self.process.started.connect(lambda: self.superimpose_btn1.setEnabled(False))
        self.process.finished.connect(lambda: self.superimpose_btn1.setEnabled(True))
        self.process.started.connect(lambda: self.superimpose_btn2.setEnabled(False))
        self.process.finished.connect(lambda: self.superimpose_btn2.setEnabled(True))

        self.superimpose_response = QTextEdit()
        self.superimpose_response.setReadOnly(True)
        self.superimpose_response.setFont(self.font4)
        self.superimpose_response.setMaximumHeight(45)
        self.superimpose_response.setStyleSheet("background-color:rgb(118, 168, 168);")

        layout = QGridLayout()
        layout.addWidget(superimpose_label, 0, 0)
        layout.addWidget(self.checkbox1, 1, 0)
        layout.addWidget(self.checkbox2, 2, 0)
        layout.addWidget(self.checkbox3, 3, 0)
        layout.addWidget(self.confirm_btn, 1, 1)
        layout.addWidget(self.superimpose_btn1, 2, 1)
        layout.addWidget(self.superimpose_btn2, 3, 1)
        layout.addWidget(self.superimpose_response, 4, 0, 1, 2)
        # layout.addWidget(self.logviewer, 5, 0, 1, 2)
        groupbox.setLayout(layout)
        return groupbox


    #############################################################
    # ANALYSE (global)
    #############################################################

    def callAnalysisBio3D(self):
        self.outputdirBio3D = self.outputdir + '/Analysis_Bio3D'
        # run the process # `start` takes the exec and a list of arguments
        self.process.start('/usr/local/bin/Rscript',
                           ['./analysis_bio3d.R', '-d', str(self.superimposed), '-o', str(self.outputdirBio3D)])
        self.Bio3d_response.setText("Results will be saved in " + self.outputdirBio3D)

    def callAnalysisProDy(self):
        self.outputdirProDy = self.outputdir + '/Analysis_ProDy'
        # run the process # `start` takes the exec and a list of arguments
        self.process.start('python3',
                           ['./analysis_prody.py', '-i', str(self.superimposed), '-o', str(self.outputdirProDy)])
        self.ProDy_response.setText("Results will be saved in " + self.outputdirProDy)


    def createGroupAnalyse(self):
        groupbox = QGroupBox("3")
        groupbox.setStyleSheet("background-color:rgb(65, 122, 122);")  # darker cyan
        lab = QLabel()
        lab.setMinimumSize(QSize(0, 40))
        lab.setMaximumSize(QSize(16777215, 55))
        lab.setFont(self.font1)
        lab.setText("Flexibility Analysis (global)")

        layout = QGridLayout()
        layout.addWidget(lab, 0, 0)
        layout.addWidget(self.createGroupAnalyseBio3D(), 1, 0)
        layout.addWidget(self.createGroupAnalyseProDy(), 2, 0)
        groupbox.setLayout(layout)
        return groupbox

    def createGroupAnalyseBio3D(self):
        groupbox = QGroupBox("")
        #groupbox.setStyleSheet("background-color:rgb(112, 189, 189);") # cyan
        lab = QLabel()
        lab.setFont(self.font2)
        # lab.setStyleSheet("background-color:rgb(255,255,255);")
        lab.setWordWrap(True)  # making it multi line
        lab.setText("A) using Bio3D:")
        lab2 = QLabel()
        lab2.setFont(self.font4)
        lab2.setWordWrap(True)  # making it multi line
        lab2.setText(" - Root Mean Square Fluctuations (RMSF)\n - Torsion/Dihedral analysis\n"
                     " - Difference distance matrix analysis (DDM)\n - Principal Component Analysis (PCA)\n"
                     " - Ensemble Normal Mode Analysis (eNMA)")
        self.runBio3D_btn = QPushButton('Run')
        # Disable the button as default - enabled when more than one file is selected in launchSelectionDialog
        self.runBio3D_btn.setEnabled(False)
        self.runBio3D_btn.clicked.connect(self.callAnalysisBio3D)

        # QProcess object for external app
        self.process = QtCore.QProcess(self)
        # QProcess emits `readyRead` when there is data to be read
        self.process.readyRead.connect(self.dataReady)
        # Just to prevent accidentally running multiple times
        # Disable the button when process starts, and enable it when it finishes
        self.process.started.connect(lambda: self.runBio3D_btn.setEnabled(False))
        self.process.finished.connect(lambda: self.runBio3D_btn.setEnabled(True))

        self.Bio3d_response = QTextEdit()
        self.Bio3d_response.setReadOnly(True)
        self.Bio3d_response.setFont(self.font4)
        self.Bio3d_response.setMaximumHeight(45)
        self.Bio3d_response.setStyleSheet("background-color:rgb(118, 168, 168);")

        layout = QGridLayout()
        layout.addWidget(lab, 0, 0)
        layout.addWidget(lab2, 1, 0)
        layout.addWidget(self.runBio3D_btn, 2, 0)
        layout.addWidget(self.Bio3d_response, 3, 0, 1, 1)
        groupbox.setLayout(layout)
        return groupbox

    def createGroupAnalyseProDy(self):
        groupbox = QGroupBox("")
        #groupbox.setStyleSheet("background-color:rgb(112, 189, 189);")  # cyan
        lab = QLabel()
        lab.setFont(self.font2)
        # lab.setStyleSheet("background-color:rgb(255,255,255);")
        lab.setText("B) using ProDy:")
        lab2 = QLabel()
        lab2.setFont(self.font4)
        lab2.setWordWrap(True)  # making it multi line
        lab2.setText(" - Root Mean Square Fluctuations (RMSF)\n - Principal Component Analysis (PCA)\n"
                     " - Anisotropic Network Model (ANM) Normal Mode Analysis (NMA)\n"
                     " - Dynamical Domain Decomposition of reference structure")

        self.runProDy_btn = QPushButton('Run')
        # Disable the button as default - enabled when more than one file is selected in launchSelectionDialog
        self.runProDy_btn.setEnabled(False)
        self.runProDy_btn.clicked.connect(self.callAnalysisProDy)

        # QProcess object for external app
        self.process = QtCore.QProcess(self)
        # QProcess emits `readyRead` when there is data to be read
        self.process.readyRead.connect(self.dataReady)
        # Just to prevent accidentally running multiple times
        # Disable the button when process starts, and enable it when it finishes
        self.process.started.connect(lambda: self.runProDy_btn.setEnabled(False))
        self.process.finished.connect(lambda: self.runProDy_btn.setEnabled(True))

        self.ProDy_response = QTextEdit()
        self.ProDy_response.setReadOnly(True)
        self.ProDy_response.setFont(self.font4)
        self.ProDy_response.setMaximumHeight(45)
        self.ProDy_response.setStyleSheet("background-color:rgb(118, 168, 168);")

        layout = QGridLayout()
        layout.addWidget(lab, 0, 0)
        layout.addWidget(lab2, 1, 0)
        layout.addWidget(self.runProDy_btn, 2, 0)
        layout.addWidget(self.ProDy_response, 3, 0, 1, 1)
        groupbox.setLayout(layout)
        return groupbox


    #############################################################
    # COMPARE (global)
    #############################################################

    def callPDBFlex(self):
        self.outputdirPDBFlex = self.outputdir + '/Compare_PDBFlex'
        # run the process # `start` takes the exec and a list of arguments
        self.process.start('python3',
                           ['./call_pdbflex.py', '-i', str(self.superimposed), '-o', str(self.outputdirPDBFlex)])
        self.PDBFlex_response.setText("Results will be saved in " + self.outputdirPDBFlex)

    def callAlphaFold(self):
        self.outputdirAlphaFold = self.outputdir + '/Compare_AlphaFold'
        # run the process # `start` takes the exec and a list of arguments
        self.process.start('python3',
                           ['./call_.py', '-i', str(self.superimposed), '-o', str(self.outputdirAlphaFold)])
        self.AlphaFold_response.setText("Results will be saved in " + self.outputdirAlphaFold)


    def createGroupCompare(self):
        groupbox = QGroupBox("4")
        groupbox.setStyleSheet("background-color:rgb(65, 122, 122);")  # darker cyan
        lab = QLabel()
        lab.setMinimumSize(QSize(0, 40))
        lab.setMaximumSize(QSize(16777215, 55))
        lab.setFont(self.font1)
        # lab4.setStyleSheet("background-color:rgb(255,255,255);")
        lab.setText("Compare Flexibility (global)")

        layout = QGridLayout()
        layout.addWidget(lab, 0, 0)
        layout.addWidget(self.createGroupComparePDBFlex(), 1, 0)
        layout.addWidget(self.createGroupCompareAlphaFold(), 2, 0)
        groupbox.setLayout(layout)
        return groupbox

    def createGroupComparePDBFlex(self):
        groupbox = QGroupBox("")
        #groupbox.setStyleSheet("background-color:rgb(112, 189, 189);") # cyan
        lab = QLabel()
        lab.setFont(self.font2)
        # lab.setStyleSheet("background-color:rgb(255,255,255);")
        lab.setText("A) Compare with PDBFlex cluster")

        lab2 = QLabel()
        lab2.setFont(self.font4)
        lab2.setWordWrap(True)  # making it multi line
        lab2.setText(" - Download PDBFlex cluster\n - Root Mean Square Fluctuations (RMSF)\n"
                     " - Principal Component Analysis (PCA)\n - Normal Mode Analysis (NMA)")

        self.runPDBFlex_btn = QPushButton('Run')
        # Disable the button as default - enabled when more than one file is selected in launchSelectionDialog
        self.runPDBFlex_btn.setEnabled(False)
        self.runPDBFlex_btn.clicked.connect(self.callPDBFlex)

        self.PDBFlex_response = QTextEdit()
        self.PDBFlex_response.setReadOnly(True)
        self.PDBFlex_response.setFont(self.font4)
        self.PDBFlex_response.setMaximumHeight(45)
        self.PDBFlex_response.setStyleSheet("background-color:rgb(118, 168, 168);")

        layout = QGridLayout()
        layout.addWidget(lab, 0, 0)
        layout.addWidget(lab2, 1, 0)
        layout.addWidget(self.runPDBFlex_btn, 2, 0)
        layout.addWidget(self.PDBFlex_response, 3, 0, 1, 1)
        groupbox.setLayout(layout)
        return groupbox

    def createGroupCompareAlphaFold(self):
        groupbox = QGroupBox("")
        #groupbox.setStyleSheet("background-color:rgb(112, 189, 189);") # cyan
        lab = QLabel()
        lab.setFont(self.font2)
        # lab.setStyleSheet("background-color:rgb(255,255,255);")
        lab.setText("B) Compare with AlphaFold prediction")

        lab2 = QLabel()
        lab2.setFont(self.font4)
        lab2.setWordWrap(True)  # making it multi line
        lab2.setText(" - Download AlphaFold prediction\n - Root Mean Square Fluctuations (RMSF)\n"
                     " - Principal Component Analysis (PCA)\n - Normal Mode Analysis (NMA)")

        self.runAlphaFold_btn = QPushButton('Run')
        # Disable the button as default - enabled when more than one file is selected in launchSelectionDialog
        self.runAlphaFold_btn.setEnabled(False)
        self.runAlphaFold_btn.clicked.connect(self.callAlphaFold)

        self.AlphaFold_response = QTextEdit()
        self.AlphaFold_response.setReadOnly(True)
        self.AlphaFold_response.setFont(self.font4)
        self.AlphaFold_response.setMaximumHeight(45)
        self.AlphaFold_response.setStyleSheet("background-color:rgb(118, 168, 168);")

        layout = QGridLayout()
        layout.addWidget(lab, 0, 0)
        layout.addWidget(lab2, 1, 0)
        layout.addWidget(self.runAlphaFold_btn, 2, 0)
        layout.addWidget(self.AlphaFold_response, 3, 0, 1, 1)
        groupbox.setLayout(layout)
        return groupbox


    #############################################################
    # Binding Site
    #############################################################

    def createGroupBindingSite(self):
        groupbox = QGroupBox("5")
        groupbox.setStyleSheet("background-color:rgb(86, 163, 163);")  # cyan
        lab = QLabel()
        lab.setMinimumSize(QSize(0, 40))
        lab.setMaximumSize(QSize(16777215, 55))
        lab.setFont(self.font1)
        # lab.setStyleSheet("background-color:rgb(255,255,255);")
        lab.setText("Focus on Binding Site")
        layout = QGridLayout()
        layout.addWidget(lab, 0, 0)
        layout.addWidget(self.createGroupSuperimposeBindingSite(), 1, 0)
        layout.addWidget(self.createGroupAnalyseBindingSite(), 0, 1, 2, 1)
        groupbox.setLayout(layout)
        return groupbox

    # SUPERIMPOSE (Binding Site)
    #############################################################
    def hide_confirm_bs_btn(self, state):
        if state == QtCore.Qt.Checked.value:
            self.confirm_bs_btn.setEnabled(True)
        else:
            self.confirm_bs_btn.setEnabled(False)

    def confirmSuper_bs(self):
        self.superimposed_bs = self.superimposed
        self.superimpose_bs_response.setText("Superimposed structures will be taken from " + self.superimposed_bs)
        self.runBio3D_bs_btn.setEnabled(True)
        self.runProDy_bs_btn.setEnabled(True)

    def hide_superimpose_bs_btn1(self, state):
        # if checkbutton is checked and input files and output folder are given
        if state == QtCore.Qt.Checked.value and self.structureFiles and len(self.structureFiles) > 1 and \
                os.path.exists(self.outputdir):
            self.superimpose_bs_btn1.setEnabled(True)
        else:
            self.superimpose_bs_btn1.setEnabled(False)

    def hide_superimpose_bs_btn2(self, state):
        # if checkbutton is checked and input files and output folder are given
        if state == QtCore.Qt.Checked.value and self.structureFiles and len(self.structureFiles) > 1 and \
                os.path.exists(self.outputdir):
            self.superimpose_bs_btn2.setEnabled(True)
        else:
            self.superimpose_bs_btn2.setEnabled(False)

    def callSuperimposeBio3D_bs(self):
        # run the process # `start` takes the exec and a list of arguments
        self.process.start('/usr/local/bin/Rscript',
                           ['./superimpose_bio3d_bs.R', '-d', str(self.inputdir), '-o', str(self.outputdir)])
        self.superimposed_bs = self.outputdir + '/superimposed_bs'
        self.superimpose_bs_response.setText("Superimposed structures are saved in " + self.superimposed_bs)
        self.runBio3D_bs_btn.setEnabled(True)
        self.runProDy_bs_btn.setEnabled(True)

    def callSuperimposeProDy_bs(self):
        # run the process # `start` takes the exec and a list of arguments
        self.process.start('python3',
                           ['./superimpose_prody_bs.py', '-i', str(self.inputdir), '-o', str(self.outputdir)])
        self.superimposed_bs = self.outputdir + '/superimposed_bs'
        self.superimpose_bs_response.setText("Superimposed structures are saved in " + self.superimposed_bs)
        self.runBio3D_bs_btn.setEnabled(True)
        self.runProDy_bs_btn.setEnabled(True)


    def createGroupSuperimposeBindingSite(self):
        groupbox = QGroupBox("")
        #groupbox.setStyleSheet("background-color:rgb(112, 189, 189);") # cyan
        lab = QLabel()
        lab.setMinimumSize(QSize(0, 40))
        lab.setMaximumSize(QSize(16777215, 55))
        lab.setFont(self.font2)
        # lab.setStyleSheet("background-color:rgb(255,255,255);")
        lab.setText("Superimpose Binding Site")

        self.checkbox_bs1 = QCheckBox("Skip and use globally superimposed structures.")
        self.checkbox_bs1.setChecked(True)
        self.checkbox_bs2 = QCheckBox("Superimpose binding site residues using Bio3D.")
        self.checkbox_bs2.setChecked(False)
        self.checkbox_bs3 = QCheckBox("Superimpose binding site residues using ProDy.")
        self.checkbox_bs3.setChecked(False)
        self.bg_bs = QButtonGroup()
        self.bg_bs.addButton(self.checkbox_bs1, 1)
        self.bg_bs.addButton(self.checkbox_bs2, 2)
        self.bg_bs.addButton(self.checkbox_bs3, 3)

        self.confirm_bs_btn = QPushButton('Confirm')
        self.confirm_bs_btn.setEnabled(True)
        self.checkbox_bs1.stateChanged.connect(self.hide_confirm_bs_btn)
        self.confirm_bs_btn.clicked.connect(self.confirmSuper_bs)

        self.superimpose_bs_btn1 = QPushButton('Run')
        # Disable the button as default - enabled when more than one file is selected in launchSelectionDialog
        self.superimpose_bs_btn1.setEnabled(False)
        self.checkbox2.stateChanged.connect(self.hide_superimpose_bs_btn1)
        self.superimpose_bs_btn1.clicked.connect(self.callSuperimposeBio3D_bs)

        self.superimpose_bs_btn2 = QPushButton('Run')
        # Disable the button as default - enabled when more than one file is selected in launchSelectionDialog
        self.superimpose_bs_btn2.setEnabled(False)
        self.checkbox_bs3.stateChanged.connect(self.hide_superimpose_bs_btn2)
        self.superimpose_bs_btn2.clicked.connect(self.callSuperimposeProDy_bs)

        self.superimpose_bs_response = QTextEdit()
        self.superimpose_bs_response.setReadOnly(True)
        self.superimpose_bs_response.setFont(self.font4)
        self.superimpose_bs_response.setMaximumHeight(45)
        self.superimpose_bs_response.setStyleSheet("background-color:rgb(118, 168, 168);")

        layout = QGridLayout()
        layout.addWidget(lab, 0, 0)
        layout.addWidget(self.checkbox_bs1, 1, 0)
        layout.addWidget(self.checkbox_bs2, 2, 0)
        layout.addWidget(self.checkbox_bs3, 3, 0)
        layout.addWidget(self.confirm_bs_btn, 1, 1)
        layout.addWidget(self.superimpose_bs_btn1, 2, 1)
        layout.addWidget(self.superimpose_bs_btn2, 3, 1)
        layout.addWidget(self.superimpose_bs_response, 4, 0, 1, 2)
        groupbox.setLayout(layout)
        return groupbox

    # ANALYSE (Binding Site)
    #############################################################
    def createGroupAnalyseBindingSite(self):
        groupbox = QGroupBox("")
        #groupbox.setStyleSheet("background-color:rgb(112, 189, 189);") # cyan
        lab = QLabel()
        lab.setMinimumSize(QSize(0, 40))
        lab.setMaximumSize(QSize(16777215, 55))
        lab.setFont(self.font2)
        # lab.setStyleSheet("background-color:rgb(255,255,255);")
        lab.setText("Binding Site Analysis")
        layout = QGridLayout()
        layout.addWidget(lab, 0, 0)
        groupbox.setLayout(layout)
        return groupbox


    #############################################################

    def initUI(self):

        # Menubar
        menubar = self.menuBar()  # QMenuBar
        file_menu = menubar.addMenu('File')  # QMenu
        file_menu.addSeparator()
        file_menu.addAction('Quit', self.close)

        scrollArea = QtWidgets.QScrollArea()  # Scroll Area which contains the widgets, set as the centralWidget
        scrollArea.setWidget(self.myCentralWidget())
        scrollArea.setWidgetResizable(True)  # Scroll Area Properties
        #scrollArea.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        #scrollArea.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)

        #self.setCentralWidget(self.myCentralWidget())
        self.setCentralWidget(scrollArea)

        self.show()


    def myCentralWidget(self):

        #self.window_width, self.window_height = 1300, 950
        #self.setMinimumSize(self.window_width, self.window_height)
        self.showMaximized()
        #self.setStyleSheet("background-color:rgb(91,90,90);")
        #self.setStyleSheet("background-color:rgb(0,124,130);")
        #self.setStyleSheet("background-color:rgb(112, 189, 189);") # cyan
        #self.setStyleSheet("background-color:rgb(255,255,255);") # white
        self.setStyleSheet("background-color:rgb(192,192,192);") # silver

        # # creating a vertical box layout
        # layout = QVBoxLayout()
        # creating a grid layout
        layout = QGridLayout()
        # # creating a form layout
        # layout = QFormLayout()

        self.setLayout(layout)

        #############################################################
        layout.addWidget(self.createGroupName(), 0, 0, 1, 2)
        layout.addWidget(self.createGroupInputOutput(), 1, 0, 2, 2)
        layout.addWidget(self.createGroupSuperimp(), 3, 0, 2, 2)
        layout.addWidget(self.logviewer, 5, 0, 2, 2)
        layout.addWidget(self.createGroupAnalyse(), 0, 3, 5, 2)
        layout.addWidget(self.createGroupCompare(), 0, 5, 5, 2)
        layout.addWidget(self.createGroupBindingSite(), 5, 3, 2, 4)


        #############################################################

        w = QWidget()
        w.setLayout(layout)
        return(w)



    # def getSaveFileName(self):
    #     file_filter = 'Data File(*.pdb *.cif *.cif.txt)'
    #     response = QFileDialog.getSaveFileName(
    #         parent=self,
    #         caption='Select a data file',
    #         dir= 'Data File.dat',
    #         filter='Data File (*.pdb *.cif *.cif.txt)'
    #     )
    #     print(response)
    #     return response[0]



        ###############################


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()

    # Open the qss styles file and read in the css-alike styling code
    with open('./styles.qss', 'r') as f:
        # Set the stylesheet of the application
        app.setStyleSheet(f.read())

    #myApp = MyApp()
    #myApp.show()

    sys.exit(app.exec())

#######################################################################################################################
