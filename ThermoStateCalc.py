# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ThermoStateCalc2.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.
'''this was first fixed by seth to more match what was shown for state 1 and 2 and state change
then irvin came along and fixed it to have intenral energy thru quailty, and to make state change say diffrence as well'''

from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(1074, 911)
        self._grp_Units = QtWidgets.QGroupBox(Dialog)
        self._grp_Units.setGeometry(QtCore.QRect(10, 20, 540, 60))
        self._grp_Units.setMinimumSize(QtCore.QSize(506, 0))
        self._grp_Units.setObjectName("_grp_Units")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self._grp_Units)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self._rdo_SI = QtWidgets.QRadioButton(self._grp_Units)
        self._rdo_SI.setChecked(True)
        self._rdo_SI.setObjectName("_rdo_SI")
        self.horizontalLayout_2.addWidget(self._rdo_SI)
        self._rdo_English = QtWidgets.QRadioButton(self._grp_Units)
        self._rdo_English.setObjectName("_rdo_English")
        self.horizontalLayout_2.addWidget(self._rdo_English)
        self._grp_StateProperties = QtWidgets.QGroupBox(Dialog)
        self._grp_StateProperties.setGeometry(QtCore.QRect(10, 80, 1051, 411))
        self._grp_StateProperties.setObjectName("_grp_StateProperties")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self._grp_StateProperties)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.groupBox = QtWidgets.QGroupBox(self._grp_StateProperties)
        self.groupBox.setObjectName("groupBox")
        self._lbl_Property1 = QtWidgets.QLabel(self.groupBox)
        self._lbl_Property1.setGeometry(QtCore.QRect(10, 20, 59, 22))
        self._lbl_Property1.setObjectName("_lbl_Property1")
        self._cmb_Property1_2 = QtWidgets.QComboBox(self.groupBox)
        self._cmb_Property1_2.setGeometry(QtCore.QRect(10, 40, 181, 22))
        self._cmb_Property1_2.setObjectName("_cmb_Property1_2")
        self._cmb_Property1_2.addItem("")
        self._cmb_Property1_2.addItem("")
        self._cmb_Property1_2.addItem("")
        self._cmb_Property1_2.addItem("")
        self._cmb_Property1_2.addItem("")
        self._cmb_Property1_2.addItem("")
        self._cmb_Property1_2.addItem("")
        self._le_Property1 = QtWidgets.QLineEdit(self.groupBox)
        self._le_Property1.setGeometry(QtCore.QRect(10, 70, 181, 31))
        self._le_Property1.setObjectName("_le_Property1")
        self._lbl_Property1_Units = QtWidgets.QLabel(self.groupBox)
        self._lbl_Property1_Units.setGeometry(QtCore.QRect(200, 70, 31, 31))
        self._lbl_Property1_Units.setObjectName("_lbl_Property1_Units")
        self._lbl_Property2 = QtWidgets.QLabel(self.groupBox)
        self._lbl_Property2.setGeometry(QtCore.QRect(240, 20, 71, 16))
        self._lbl_Property2.setObjectName("_lbl_Property2")
        self._cmb_Property2 = QtWidgets.QComboBox(self.groupBox)
        self._cmb_Property2.setGeometry(QtCore.QRect(240, 40, 239, 22))
        self._cmb_Property2.setObjectName("_cmb_Property2")
        self._cmb_Property2.addItem("")
        self._cmb_Property2.addItem("")
        self._cmb_Property2.addItem("")
        self._cmb_Property2.addItem("")
        self._cmb_Property2.addItem("")
        self._cmb_Property2.addItem("")
        self._cmb_Property2.addItem("")
        self._le_Property2 = QtWidgets.QLineEdit(self.groupBox)
        self._le_Property2.setGeometry(QtCore.QRect(240, 70, 239, 31))
        self._le_Property2.setObjectName("_le_Property2")
        self._lbl_Property2_Units = QtWidgets.QLabel(self.groupBox)
        self._lbl_Property2_Units.setGeometry(QtCore.QRect(480, 70, 31, 31))
        self._lbl_Property2_Units.setObjectName("_lbl_Property2_Units")
        self.verticalLayout_2.addWidget(self.groupBox)
        self.groupBox_2 = QtWidgets.QGroupBox(self._grp_StateProperties)
        self.groupBox_2.setObjectName("groupBox_2")
        self.label = QtWidgets.QLabel(self.groupBox_2)
        self.label.setGeometry(QtCore.QRect(10, 20, 101, 16))
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(self.groupBox_2)
        self.label_2.setGeometry(QtCore.QRect(250, 20, 71, 16))
        self.label_2.setObjectName("label_2")
        self._cmb_Property1_3 = QtWidgets.QComboBox(self.groupBox_2)
        self._cmb_Property1_3.setGeometry(QtCore.QRect(10, 50, 181, 22))
        self._cmb_Property1_3.setObjectName("_cmb_Property1_3")
        self._cmb_Property1_3.addItem("")
        self._cmb_Property1_3.addItem("")
        self._cmb_Property1_3.addItem("")
        self._cmb_Property1_3.addItem("")
        self._cmb_Property1_3.addItem("")
        self._cmb_Property1_3.addItem("")
        self._cmb_Property1_3.addItem("")
        self._cmb_Property2_3 = QtWidgets.QComboBox(self.groupBox_2)
        self._cmb_Property2_3.setGeometry(QtCore.QRect(250, 50, 211, 22))
        self._cmb_Property2_3.setObjectName("_cmb_Property2_3")
        self._cmb_Property2_3.addItem("")
        self._cmb_Property2_3.addItem("")
        self._cmb_Property2_3.addItem("")
        self._cmb_Property2_3.addItem("")
        self._cmb_Property2_3.addItem("")
        self._cmb_Property2_3.addItem("")
        self._cmb_Property2_3.addItem("")
        self._le_Property1_2 = QtWidgets.QLineEdit(self.groupBox_2)
        self._le_Property1_2.setGeometry(QtCore.QRect(10, 80, 181, 31))
        self._le_Property1_2.setObjectName("_le_Property1_2")
        self._lbl_Property1_Units_2 = QtWidgets.QLabel(self.groupBox_2)
        self._lbl_Property1_Units_2.setGeometry(QtCore.QRect(200, 80, 31, 31))
        self._lbl_Property1_Units_2.setObjectName("_lbl_Property1_Units_2")
        self._le_Property2_2 = QtWidgets.QLineEdit(self.groupBox_2)
        self._le_Property2_2.setGeometry(QtCore.QRect(250, 80, 211, 31))
        self._le_Property2_2.setObjectName("_le_Property2_2")
        self._lbl_Property2_Units_2 = QtWidgets.QLabel(self.groupBox_2)
        self._lbl_Property2_Units_2.setGeometry(QtCore.QRect(470, 80, 31, 31))
        self._lbl_Property2_Units_2.setObjectName("_lbl_Property2_Units_2")
        self.verticalLayout_2.addWidget(self.groupBox_2)
        self._pb_Calculate = QtWidgets.QPushButton(self._grp_StateProperties)
        self._pb_Calculate.setObjectName("_pb_Calculate")
        self.verticalLayout_2.addWidget(self._pb_Calculate)
        self._grp_StateProperties_2 = QtWidgets.QGroupBox(Dialog)
        self._grp_StateProperties_2.setGeometry(QtCore.QRect(10, 500, 1041, 391))
        self._grp_StateProperties_2.setObjectName("_grp_StateProperties_2")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self._grp_StateProperties_2)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.groupBox_3 = QtWidgets.QGroupBox(self._grp_StateProperties_2)
        self.groupBox_3.setObjectName("groupBox_3")
        self._lbl_StateProperties_2 = QtWidgets.QLabel(self.groupBox_3)
        self._lbl_StateProperties_2.setGeometry(QtCore.QRect(10, 40, 482, 121))
        self._lbl_StateProperties_2.setObjectName("_lbl_StateProperties_2")
        self._lbl_State_2 = QtWidgets.QLabel(self.groupBox_3)
        self._lbl_State_2.setGeometry(QtCore.QRect(10, 20, 482, 16))
        self._lbl_State_2.setObjectName("_lbl_State_2")
        self.horizontalLayout.addWidget(self.groupBox_3)
        self.groupBox_4 = QtWidgets.QGroupBox(self._grp_StateProperties_2)
        self.groupBox_4.setObjectName("groupBox_4")
        self._lbl_State = QtWidgets.QLabel(self.groupBox_4)
        self._lbl_State.setGeometry(QtCore.QRect(10, 20, 482, 16))
        self._lbl_State.setObjectName("_lbl_State")
        self._lbl_StateProperties = QtWidgets.QLabel(self.groupBox_4)
        self._lbl_StateProperties.setGeometry(QtCore.QRect(10, 40, 482, 111))
        self._lbl_StateProperties.setObjectName("_lbl_StateProperties")
        self.horizontalLayout.addWidget(self.groupBox_4)
        self.groupBox_5 = QtWidgets.QGroupBox(self._grp_StateProperties_2)
        self.groupBox_5.setObjectName("groupBox_5")
        self._lbl_State_3 = QtWidgets.QLabel(self.groupBox_5)
        self._lbl_State_3.setGeometry(QtCore.QRect(10, 20, 482, 16))
        self._lbl_State_3.setObjectName("_lbl_State_3")
        self._lbl_StateProperties_3 = QtWidgets.QLabel(self.groupBox_5)
        self._lbl_StateProperties_3.setGeometry(QtCore.QRect(10, 40, 482, 131))
        self._lbl_StateProperties_3.setObjectName("_lbl_StateProperties_3")
        self.horizontalLayout.addWidget(self.groupBox_5)

        self.retranslateUi(Dialog)
        self._cmb_Property2.setCurrentIndex(1)
        self._cmb_Property2_3.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self._grp_Units.setTitle(_translate("Dialog", "System of Units"))
        self._rdo_SI.setText(_translate("Dialog", "SI"))
        self._rdo_English.setText(_translate("Dialog", "English"))
        self._grp_StateProperties.setTitle(_translate("Dialog", "State Properties"))
        self.groupBox.setTitle(_translate("Dialog", "State 1"))
        self._lbl_Property1.setText(_translate("Dialog", "Property 1"))
        self._cmb_Property1_2.setCurrentText(_translate("Dialog", "Pressure (p)"))
        self._cmb_Property1_2.setItemText(0, _translate("Dialog", "Pressure (p)"))
        self._cmb_Property1_2.setItemText(1, _translate("Dialog", "Temperature (T)"))
        self._cmb_Property1_2.setItemText(2, _translate("Dialog", "Quality (x)"))
        self._cmb_Property1_2.setItemText(3, _translate("Dialog", "Specific Internal Energy (u)"))
        self._cmb_Property1_2.setItemText(4, _translate("Dialog", "Specific Enthalpy (h)"))
        self._cmb_Property1_2.setItemText(5, _translate("Dialog", "Specific Volume (v)"))
        self._cmb_Property1_2.setItemText(6, _translate("Dialog", "Specific Entropy (s)"))
        self._le_Property1.setText(_translate("Dialog", "1.0"))
        self._lbl_Property1_Units.setText(_translate("Dialog", "Bar"))
        self._lbl_Property2.setText(_translate("Dialog", "Property 2"))
        self._cmb_Property2.setItemText(0, _translate("Dialog", "Pressure (p)"))
        self._cmb_Property2.setItemText(1, _translate("Dialog", "Temperature (T)"))
        self._cmb_Property2.setItemText(2, _translate("Dialog", "Quality (x)"))
        self._cmb_Property2.setItemText(3, _translate("Dialog", "Specific Internal Energy (u)"))
        self._cmb_Property2.setItemText(4, _translate("Dialog", "Specific Enthalpy (h)"))
        self._cmb_Property2.setItemText(5, _translate("Dialog", "Specific Volume (v)"))
        self._cmb_Property2.setItemText(6, _translate("Dialog", "Specific Entropy (s)"))
        self._le_Property2.setText(_translate("Dialog", "100.0"))
        self._lbl_Property2_Units.setText(_translate("Dialog", "C"))
        self.groupBox_2.setTitle(_translate("Dialog", "State 2"))
        self.label.setText(_translate("Dialog", "Property 1"))
        self.label_2.setText(_translate("Dialog", "Property 2"))
        self._cmb_Property1_3.setCurrentText(_translate("Dialog", "Pressure (p)"))
        self._cmb_Property1_3.setItemText(0, _translate("Dialog", "Pressure (p)"))
        self._cmb_Property1_3.setItemText(1, _translate("Dialog", "Temperature (T)"))
        self._cmb_Property1_3.setItemText(2, _translate("Dialog", "Quality (x)"))
        self._cmb_Property1_3.setItemText(3, _translate("Dialog", "Specific Internal Energy (u)"))
        self._cmb_Property1_3.setItemText(4, _translate("Dialog", "Specific Enthalpy (h)"))
        self._cmb_Property1_3.setItemText(5, _translate("Dialog", "Specific Volume (v)"))
        self._cmb_Property1_3.setItemText(6, _translate("Dialog", "Specific Entropy (s)"))
        self._cmb_Property2_3.setItemText(0, _translate("Dialog", "Pressure (p)"))
        self._cmb_Property2_3.setItemText(1, _translate("Dialog", "Temperature (T)"))
        self._cmb_Property2_3.setItemText(2, _translate("Dialog", "Quality (x)"))
        self._cmb_Property2_3.setItemText(3, _translate("Dialog", "Specific Internal Energy (u)"))
        self._cmb_Property2_3.setItemText(4, _translate("Dialog", "Specific Enthalpy (h)"))
        self._cmb_Property2_3.setItemText(5, _translate("Dialog", "Specific Volume (v)"))
        self._cmb_Property2_3.setItemText(6, _translate("Dialog", "Specific Entropy (s)"))
        self._le_Property1_2.setText(_translate("Dialog", "1.0"))
        self._lbl_Property1_Units_2.setText(_translate("Dialog", "Bar"))
        self._le_Property2_2.setText(_translate("Dialog", "100.0"))
        self._lbl_Property2_Units_2.setText(_translate("Dialog", "C"))
        self._pb_Calculate.setText(_translate("Dialog", "Calculate"))
        self._grp_StateProperties_2.setTitle(_translate("Dialog", "Specified Properties"))
        self.groupBox_3.setTitle(_translate("Dialog", "State 2 Properties"))
        self._lbl_StateProperties_2.setText(_translate("Dialog", "Pressure = 1000 kPa\n"
"Temperature = 100 C [TSat=99.606(C)] \n"
" internal Energy = 2506.171 (W/kg) \n"
" Enthalpy= 2675.767 (W/kg) \n"
" Entropy = 7.361 (W/kg*C) \n"
" Specfic Volume = 1.696 (m^3/kg) \n"
" Quality =1.0"))
        self._lbl_State_2.setText(_translate("Dialog", "Region:  saturated"))
        self.groupBox_4.setTitle(_translate("Dialog", "State 1 Properties"))
        self._lbl_State.setText(_translate("Dialog", "Region:  saturated"))
        self._lbl_StateProperties.setText(_translate("Dialog", "Pressure = 1000 kPa\n"
"Temperature = 100 C [TSat=99.606(C)] \n"
" internal Energy = 2506.171 (W/kg) \n"
" Enthalpy= 2675.767 (W/kg) \n"
" Entropy = 7.361 (W/kg*C) \n"
" Specfic Volume = 1.696 (m^3/kg) \n"
" Quality =1.0"))
        self.groupBox_5.setTitle(_translate("Dialog", "State Changes"))
        self._lbl_State_3.setText(_translate("Dialog", "Region:  saturated to saturated"))
        self._lbl_StateProperties_3.setText(_translate("Dialog", "Pressure diffrence = 0 \n"
"Temperature diffrence = 0 C \n"
" temp sat diffrence = 0 \n"
" internal diff Energy = 0 (W/kg) \n"
" Enthalpy= 2675.767 (W/kg) \n"
" Entropy diffrence = 0 (W/kg*C) \n"
" Specfic Volume diffrence = 0 (m^3/kg) \n"
" Quality diffrence = 1.0"))
