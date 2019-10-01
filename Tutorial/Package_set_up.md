# Package set up

The `BoXFP.py` should be stored in the same folder as the QuShape packages and the path to the directory for each script should be set to:

`sys.path.append(os.path.abspath('package_directory_path'))`

The following line in the QuShape `funcByRef` and `funcGeneral` packages should also be commented out:

`from imports import QtGui, QtCore`
