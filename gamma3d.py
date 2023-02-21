'''
Tested with Python 3.6.3 & pymedphys 0.36.1
'''

import pymedphys, pydicom, sys, os
from pymedphys import gamma
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5 import QtWidgets, QtCore

# Create QApplication instance
app = QtWidgets.QApplication(sys.argv)

def gamma_3d(dpt, dta, path_to_ref, path_to_eval, cutoff=20, interp=10):
    # Load the two DICOM dose files to be compared
    ref = pymedphys.dicom.zyx_and_dose_from_dataset(pydicom.read_file(path_to_ref))
    eval = pymedphys.dicom.zyx_and_dose_from_dataset(pydicom.read_file(path_to_eval))

    # Get the dose grids and coordinate arrays from the DICOM files
    dose_ref_grid, coord_ref = ref[:2]
    dose_eval_grid, coord_eval = eval[:2]

    # Define the gamma analysis parameters
    dose_percent_threshold = dpt
    distance_mm_threshold = dpt
    lower_percent_dose_cutoff = cutoff
    interp_fraction = interp

    # Calculate the gamma index matrix
    gamma_analysis = gamma(
        dose_ref_grid, coord_ref,
        dose_eval_grid, coord_eval,
        dose_percent_threshold,
        distance_mm_threshold,
        lower_percent_dose_cutoff,
        interp_fraction
    )

    # Filter the gamma index matrix to remove noise
    filtered_gamma = gamma_analysis[~np.isnan(gamma_analysis)]

    # Generate return values
    passed = np.sum(filtered_gamma <= 1) / len(filtered_gamma)
    max_gamma = np.max(filtered_gamma)
    mean_gamma = np.mean(filtered_gamma)
    std_gamma = np.std(filtered_gamma)

    return passed, max_gamma, mean_gamma, std_gamma


# Create a class for the gamma map plotting tool
class GammaViewer(QtWidgets.QWidget):
    def __init__(self, gamma_map):
        super().__init__()

        self.gamma_map = gamma_map.reshape(dose_ref_grid.shape)

        # Set up the figure canvas for the gamma map plot
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)

        # Add the gamma map plot to the figure canvas
        self.ax = self.figure.add_subplot(111)
        self.plot_gamma_map()

        # Set up the scrollbar for scrolling through the dataset
        self.scrollbar = QtWidgets.QScrollBar(QtCore.Qt.Horizontal)
        self.scrollbar.setMaximum(gamma.shape[0] - 1)
        self.scrollbar.valueChanged.connect(self.update_slice)

        # Set up the layout for the GUI
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.canvas)
        layout.addWidget(self.scrollbar)
        self.setLayout(layout)

    def plot_gamma_map(self):
        self.ax.imshow(self.gamma_map, cmap='coolwarm', vmin=0, vmax=2)
        self.ax.set_title('Gamma Map')
        self.ax.set_xlabel('X (mm)')
        self.ax.set_ylabel('Y (mm)')

    def update_slice(self, slice_num):
        self.ax.cla()
        self.ax.imshow(self.gamma_map[slice_num], cmap='coolwarm', vmin=0, vmax=2)
        self.ax.set_title(f'Gamma Map (Slice {slice_num+1}/{self.gamma_map.shape[0]})')
        self.ax.set_xlabel('X (mm)')
        self.ax.set_ylabel('Y (mm)')
        self.canvas.draw()


if __name__ == '__main__':
    root_dir = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\Logfiles\DeliveredPlans\1663630\Doses'
    eval_dir = os.path.join(root_dir, 'eval')
    dpt, dta = 1, 0  # percent, mm
    ref = os.path.join(root_dir, 'RD_R1_HNO.dcm')
    evals = [os.path.join(eval_dir, dcm) for dcm in os.listdir(eval_dir) if dcm.__contains__('RD') and dcm.endswith('.dcm')]
    data = ['Fx\tPass\tMax\tMean\tStd\n']
    for fx, eval in enumerate(evals):
        print(f'>> STARTING 3D gamma evaluation {fx+1}/{len(evals)} with {dpt}%/{dta}mm..')
        passed, max, mean, std = gamma_3d(dpt, dta, ref, eval)
        data.append(f'{fx+1}\t{np.round(passed, 5)}\t{np.round(max, 5)}\t{np.round(mean, 5)}\t{np.round(std, 5)}\n')
    
    with open(os.path.join(root_dir, f'gamma3D_{dpt}p_{dta}mm.txt'), 'w+') as file:
        file.writelines(data)
        file.close()
