'''
Tested with Python 3.11.1 & pymedphys 0.36.1

'''

import pymedphys, pydicom
from pymedphys import gamma
import numpy as np


def gamma_3d(dpt, dta, path_to_ref, path_to_eval, cutoff=10, interp=10, max=1.1, norm=None):
    '''
    dpt .. dose percentage thresh (e.g. 2%), dta .. distance to agreement (e.g. 2mm), 
    path_to_ref/eval .. file paths to reference and evaluation RTDOSE file
    cutoff .. percentage underneath which no gamma will be evaluated
    interp .. fraction of dta used for distance interpolation
    max .. maximum allowed for gamma value (calculation will stop for a voxel if this mark is passed)
    norm .. Normalization point (set to maximum of reference grid by default), can be set to e.g. prescription
    '''

    # Load the two DICOM dose files to be compared
    try:
        ref = pymedphys.dicom.zyx_and_dose_from_dataset(pydicom.read_file(path_to_ref))
        eval = pymedphys.dicom.zyx_and_dose_from_dataset(pydicom.read_file(path_to_eval))
    except:  # if datasets are already loaded by pydicom
        ref = pymedphys.dicom.zyx_and_dose_from_dataset(path_to_ref)
        eval = pymedphys.dicom.zyx_and_dose_from_dataset(path_to_eval)

    # Get the dose grids and coordinate arrays from the DICOM files
    dose_ref_grid, coord_ref = ref[:2]
    dose_eval_grid, coord_eval = eval[:2]

    # Define the gamma analysis parameters
    dose_percent_threshold = dpt
    distance_mm_threshold = dta
    lower_percent_dose_cutoff = cutoff
    interp_fraction = interp
    maximum_gamma = max

    # Calculate the gamma index matrix
    gamma_analysis = gamma(
        dose_ref_grid, coord_ref,
        dose_eval_grid, coord_eval,
        dose_percent_threshold,
        distance_mm_threshold,
        lower_percent_dose_cutoff,
        interp_fraction,
        maximum_gamma,
        global_normalisation=norm,
        quiet=False
    )

    # Filter the gamma index matrix to remove NaNs
    filtered_gamma = gamma_analysis[~np.isnan(gamma_analysis)]

    # Generate return values
    passed = np.sum(filtered_gamma <= 1) / len(filtered_gamma)
    max_gamma = np.max(filtered_gamma)
    mean_gamma = np.mean(filtered_gamma)
    std_gamma = np.std(filtered_gamma)

    return passed, max_gamma, mean_gamma, std_gamma

if __name__=='__main__':
    # Example script use with 3%/3mm and global maximum normalization
    ref = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\Logfiles\DeliveredPlans\280735\Doses\RD1.2.752.243.1.1.20230213142505978.3800.17030.dcm' 
    eval = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\Logfiles\DeliveredPlans\280735\Doses\eval\RD1.2.752.243.1.1.20230303121848345.8620.84151.dcm'
    gammapass, _, _, _ = gamma_3d(3, 3, ref, eval)  # 4 return values, mostly pass rate is of interest
    print(gammapass)
