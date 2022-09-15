'''Encoding: UTF-8'''
__author__ = 'Lukas C. Wolter, OncoRay ZIK, Dresden, Germany'
__project__ = 'Logfile-based dose calculation & beam statistics'
__version__ = 1.0

import os
import pydicom
from tkinter import filedialog


def spot_sorter(self, fraction_id, beam_id):
    to_be_sorted = self.patient_record_df.loc[(self.patient_record_df['FRACTION_ID'] == fraction_id) & (self.patient_record_df['BEAM_ID'] == beam_id)]
    gtr_angle = to_be_sorted['GANTRY_ANGLE'].iloc[0]
    n_layers = to_be_sorted['TOTAL_LAYERS'].iloc[0]
    found = False

    print('\nTrying to auto-locate patient plan dicoms..')  # read RT plan dicom via filedialog
    for path, dirnames, filenames in os.walk(os.path.join(self.logfile_dir, '..')):
        for fname in filenames:
            if fname.__contains__('RP') and fname.endswith('.dcm'):
                ds = pydicom.read_file(os.path.join(path, fname))
                for i, beam in enumerate(ds.IonBeamSequence):
                    if (beam.BeamName == beam_id or beam.BeamDescription == beam_id) and float(beam.IonControlPointSequence[0].GantryAngle) == gtr_angle and len(beam.IonControlPointSequence) == n_layers * 2:
                        beam_dcm = beam
                        found = True
                        
    while not found:
        print(f'  /!\ Auto-location failed for beam-ID {beam_id} in fraction-ID {fraction_id}, select plan DICOM..')
        target = filedialog.askopenfilename(initialdir=os.path.join(self.logfile_dir, '..'))
        if target == '':
            print('  /x\ Process cancelled by user')
            return None
            
        ds = pydicom.read(target)
        for beam in ds.IonBeamSequence:
            if (beam.BeamName == beam_id or beam.BeamDescription == beam_id) and float(beam.IonControlPointSequence[0].GantryAngle) == gtr_angle and len(beam.IonControlPointSequence) == n_layers * 2:
                beam_dcm = beam
                found = True
    
    return beam_dcm
