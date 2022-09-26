'''Encoding: UTF-8'''
__author__ = 'Lukas C. Wolter, OncoRay ZIK, Dresden, Germany'
__project__ = 'Logfile-based dose calculation & beam statistics'
__version__ = 1.0

import os
import numpy as np
import pydicom
from tkinter import filedialog
from matplotlib import pyplot as plt


def spot_sorter(self, fraction_id, beam_id):
    to_be_sorted = self.patient_record_df.loc[(self.patient_record_df['FRACTION_ID'] == fraction_id) & (self.patient_record_df['BEAM_ID'] == beam_id)]
    gtr_angle = to_be_sorted['GANTRY_ANGLE'].iloc[0]
    n_layers = to_be_sorted['TOTAL_LAYERS'].iloc[0]
    found = False
    sorted_positions, sorted_metersets = [], []

    # auto-location of plan DICOM
    print('\nTrying to auto-locate patient plan dicoms..')  # read RT plan dicom via filedialog
    for path, dirnames, filenames in os.walk(os.path.join(self.logfile_dir, '..')):
        for fname in filenames:
            if fname.__contains__('RP') and fname.endswith('.dcm'):
                ds = pydicom.read_file(os.path.join(path, fname))
                for i, beam in enumerate(ds.IonBeamSequence):
                    if (beam.BeamName == beam_id or beam.BeamDescription == beam_id) and float(beam.IonControlPointSequence[0].GantryAngle) == gtr_angle and len(beam.IonControlPointSequence) == n_layers * 2:
                        beam_dcm = beam
                        found = True
                        
    while not found:  # open dicom file manually if failed
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
    
    # sorting
    for layer_id in to_be_sorted['LAYER_ID'].drop_duplicates():
        try:
            plan_layer = beam_dcm.IonControlPointSequence[layer_id * 2]
        except IndexError:
            print(f'''  /!\ Layer-ID mismatch, skipping layer #{layer_id + 1} in beam {beam_id}, fraction {fraction_id}''')
            continue

        plan_spotmap = plan_layer.ScanSpotPositionMap
        plan_x, plan_y = [], []
        for i, spot in enumerate(plan_spotmap):
            if i % 2 == 0:
                plan_x.append(spot)
            else:
                plan_y.append(spot)
        plan_xy = [tup for tup in zip(plan_x, plan_y)]
        try:
            plan_mu = [mu for mu in plan_layer.ScanSpotMetersetWeights]
        except TypeError:
            plan_mu = [plan_layer.ScanSpotMetersetWeights]

        log_layer = to_be_sorted.loc[(to_be_sorted['LAYER_ID'] == layer_id) & (to_be_sorted['FRACTION_ID'] == fraction_id)]
        log_xy = [tup for tup in zip(log_layer['X_POSITION(mm)'], log_layer['Y_POSITION(mm)'])]
        log_mu  = log_layer['MU']
        log_xy_sorted = [(np.nan, np.nan) for _ in range(max(len(log_xy), len(plan_xy)))]
        log_mu_sorted = [np.nan for _ in log_xy_sorted]

        # match (x,y)-positions to plan, transform MU list equally
        for i, log_spot in enumerate(log_xy):
            shifts = [np.array(plan_spot) - np.array(log_spot) for plan_spot in plan_xy]
            dists = [np.abs((shift).dot(shift)) for shift in shifts]
            index = dists.index(min(dists))
            log_xy_sorted[index] = log_xy[i]
            log_mu_sorted[index] = log_mu[i]

        dropped = 0
        for i, xy in enumerate(log_xy_sorted):
            if xy == (np.nan, np.nan) :
                log_xy_sorted.remove(xy)
                log_mu_sorted.remove(log_mu_sorted[i])
                dropped += 1
        
        if dropped > 1:
            print(f'Dropped {dropped} spots | fx-ID {fraction_id} | beam-ID {beam_id} | layer-ID {layer_id}')
            plt.plot(*zip(*plan_xy), marker='x', ls='-', color='tab:blue', label='plan')
            plt.plot(*zip(*log_xy), marker='o', ls='--', color='tab:grey', label='log')
            plt.plot(*zip(*log_xy_sorted), marker='o', ls='-', color='black', label='sorted')
            plt.legend()
            plt.draw()
            raise ValueError(f'  /!\ Log beam {beam_id} does not match plan beam {beam_dcm.BeamName}, retrying..')
    
        sorted_positions.append(log_xy_sorted), sorted_metersets.append(log_mu_sorted)
    
    # think about returning a dict for re-indexing
    return sorted_positions, sorted_metersets