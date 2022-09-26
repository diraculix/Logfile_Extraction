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
    n_layers = to_be_sorted['TOTAL_LAYERS'].iloc[0]
    sorting_dict = {lid:{} for lid in range(n_layers)}

    # auto-locate plan dicoms
    plan_dcm, beam_ds = self.dicom_finder(fraction_id, beam_id)
    
    # sorting 
    print(f'  Sorting spots for beam-ID {beam_id}..')
    for layer_id in to_be_sorted['LAYER_ID'].drop_duplicates():
        try:
            plan_layer = beam_ds.IonControlPointSequence[layer_id * 2]
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

        log_layer = to_be_sorted.loc[(to_be_sorted['LAYER_ID'] == layer_id) & (to_be_sorted['FRACTION_ID'] == fraction_id)]
        log_xy = [tup for tup in zip(log_layer['X_POSITION(mm)'].to_list(), log_layer['Y_POSITION(mm)'].to_list())]
        log_xy_sorted = [(np.nan, np.nan) for _ in range(max(len(log_xy), len(plan_xy)))]

        # match (x,y)-positions to plan
        for i, log_spot in enumerate(log_xy):
            shifts = [np.array(plan_spot) - np.array(log_spot) for plan_spot in plan_xy]
            dists = [abs((shift).dot(shift)) for shift in shifts]
            index = dists.index(min(dists))
            log_xy_sorted[index] = log_xy[i]

        # remove omitted spots (due to high-weighted tuning)
        dropped = 0
        for xy in log_xy_sorted:
            if xy == (np.nan, np.nan):
                log_xy_sorted.remove(xy)
                dropped += 1
        
        # alert if more than one spot is dropped
        if dropped > 1:
            print(f'Dropped {dropped} spots | fx-ID {fraction_id} | beam-ID {beam_id} | layer-ID {layer_id}')
            plt.plot(*zip(*plan_xy), marker='x', ls='-', color='tab:blue', label='plan')
            plt.plot(*zip(*log_xy), marker='o', ls='--', color='tab:grey', label='log')
            plt.plot(*zip(*log_xy_sorted), marker='o', ls='-', color='black', label='sorted')
            plt.legend()
            plt.draw()
            raise ValueError(f'  /!\ Log beam {beam_id} does not match plan beam {beam_ds.BeamName}, retrying..')

        # define sorting dictionary to be used for position and MU sorting
        sorting_dict[layer_id] = {log_xy.index(log_xy[i]) : log_xy.index(log_xy_sorted[i]) for i in range(len(log_xy))}  # do not iterate over elements directly, index() fails in this case

    return sorting_dict       