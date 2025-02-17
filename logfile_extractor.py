﻿'''
Author:     Lukas C. Wolter, OncoRay ZIK, Dresden, Germany
Project:    Logfile-based patient-specific quality assurance
Encoding:   UTF-8
Synopsis:   The main component of this script is the 'MachineLog'-class,
            which operates on a specific file hierarchy representing the
            patient in terms of treatment fractions, each containing the 
            respective log-files. MachineLog() takes one argument for
            __init()__, which is the parent directory of all fractions
            and has to be submitted to the class by the user. The class
            methods prepare_*() are mandatory prerequisites for further 
            usage, as they compress single logs into one .csv files to 
            enable dataframe operations.
User input: -   root_dir: parent directory containing log-file sorted by
                treatment fraction, submitted to MachineLog-class call
            -   output_dir: directory for output graphs, to be set below
                imports  
            -   self.df_destination: directory for dataframe (compressed
                log-file) storage, to be set in MachineLog.__init__()
'''

import os, sys
import pydicom
import pandas as pd
import numpy as np
import seaborn as sns
from time import sleep
from scipy import optimize
from tkinter import Tk, filedialog
from matplotlib import pyplot as plt

output_dir = r'N:\fs4-HPRT\HPRT-Docs\Lukas\Logfile_Extraction\output'  # set output dir, fallback alternative below (if code is used between multiple systems)
if not os.path.isdir(os.path.join(output_dir, '..')):
    output_dir = r'/home/luke/Scripts/Logfile_Extraction/output'

if not os.path.isdir(output_dir): os.mkdir(output_dir)


'''
Helper functions outside of class, use cases are:
    -   map_spot functions: project PBS spot position recorded from log
        to isocenter (TPS) coordinates, convert charge released in
        nozzle ionizaiton chamber (IC) to monitor units (MU)
    -   fit_sin: returns a correction function for systematic shift in
        spot position (experimental)
'''
def map_spot_pos(pos_arr, ic_offset, sad, ictoiso):
    return np.multiply(np.subtract(pos_arr, ic_offset), np.divide(sad, np.subtract(sad, ictoiso)))


def map_spot_width(width_arr, sad, ictoiso):
    return np.multiply(width_arr, np.divide(sad, np.subtract(sad, ictoiso)))


def map_spot_mu(charge_arr, correction_factor, charge_per_mu):
    return np.divide(np.multiply(charge_arr, correction_factor), charge_per_mu)

def fit_sin(t, y):
    tt, yy = np.array(t), np.array(y)
    
    def sinfunc(t, A, w, p, c):  return A * np.sin(w*t + p) + c

    guess_freq = 1 / 360  # full gantry rotation 360 degrees
    guess_amp = np.std(yy) * 2.**0.5
    guess_offset = np.mean(yy)
    guess = np.array([guess_amp, 2.*np.pi*guess_freq, 0., guess_offset])

    popt, pcov = optimize.curve_fit(sinfunc, tt, yy, p0=guess)
    A, w, p, c = popt
    fitfunc = lambda t: A * np.sin(w*t + p) + c

    return fitfunc


class MachineLog():
    def __init__(self, root_dir):
        valid_dir = False

        # open root directory containing logs, very sensitive to dir structure (only directories must be contained)
        while not valid_dir:
            self.logfile_dir = root_dir
            if self.logfile_dir == '' or str(self.logfile_dir) == '()':
                sys.exit('/x\ No directory selected, exiting..')
            for index, element in enumerate(os.listdir(self.logfile_dir)):
                if not os.path.isdir(os.path.join(self.logfile_dir, element)):
                    print(f'''  /!\ Chosen path '{self.logfile_dir}' should only contain directories (one per fraction). Beware..''')
                if index == len(os.listdir(self.logfile_dir)) - 1:
                    valid_dir = True
                    break
        
        # set dataframe storage directory, fallback alternative below
        # self.df_destination = r'N:\fs4-HPRT\HPRT-Docs\Lukas\Logfile_Extraction\dataframes'
        self.df_destination = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\Logfiles\extracted'
        # self.df_destination = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\04_MA_PSQA_2024\Logfiles\LogfilesExtracted'
        # self.df_destination = r'N:\fs4-HPRT\HPRT-Docs\Lukas\Logfile_Extraction\dataframes'
        if not os.path.isdir(os.path.join(self.df_destination, '..')):            
            self.df_destination = r'/home/luke/Scripts/Logfile_Extraction/dataframes'
        
        if not os.path.isdir(self.df_destination):
            os.mkdir(self.df_destination)

        # initialize fraction and beam list by directory names
        self.fraction_list = sorted(os.listdir(self.logfile_dir))
        self.num_fractions = len(self.fraction_list)
        self.beam_list = []
        for f in self.fraction_list:
            beams_in_frac = []
            if os.path.isdir(os.path.join(self.logfile_dir, f)):
                for dir in sorted(os.listdir(os.path.join(self.logfile_dir, f))):
                    if os.path.isdir(os.path.join(self.logfile_dir, f, dir)):
                        beams_in_frac.append(dir)
            self.beam_list.append(beams_in_frac)

        # scan fraction dirs for missing beams (empy dirs) or invalid beams (e.g. from qa measures)
        self.missing_beams = []  # list of tuples (fraction_id, beam_id)
        for fraction_no, fraction_id in enumerate(self.fraction_list):
            for beam_no, beam_id in enumerate(self.beam_list[fraction_no]):
                current_beam_path = os.path.join(self.logfile_dir, fraction_id, beam_id)
                os.chdir(current_beam_path)
                valid_beam = True
                while len(os.listdir('.')) <= 3:  # log-file dirs are nested
                    try:
                        os.chdir(os.listdir('.')[0])
                    except:
                        print(f'  /x\ No log-files detected for beam {beam_id} in fraction {fraction_id}, continuing..')
                        self.missing_beams.append((fraction_id, beam_id))
                        valid_beam = False
                        break
                
                if not valid_beam: continue

                for file in os.listdir('.'):
                    if file.__contains__('beam.'):
                        beam_file = file

                try:
                    with open(beam_file, 'r') as beam_file:
                        for line in beam_file.readlines():
                                if line.__contains__('PatientId'):
                                    self.patient_id = line.split('>')[1].split('<')[0]
                        
                        beam_file.close()
                except TypeError:
                    print(f'  /!\ Invalid beam {beam_id} in fraction {fraction_id}')
                    continue

        # check for already existent dataframe for patient, read in and prevent re-creation in this case
        record_df_exists, tuning_df_exists = False, False
        for fname in os.listdir(self.df_destination):  # read existing dataframes
            if fname.__contains__(f'{self.patient_id}_records') and fname.endswith('.csv'):
                self.record_df_name = fname
                print(f'''Found patient record dataframe '{self.record_df_name}', reading in..''')
                self.patient_record_df = pd.read_csv(os.path.join(self.df_destination, fname), parse_dates=['TIME'], index_col='TIME', dtype={'BEAM_ID':str, 'FRACTION_ID':str})
                record_df_exists = True
            elif fname.__contains__(f'{self.patient_id}_tuning') and fname.endswith('.csv'):
                self.tuning_df_name = fname
                print(f'''Found patient tuning dataframe '{self.tuning_df_name}', reading in..''')
                self.patient_tuning_df = pd.read_csv(os.path.join(self.df_destination, fname), parse_dates=['TIME'], index_col='TIME', dtype={'BEAM_ID':str, 'FRACTION_ID':str})
                tuning_df_exists = True
            
            if record_df_exists and tuning_df_exists: break

        # if no patient dataframes are found, create them by default
        if not record_df_exists or not tuning_df_exists:
            print(f'''\nUnable to locate patient record/tuning dataframes for patient-ID {self.patient_id}. Please run prepare_dataframe()..''')
            self.patient_record_df, self.patient_tuning_df = pd.DataFrame(), pd.DataFrame()
        
        os.chdir(self.logfile_dir)    


    '''
    Main data extraction method:
    Requires:   Initialized, valid log-file directory
    Arguments:  None
    Operation:  Loop over fraction directories, collect relevant beam data in pandas dataframe, write to .csv
    Returns:    None
    '''
    def prepare_dataframe(self) -> None:
        if not self.patient_record_df.empty and not self.patient_tuning_df.empty:  # function call obsolete if record/tuning dataframes exist, re-init possible
            print('Already read existing patient record/tuning dataframes:')
            print(f'  {self.record_df_name}\n  {self.tuning_df_name}\n')
            re_init = input('Re-initialize patient dataframe [y/n]? ')
            if re_init != 'y':
                return None
        
        # try:  # set up function for gtr-angle-dependent spot position correction (based on measured QA data from 2017-2022)
        #     qa_df = pd.read_csv(f'{output_dir}/QA_angular_dependence.csv')
        #     qa_angles = qa_df.GANTRY_ANGLE.to_list()
        #     qa_angles.append(360.)
        #     qa_median_x = qa_df['DELTA_X[mm]'].to_list()
        #     qa_median_x.append(qa_median_x[0])
        #     qa_median_y = qa_df['DELTA_Y[mm]'].to_list()
        #     qa_median_y.append(qa_median_y[0])
        #     correct_x = fit_sin(qa_angles, qa_median_x)  # correct_x/y are functions
        #     correct_y = fit_sin(qa_angles, qa_median_y)
        # except FileNotFoundError:
        #     print('  /!\ Spot position correction data not found, no funciton will be applied..')
        #     def temp(x): return x
        #     correct_x = temp, correct_y = temp  # <-- function

        # sweep over all fractions and beams, extract information
        print(f'Initializing dataframe for patient-ID {self.patient_id}..')
        self.patient_record_df, self.patient_tuning_df = pd.DataFrame(), pd.DataFrame()  # empty target dataframes
        for fraction_no, fraction_id in enumerate(self.fraction_list):
            for beam_id in self.beam_list[fraction_no]:
                if (fraction_id, beam_id) in self.missing_beams:
                    print(f'  /!\ Skipping missing beam {beam_id} in fraction {fraction_id}')
                    continue

                current_beam_path = os.path.join(self.logfile_dir, fraction_id, beam_id)
                os.chdir(current_beam_path)
                while len(os.listdir('.')) <= 3:
                    try:
                        os.chdir(sorted(os.listdir('.'))[0])  # navigate through nested logfile dir structure
                    except OSError:
                        print(f'  /!\ No directory to enter in {os.getcwd()}, staying here..')
                        break
                
                # identify files in beam directory
                map_records, tunings, record_specifs, tuning_specifs = [], [], [], []
                for file in sorted(os.listdir('.')):
                    if file.__contains__('beam.'):
                        beam_file = file
                    elif file.__contains__('beam_config.'):
                        beam_config = file
                    elif file.__contains__('map_record') and file.__contains__('part') and not file.__contains__('temp'):
                        map_records.append(file)
                    elif file.__contains__('map_record') and file.__contains__('tuning') and not file.__contains__('temp'):
                        tunings.append(file)
                    elif file.__contains__('map_specif') and file.__contains__('part'):
                        record_specifs.append(file)
                    elif file.__contains__('map_specif') and file.__contains__('tuning'):
                        tuning_specifs.append(file)
                    
                if len(map_records) == 0 or len(tunings) == 0:
                        print(f'  /x\ No logfiles found for beam-ID {beam_id} in fraction {fraction_id}')
                        continue

                num_layers = max([int(fname.split('_')[2].split('_')[0]) + 1 for fname in tunings])  # every planned layer must have a tuning file, but not necessarily a record file

                # draw beam specs from *_beam.csv
                with open(beam_file, 'r') as beam_file:  
                    for line in beam_file.readlines():
                        if line.__contains__('GantryAngle'):
                            gantry_angle = float(line.split('>')[1].split('<')[0])
                        elif line.__contains__('pressure'):
                            pressure = float(line.split(',')[-1])
                        elif line.__contains__('temperature'):
                            temperature = float(line.split(',')[-1])
                        elif line.__contains__('doseCorrectionFactor'):
                            chamber_correction = float(line.split(',')[-1])

                    beam_file.close()

                # draw machine parameters from *beam_config.csv
                with open(beam_config, 'r') as beam_config:
                    iso_disp_x, iso_disp_y = [], []
                    for no, line in enumerate(beam_config.readlines()):
                        if line.__contains__('SAD parameter'):
                            sad_x, sad_y = float(line.split(';')[-1].split(',')[1]), float(line.split(';')[-1].split(',')[0])  # coordinate system switch x <--> y
                        elif line.__contains__('distanceFromIcToIsocenter'):
                            ictoiso_x, ictoiso_y = float(line.split(';')[-1].split(',')[1]), float(line.split(';')[-1].split(',')[0])  # coordinate system switch x <--> y
                        elif line.__contains__('chargePerMUIC2'):
                            charge_per_mu = float(line.split(';')[-1])
                        elif line.__contains__('Nozzle WET polynomial'):
                            nozzle_wet_poly = [float(x) for x in line.split(';')[-1].split(',')]
                        
                    ref_pressure, ref_temperature = 1013., 293.  # [hPa, K] standard reference, can also be read from same file
                    correction_factor = (1 / chamber_correction) * (ref_pressure * temperature) / (ref_temperature * pressure)  # air density correction for IC readout

                    beam_config.close()
                
                # Source: IBA Particle Therapy GmbH 08/22 (JBO), UPTD universal nozzle WET polynomial coefficients
                iba_gtr2_poly = [0.001684756748152, -0.00490089228886989, 0.561372013469097, 3.46404838890297]

                # read in irradiation records layer-wise, process, then append to final dataframe
                layer_exceptions, finalized_layers, finalized_tunings = [], [self.patient_record_df], [self.patient_tuning_df]
                for layer_id in range(num_layers):
                    to_do_layers, to_do_tunings = [], []
                    no_exceptions = True
                    for rec_idx, record_file in enumerate(map_records):  # extraction of unpacked log-file data
                        if int(record_file.split('_')[2]) == layer_id:
                            try:
                                record_file_df = pd.read_csv(record_file, delimiter=',', skiprows=10, skipfooter=11, engine='python')

                            except OSError:
                                print(f'  /!\ OSError: Check that {record_file} path is valid')
                                print(f'      File exists -> {os.path.exists(os.path.join(os.getcwd(), record_file))}')
                                for s in range(5):
                                    print(f'      Retrying in {4 - s}s ...', end='\r')
                                    sleep(1)

                            except pd.errors.ParserError:  # log-file may not be read if interlock caused restart of logging in same file 
                                print('  /!\ Read CSV error:', record_file)
                                print('      Cleaning record file..')

                                # workaround: manually split record file in temporary files, then add them to queue
                                with open(record_file, 'r') as record:
                                    lines, splits = record.readlines(), []
                                    for nr, line in enumerate(lines):
                                        if line.__contains__('beamline') and lines[nr - 1] in ['\n', '\r\n']:
                                            splits.append(nr)
                                    record.close()
                                
                                temp_record_files = []
                                for nr, split_at in enumerate(splits):
                                    temp_record_name = record_file.split('.')[:-1]
                                    temp_record_name.append(f'temp_{nr + 1}.csv')
                                    temp_record_name = '_'.join(temp_record_name)
                                    with open(temp_record_name, 'w+') as temp_record:  # initialize temporary record file
                                        if nr == 0:
                                            temp_record.writelines(lines[:splits[nr + 1]])
                                        else:
                                            try:
                                                temp_record.writelines(lines[split_at:splits[nr + 1]])
                                            except:
                                                temp_record.writelines(lines[split_at:])
                                        
                                        temp_record_files.append(temp_record.name)
                                        temp_record.close()
                                
                                for i, temp_record_file in enumerate(temp_record_files):
                                    if temp_record_file.__contains__('_temp_'):
                                        map_records.insert(rec_idx + (i+1), temp_record_file)
                                    
                                continue  # move on to appended file in queue
                                        
                            try:
                                record_file_df['TIME'] = pd.to_datetime(record_file_df['TIME'], dayfirst=True)  # datetime index --> chronological order
                                charge_col = pd.Series(record_file_df['DOSE_PRIM(C)'])                          # ion dose [C], to be converted in MU
                                current_col = pd.Series(record_file_df['DOSE_RATE_PRIM(A)'])                       # beam current, measured as voltage at ... ?
                                record_file_df = record_file_df.loc[:, :'Y_POSITION(mm)']                       # slice dataframe, drop redundant columns
                                record_file_df['DOSE_PRIM(C)'] = charge_col
                                record_file_df['BEAMCURRENT(A)'] = current_col
                                try:
                                    record_file_df.drop(record_file_df[record_file_df['SUBMAP_NUMBER'] < 0].index, inplace=True)  # if split record file, nothing serious
                                except: pass

                                record_file_df = record_file_df[record_file_df.groupby('SUBMAP_NUMBER')['SUBMAP_NUMBER'].transform('count') > 1]  # drop all rows without plan-relevant data
                            
                            except:  # if unusable information in log-file (possible if split into multiple files due to interlock)
                                no_exceptions = False
                                layer_exceptions.append(layer_id)
                                continue
                            
                            if record_file_df.loc[(record_file_df['X_POSITION(mm)'] != -10000.) & (record_file_df['Y_POSITION(mm)'] != -10000.)].empty:  # in case of split record file containing only invalid entries
                                no_exceptions = False
                                layer_exceptions.append(layer_id)
                                continue

                            # submap number: running index of table --> enables filtering of actual irradiations
                            current_spot_submap = record_file_df['SUBMAP_NUMBER'].min()
                            current_spot_id = 0
                            record_file_df['SPOT_ID'] = 0
                            previous_xpos, previous_ypos, previous_spot_submap = None, None, None
                            record_file_df.reindex()

                            # sweep over all spots in layer, SUBMAP_NUMBER is locked whenever a spot is active
                            while current_spot_submap <= record_file_df['SUBMAP_NUMBER'].max():
                                split_spot = False

                                # get spot drill time and released charge
                                spot_drill_time = len(record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == current_spot_submap]) * 0.25  # in [ms]
                                accumulated_charge = record_file_df.loc[(record_file_df['SUBMAP_NUMBER'] == current_spot_submap) & (record_file_df['DOSE_PRIM(C)'] != -10000.0), ['DOSE_PRIM(C)']].abs().sum().iloc[0]
                                if accumulated_charge == 0.:
                                    print('  /!\ ZERO MU ERROR @ spot', current_spot_id, '(', current_spot_submap, ') in', record_file)
                                
                                # drop unusable rows AFTER charge extraction
                                record_file_df.drop(record_file_df.loc[(record_file_df['SUBMAP_NUMBER'] == current_spot_submap) & ((record_file_df['X_POSITION(mm)'] == -10000.0) | (record_file_df['Y_POSITION(mm)'] == -10000.0))].index, inplace=True)  # drop unusable rows for current submap
                                
                                # average over all spot entries for most accurate position/shape (recommended by IBA)                                
                                mean_xpos = record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['X_POSITION(mm)']].mean().iloc[0]
                                mean_ypos = record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['Y_POSITION(mm)']].mean().iloc[0]
                                beam_curr = record_file_df.loc[(record_file_df['SUBMAP_NUMBER'] == current_spot_submap) & (record_file_df['BEAMCURRENT(A)'] > 0), ['BEAMCURRENT(A)']].median().iloc[0]     # last entry always drops to 0, distorts mean
                                                                                                                                                                # current is also reduced (~0.5) in case of light (<4ms) spots
                                if current_spot_id != 0:
                                    if abs(mean_xpos - previous_xpos) < 1 and abs(mean_ypos - previous_ypos) < 1:  # proximity check for split spots (in case of high MU), custom tolerance 1mm
                                        print(f'  /!\ Split spot detected in layer-ID {layer_id}, merging..')
                                        split_spot = True

                                # in case of split spot, add time and charge to previous spot
                                if split_spot:
                                    record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == previous_spot_submap, ['DRILL_TIME(ms)']] += spot_drill_time
                                    record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == previous_spot_submap, ['CHARGE(C)']] += accumulated_charge
                                    record_file_df.drop(record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == current_spot_submap].index, inplace=True)
                                    record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == previous_spot_submap, ['SUBMAP_NUMBER']] = current_spot_submap
                                else:
                                    record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['SPOT_ID']] = current_spot_id  # assign new column
                                    record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['DRILL_TIME(ms)']] = spot_drill_time  # assign new column for drill time
                                    record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['CHARGE(C)']] = accumulated_charge  # accumulate charge released per spot
                                    record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['X_POS_IC23(mm)']] = mean_xpos
                                    record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['Y_POS_IC23(mm)']] = mean_ypos
                                    record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['X_WID_IC23(mm)']] = record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['X_WIDTH(mm)']].mean().iloc[0]
                                    record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['Y_WID_IC23(mm)']] = record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['Y_WIDTH(mm)']].mean().iloc[0]
                                    record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['BEAM_CURRENT(A)']] = beam_curr

                                previous_xpos, previous_ypos, previous_spot_submap = mean_xpos, mean_ypos, current_spot_submap
                                current_spot_submap = record_file_df.loc[record_file_df['SUBMAP_NUMBER'] > current_spot_submap]['SUBMAP_NUMBER'].min()  # proceed to next submap
                                if not split_spot: current_spot_id += 1  # keep track of spot id
                            
                            record_file_df.drop(columns=['DOSE_PRIM(C)', 'BEAMCURRENT(A)'], inplace=True)
                            record_file_df.drop_duplicates(subset=['SUBMAP_NUMBER'], keep='last', inplace=True)  # keep only last entry for each spot
                            record_file_df.index = record_file_df['TIME']  # change to datetime index AFTER filtering, timestamp is NOT unique!
                            record_file_df.drop(columns=['TIME', 'SUBMAP_NUMBER'], inplace=True)
                            record_file_df.reindex()

                            # draw machine parameters from *map_specif*.csv
                            for specif_file in record_specifs:  
                                if int(specif_file.split('_')[2].split('_')[0]) == layer_id:
                                    try:
                                        with open(specif_file, 'r') as specif_file:
                                            lines = specif_file.readlines()
                                            ic_offsets = lines[3]
                                            ic_offset_x, ic_offset_y = float(ic_offsets.split(',')[3]), float(ic_offsets.split(',')[2])  # coordinate system switch x <--> y
                                            range_at_degrader = float(lines[1].split(',')[1])  # beam range after energy selection system
                                    except OSError:
                                        print(f'  /!\ OSError: Check that {specif_file} path is valid')
                                        print(f'      File exists -> {os.path.exists(os.path.join(os.getcwd(), specif_file))}')
                            
                            # calculate energy at isocenter from beam range
                            nozzle_wet = np.polyval(nozzle_wet_poly, range_at_degrader)  # [mm]
                            range_at_iso = range_at_degrader - nozzle_wet
                            layer_energy = np.exp(np.polyval(iba_gtr2_poly, np.log(range_at_iso)))  # [MeV]
                            record_file_df['LAYER_ENERGY(MeV)'] = layer_energy  # set new column

                            # coordinate system transform IBA <-> RayStation/ISO (x <-> y)
                            record_file_df['X_POSITION(mm)'] = record_file_df[['Y_POS_IC23(mm)']].apply(map_spot_pos, args=(ic_offset_x, sad_x, ictoiso_x))
                            record_file_df['Y_POSITION(mm)'] = record_file_df[['X_POS_IC23(mm)']].apply(map_spot_pos, args=(ic_offset_y, sad_y, ictoiso_y))
                            # record_file_df['X_POSITION_CORR(mm)'] = record_file_df['X_POSITION(mm)'] - correct_x(gantry_angle)  # if correction data is present (experimental)
                            # record_file_df['Y_POSITION_CORR(mm)'] = record_file_df['Y_POSITION(mm)'] - correct_y(gantry_angle)
                            record_file_df['X_WIDTH(mm)'] = record_file_df[['Y_WID_IC23(mm)']].apply(map_spot_width, args=(sad_x, ictoiso_x))
                            record_file_df['Y_WIDTH(mm)'] = record_file_df[['X_WID_IC23(mm)']].apply(map_spot_width, args=(sad_y, ictoiso_y))
                            # record_file_df['SQDIST_TO_ISO(mm)'] = np.square(record_file_df['X_POSITION(mm)']) + np.square(record_file_df['Y_POSITION(mm)'])
                        
                            # charge to MU conversion using correction factor
                            record_file_df['CORRECTION_FACTOR'] = correction_factor
                            record_file_df['MU'] = record_file_df[['CHARGE(C)']].apply(map_spot_mu, args=(correction_factor, charge_per_mu))
                            record_file_df.reindex()  # make sure modified layer df is consistent with indexing

                            # in the case of split record file
                            if len(to_do_layers) > 0:
                                previous_record_file_df = to_do_layers[-1]
                                prev_x_last, prev_y_last = previous_record_file_df['X_POSITION(mm)'].iloc[-1], previous_record_file_df['Y_POSITION(mm)'].iloc[-1]
                                this_x_first, this_y_first = record_file_df['X_POSITION(mm)'].iloc[0], record_file_df['Y_POSITION(mm)'].iloc[0]

                                # check proximity to last spot entry of previous file (if same coordinates, spot is split over 2 files)
                                if abs(this_x_first - prev_x_last) < 1 and abs(this_y_first - prev_y_last) < 1:  
                                    print(f'''  /!\ Spot-ID {previous_record_file_df['SPOT_ID'].max()} of layer-ID {layer_id} split over 2 files, merging..''')
                                    this_charge, this_dose, this_drill = record_file_df['CHARGE(C)'].iloc[0] ,record_file_df['MU'].iloc[0], record_file_df['DRILL_TIME(ms)'].iloc[0]
                                    # add dose and time to previous
                                    previous_record_file_df.at[previous_record_file_df.index[-1], 'CHARGE(C)'] += this_charge
                                    previous_record_file_df.at[previous_record_file_df.index[-1], 'MU'] += this_dose  
                                    previous_record_file_df.at[previous_record_file_df.index[-1], 'DRILL_TIME(ms)'] += this_drill
                                    to_do_layers[-1] = previous_record_file_df  # modify element in queue
                                    record_file_df = record_file_df.iloc[1:]  # drop first entry of current file
                                    record_file_df['SPOT_ID'] -= 1  # consequence: reduce running spot-id of current file

                            if not record_file_df.empty: to_do_layers.append(record_file_df)  # if only spot was deleted in previous step
                    
                    # same procedure for all tuning spots (max. 3 per layer possible)
                    for tune_idx, tuning_file in enumerate(tunings): 
                        if int(tuning_file.split('_')[2]) == layer_id:
                            try:
                                tuning_file_df = pd.read_csv(tuning_file, delimiter=',', skiprows=10, skipfooter=11, engine='python')
                            except:
                                print('  /!\ Read CSV error:', tuning_file)
                                print('      Cleaning tuning file..')
                                with open(tuning_file, 'r') as tuning:
                                    lines, splits = tuning.readlines(), []
                                    for nr, line in enumerate(lines):
                                        if line.__contains__('beamline') and lines[nr - 1] in ['\n', '\r\n']:
                                            splits.append(nr)
                                    tuning.close()
                                
                                temp_tuning_files = []
                                for nr, split_at in enumerate(splits):
                                    temp_tuning_name = tuning_file.split('.')[:-1]
                                    temp_tuning_name.append(f'temp_{nr + 1}.csv')
                                    temp_tuning_name = '_'.join(temp_tuning_name)
                                    with open(temp_tuning_name, 'w+') as temp_tuning:
                                        if nr == 0:
                                            temp_tuning.writelines(lines[:splits[nr + 1]])
                                        else:
                                            try:
                                                temp_tuning.writelines(lines[split_at:splits[nr + 1]])
                                            except:
                                                temp_tuning.writelines(lines[split_at:])
                                        
                                        temp_tuning_files.append(temp_tuning.name)
                                        temp_tuning.close()
                                
                                for i, temp_tuning_file in enumerate(temp_tuning_files):
                                    if temp_tuning_file.__contains__('_temp_'):
                                        tunings.insert(tune_idx + (i+1), temp_tuning_file)
                                
                                continue

                            try:
                                tuning_file_df['TIME'] = pd.to_datetime(tuning_file_df['TIME'], dayfirst=True)
                                charge_col = pd.Series(tuning_file_df['DOSE_PRIM(C)'])
                                tuning_file_df = tuning_file_df.loc[:, :'Y_POSITION(mm)']
                                tuning_file_df['DOSE_PRIM(C)'] = charge_col
                                try:
                                    tuning_file_df.drop(tuning_file_df[tuning_file_df['SUBMAP_NUMBER'] < 0].index, inplace=True)
                                except: pass
                                    
                                tuning_file_df = tuning_file_df[tuning_file_df.groupby('SUBMAP_NUMBER')['SUBMAP_NUMBER'].transform('count') > 1]
                                
                            except KeyError:
                                print(f'''\n  /!\ Key error occured while handling '{tuning_file}' (layer {str(layer_id).zfill(2)}), continuing..''')
                                no_exceptions = False
                                layer_exceptions.append(layer_id)
                                continue

                            if tuning_file_df.empty:  # in case of low-weighted tuning spots, no usable information will be left in log-file, skip these occasions
                                no_exceptions = False
                                layer_exceptions.append(layer_id)
                                continue

                            current_spot_submap = tuning_file_df['SUBMAP_NUMBER'].min()
                            current_spot_id = 0
                            tuning_file_df['SPOT_ID'] = 0
                            tuning_file_df.reindex()
                            
                            while current_spot_submap <= tuning_file_df['SUBMAP_NUMBER'].max():
                                tuning_file_df.loc[tuning_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['SPOT_ID']] = current_spot_id  # assign new column
                                spot_drill_time = len(tuning_file_df.loc[tuning_file_df['SUBMAP_NUMBER'] == current_spot_submap]) * 0.25  # in [ms]
                                tuning_file_df.loc[tuning_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['DRILL_TIME(ms)']] = spot_drill_time  # assign new column for drill time
                                accumulated_charge = tuning_file_df.loc[(tuning_file_df['SUBMAP_NUMBER'] == current_spot_submap) & (tuning_file_df['DOSE_PRIM(C)'] != -10000.0), ['DOSE_PRIM(C)']].abs().sum().iloc[0]
                                tuning_file_df.loc[tuning_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['CHARGE(C)']] = accumulated_charge  # accumulate charge released per spot
                                
                                # average over all spot entries for most accurate position/shape (recommended by IBA)
                                tuning_file_df.drop(tuning_file_df.loc[(tuning_file_df['SUBMAP_NUMBER'] == current_spot_submap) & ((tuning_file_df['X_POSITION(mm)'] == -10000.0) | (tuning_file_df['Y_POSITION(mm)'] == -10000.0))].index, inplace=True)  # drop unusable rows for current submap
                                tuning_file_df.loc[tuning_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['X_POS_IC23(mm)']] = tuning_file_df.loc[tuning_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['X_POSITION(mm)']].mean().iloc[0]
                                tuning_file_df.loc[tuning_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['Y_POS_IC23(mm)']] = tuning_file_df.loc[tuning_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['Y_POSITION(mm)']].mean().iloc[0]
                                tuning_file_df.loc[tuning_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['X_WID_IC23(mm)']] = tuning_file_df.loc[tuning_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['X_WIDTH(mm)']].mean().iloc[0]
                                tuning_file_df.loc[tuning_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['Y_WID_IC23(mm)']] = tuning_file_df.loc[tuning_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['Y_WIDTH(mm)']].mean().iloc[0]

                                current_spot_submap = tuning_file_df.loc[tuning_file_df['SUBMAP_NUMBER'] > current_spot_submap]['SUBMAP_NUMBER'].min()  # proceed to next submap
                                current_spot_id += 1  # keep track of spot id

                            tuning_file_df.drop(columns=['DOSE_PRIM(C)'], inplace=True)
                            tuning_file_df.drop_duplicates(subset=['SUBMAP_NUMBER'], keep='last', inplace=True)  # keep only last entry for each spot
                            tuning_file_df.index = tuning_file_df['TIME']
                            tuning_file_df.drop(columns=['TIME', 'SUBMAP_NUMBER'], inplace=True)
                            
                            for specif_file in tuning_specifs:
                                if int(specif_file.split('_')[2].split('_')[0]) == layer_id:
                                    with open(specif_file, 'r') as specif_file:
                                        lines = specif_file.readlines()
                                        ic_offsets = lines[3]
                                        ic_offset_x, ic_offset_y = float(ic_offsets.split(',')[3]), float(ic_offsets.split(',')[2])
                                        range_at_degrader = float(lines[1].split(',')[1])

                            nozzle_wet = np.polyval(nozzle_wet_poly, range_at_degrader)  # [mm]
                            range_at_iso = range_at_degrader - nozzle_wet
                            layer_energy = np.exp(np.polyval(iba_gtr2_poly, np.log(range_at_iso)))  # [MeV]
                            tuning_file_df['LAYER_ENERGY(MeV)'] = layer_energy

                            # tuning_file_df['X_POS_IC23(mm)'], tuning_file_df['Y_POS_IC23(mm)'] = tuning_file_df['X_POSITION(mm)'], tuning_file_df['Y_POSITION(mm)']
                            # tuning_file_df['X_WID_IC23(mm)'], tuning_file_df['Y_WID_IC23(mm)'] = tuning_file_df['X_WIDTH(mm)'], tuning_file_df['Y_WIDTH(mm)']
                            tuning_file_df['X_POSITION(mm)'] = tuning_file_df[['Y_POS_IC23(mm)']].apply(map_spot_pos, args=(ic_offset_x, sad_x, ictoiso_x))
                            tuning_file_df['Y_POSITION(mm)'] = tuning_file_df[['X_POS_IC23(mm)']].apply(map_spot_pos, args=(ic_offset_y, sad_y, ictoiso_y))
                            tuning_file_df['X_WIDTH(mm)'] = tuning_file_df[['Y_WID_IC23(mm)']].apply(map_spot_width, args=(sad_x, ictoiso_x))
                            tuning_file_df['Y_WIDTH(mm)'] = tuning_file_df[['X_WID_IC23(mm)']].apply(map_spot_width, args=(sad_y, ictoiso_y))
                            tuning_file_df['MU'] = tuning_file_df[['CHARGE(C)']].apply(map_spot_mu, args=(correction_factor, charge_per_mu))

                            tuning_file_df.reindex()
                            to_do_tunings.append(tuning_file_df)
                    
                    # in case of multiple layer parts: enable continuous spot indexing
                    for i in range(len(to_do_tunings)):
                        if i > 0:
                            to_do_tunings[i]['SPOT_ID'] += (to_do_tunings[i - 1]['SPOT_ID'].max() + 1)
                    for j in range(len(to_do_layers)):
                        if j > 0:
                            to_do_layers[j]['SPOT_ID'] += (to_do_layers[j - 1]['SPOT_ID'].max() + 1)
                    
                    # concatenate layers, assign additional columns
                    omitted_by_tuning = False

                    if len(to_do_layers) > 0:  # can be zero, if only one spot in layer and omitted by high-weighted tuning
                        layer_df = pd.concat(to_do_layers)
                        layer_df['LAYER_ID'] = layer_id
                        layer_df['TOTAL_LAYERS'] = num_layers
                        layer_df['BEAM_ID'] = beam_id
                        layer_df['GANTRY_ANGLE'] = gantry_angle
                        layer_df['TEMPERATURE(K)'] = temperature
                        layer_df['PRESSURE(hPa)'] = pressure
                        layer_df['FRACTION_ID'] = fraction_id
                        layer_df['PATIENT_ID'] = self.patient_id
                        layer_df = layer_df[~layer_df.index.duplicated(keep='first')]
                    else:
                        print(f'  /!\ No record for layer-ID {layer_id} in beam {beam_id}, only spot replaced by tuning')
                        omitted_by_tuning = True
                                   
                    if len(to_do_tunings) > 0:
                        tuning_df = pd.concat(to_do_tunings)
                        tuning_df['LAYER_ID'] = layer_id
                        tuning_df['TOTAL_LAYERS'] = num_layers
                        tuning_df['BEAM_ID'] = beam_id
                        tuning_df['GANTRY_ANGLE'] = gantry_angle
                        tuning_df['TEMPERATURE(K)'] = temperature
                        tuning_df['PRESSURE(hPa)'] = pressure
                        tuning_df['FRACTION_ID'] = fraction_id
                        tuning_df['PATIENT_ID'] = self.patient_id
                        tuning_df = tuning_df[~tuning_df.index.duplicated(keep='first')]
                    else:
                        print(f'  /!\ No tunings found for layer-ID {layer_id} and beam {beam_id}, skipping this layer..')
                        continue
                    
                    del to_do_layers, to_do_tunings

                    if not omitted_by_tuning:
                        finalized_layers.append(layer_df)
                    finalized_tunings.append(tuning_df)

                    if no_exceptions:
                        char = '#'
                    else:
                        char = '_'

                    if layer_id == (num_layers - 1):  # progress visualization
                        if no_exceptions:
                            print('  ', '[' + (layer_id + 1) * char + (num_layers - layer_id - 1) * '-' + ']', end=f' Beam {beam_id} complete\n')
                        else:    
                            print('  ', '[' + (layer_id + 1) * char + (num_layers - layer_id - 1) * '-' + ']', end=f' Beam {beam_id} complete (empty dataframe exception in layer(s) {layer_exceptions}\n')
                    else:
                        print('  ', '[' + (layer_id + 1) * char + (num_layers - layer_id - 1) * '-' + ']', end=f' Layer {str(layer_id + 1).zfill(2)}/{str(num_layers).zfill(2)}\r')
                    
                    no_exceptions = True

                # remove temporary files
                for file in os.listdir('.'):
                    if file.__contains__('temp'):
                        os.remove(file)

                # concatenate single layer dataframes
                self.patient_record_df = pd.concat(finalized_layers, sort=True)
                self.patient_tuning_df = pd.concat(finalized_tunings, sort=True)
                os.chdir(self.logfile_dir)

            print(f'  ..Fraction {str(fraction_no + 1).zfill(2)}/{str(self.num_fractions).zfill(2)} complete..')

        # write out dataframe as .csv
        os.chdir(self.df_destination)
        self.record_df_name = f'patient_{self.patient_id}_records_data.csv'
        self.tuning_df_name = f'patient_{self.patient_id}_tuning_data.csv'

        # drop unnecessary cols
        self.patient_record_df.drop(columns=['CHARGE(C)', 'DRILL_TIME(ms)', 'X_POS_IC23(mm)', 'Y_POS_IC23(mm)', 'X_WIDTH(mm)', 'Y_WIDTH(mm)', 'X_WID_IC23(mm)', 'Y_WID_IC23(mm)'], inplace=True)
        self.patient_tuning_df.drop(columns=['CHARGE(C)', 'DRILL_TIME(ms)', 'X_POS_IC23(mm)', 'Y_POS_IC23(mm)', 'X_WIDTH(mm)', 'Y_WIDTH(mm)', 'X_WID_IC23(mm)', 'Y_WID_IC23(mm)'], inplace=True)
        
        print(f'''  ..Writing dataframe to '{self.df_destination}' as .CSV.. ''')
        while True:   
            try:
                self.patient_record_df.to_csv(self.record_df_name)
                self.patient_tuning_df.to_csv(self.tuning_df_name)
                break
            except PermissionError:
                input('  Permission denied, close target file and press ENTER.. ')

        # adjust column types if script runs further 
        self.patient_record_df = self.patient_record_df.astype(dtype={'BEAM_ID':str, 'FRACTION_ID':str})
        self.patient_tuning_df = self.patient_tuning_df.astype(dtype={'BEAM_ID':str, 'FRACTION_ID':str})
        print('Complete')


    '''
    Single-purpose function to filter for and extract log-file data from QA measurements, adopted from prepare_dataframe():
    Requires:   Initialized valid log-file directory
    Arguments:  None
    Operation:  Restricted log-file data mining using only allowed beam parameters from QA protocol, collect and write to .csv
    Returns:    None
    '''
    def prepare_qa_dataframe(self) -> None:
        # allowed beam parameters
        qa_energies = [100., 120., 140., 165., 185., 205., 226.7]
        qa_angles = np.linspace(0., 360., 8, endpoint=False)
        qa_xy = [0., 30., 60., 120.]
        qa_xy += [- xy for xy in qa_xy]

        print(f'\nMining biannual spot-QA data..')
        self.qa_record_df = pd.DataFrame()  # overwrite stored df's
        for fraction_no, fraction_id in enumerate(self.fraction_list):
            num_beams = len(self.beam_list[fraction_no])
            for beam_no, beam_id in enumerate(self.beam_list[fraction_no]):
                current_beam_path = os.path.join(self.logfile_dir, fraction_id, beam_id)
                os.chdir(current_beam_path)
                while len(os.listdir('.')) <= 3:
                    try:
                        os.chdir(sorted(os.listdir('.'))[0])  # navigate through nested logfile dir structure
                    except OSError:
                        print(f'  /!\ No directory to enter in {os.getcwd()}, staying here..')
                        break

                map_records, record_specifs = [], []  # get file lists
                for file in sorted(os.listdir('.')):
                    if file.__contains__('beam.'):
                        beam_file = file
                    elif file.__contains__('beam_config.'):
                        beam_config = file
                    elif file.__contains__('map_record') and not file.__contains__('tuning') and not file.__contains__('temp'):
                        map_records.append(file)
                    elif file.__contains__('map_specif') and not file.__contains__('tuning') and not file.__contains__('temp'):
                        record_specifs.append(file)
                    
                if len(map_records) == 0:
                    print(f'  /!\ Skipping fraction-ID {fraction_id} beam-ID {beam_id}')
                    continue
                                
                num_layers = max([int(fname.split('_')[2].split('.')[0].split('_')[0]) + 1 for fname in map_records])

                with open(beam_file, 'r') as beam_file:  # draw beam specs from *_beam.csv
                    for line in beam_file.readlines():
                        if line.__contains__('GantryAngle'):
                            gantry_angle = float(line.split('>')[1].split('<')[0])
                            break

                    beam_file.close()
                
                if not np.round(gantry_angle, 0) in qa_angles:
                    continue

                with open(beam_config, 'r') as beam_config:  # draw machine parameters from *beam_config.csv
                    for line in beam_config.readlines():
                        if line.__contains__('SAD parameter'):
                            sad_x, sad_y = float(line.split(';')[-1].split(',')[1]), float(line.split(';')[-1].split(',')[0])  # coordinate system switch x <--> y
                        elif line.__contains__('distanceFromIcToIsocenter'):
                            ictoiso_x, ictoiso_y = float(line.split(';')[-1].split(',')[1]), float(line.split(';')[-1].split(',')[0])  # coordinate system switch x <--> y
                        elif line.__contains__('Nozzle WET polynomial'):
                            nozzle_wet_poly = [float(x) for x in line.split(';')[-1].split(',')]
                
                # Source: IBA Particle Therapy 08/22 (Jozef Bokor), universal nozzle WET polynomial coefficients
                iba_gtr2_poly = [0.001684756748152, -0.00490089228886989, 0.561372013469097, 3.46404838890297]
                
                layer_exceptions, finalized_layers= [], [self.qa_record_df]
                for layer_id in range(num_layers):
                    to_do_layers = []
                    no_exceptions = True
                    for file_id, record_file in enumerate(map_records):  # actual (processed) log-file extraction
                        if int(record_file.split('_')[2].split('.')[0].split('_')[0]) == layer_id:
                            try:
                                record_file_df = pd.read_csv(record_file, delimiter=',', skiprows=10, skipfooter=11, engine='python')

                            except:
                                print('  /!\ Read CSV error:', record_file)
                                print('      Cleaning record file..')
                                with open(record_file, 'r') as record:
                                    lines, splits = record.readlines(), []
                                    for nr, line in enumerate(lines):
                                        if line.__contains__('beamline') and lines[nr - 1] in ['\n', '\r\n']:
                                            splits.append(nr)
                                    record.close()
                                
                                temp_record_files = []
                                for nr, split_at in enumerate(splits):
                                    temp_record_name = record_file.split('.')[:-1]
                                    temp_record_name.append(f'temp_{nr + 1}.csv')
                                    temp_record_name = '_'.join(temp_record_name)
                                    with open(temp_record_name, 'w+') as temp_record:
                                        if nr == 0:
                                            temp_record.writelines(lines[:splits[nr + 1]])
                                        else:
                                            try:
                                                temp_record.writelines(lines[split_at:splits[nr + 1]])
                                            except:
                                                temp_record.writelines(lines[split_at:])
                                        
                                        temp_record_files.append(temp_record.name)
                                        temp_record.close()
                                
                                for i, temp_record_file in enumerate(temp_record_files):
                                    if temp_record_file.__contains__('_temp_'):
                                        map_records.insert(file_id + (i+1), temp_record_file)
                                    
                                continue
                                        
                            try:
                                record_file_df['TIME'] = pd.to_datetime(record_file_df['TIME'], dayfirst=True)     # datetime index --> chronological order
                                record_file_df = record_file_df.loc[:, :'Y_POSITION(mm)']           # slice dataframe, drop redundant columns
                                try:
                                    record_file_df.drop(record_file_df[record_file_df['SUBMAP_NUMBER'] < 0].index, inplace=True)
                                except:
                                    pass

                                record_file_df = record_file_df[record_file_df.groupby('SUBMAP_NUMBER')['SUBMAP_NUMBER'].transform('count') > 1]  # drop all rows without plan-relevant data
                            
                            except:  # unlikely event of unusable information in log-file (possible if split into parts)
                                no_exceptions = False
                                layer_exceptions.append(layer_id)
                                continue
                            
                            record_file_df = record_file_df.loc[(record_file_df['X_POSITION(mm)'] != -10000.0) & (record_file_df['Y_POSITION(mm)'] != -10000.0)]  # drop redundant rows
                            for i, submap in enumerate(record_file_df['SUBMAP_NUMBER'].drop_duplicates()):
                                submap_df = record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == submap]
                                x_pos, y_pos, x_wid, y_wid = submap_df['X_POSITION(mm)'].mean(), submap_df['Y_POSITION(mm)'].mean(), submap_df['X_WIDTH(mm)'].mean(), submap_df['Y_WIDTH(mm)'].mean()
                                if i > 0 and abs(prev_x - x_pos) < 2 and abs(prev_y - y_pos) < 2:
                                    # print('  /!\ Split spot detected', abs(record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == prev_submap, 'X_POSITION(mm)'].mean() - x_pos))
                                    record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == prev_submap, ['X_POS_IC23(mm)']] = x_pos
                                    record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == prev_submap, ['Y_POS_IC23(mm)']] = y_pos
                                    record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == prev_submap, ['X_WID_IC23(mm)']] = x_wid
                                    record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == prev_submap, ['Y_WID_IC23(mm)']] = y_wid
                                    record_file_df.drop(record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == submap].index, inplace=True)
                                else:
                                    record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == submap, ['X_POS_IC23(mm)']] = x_pos
                                    record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == submap, ['Y_POS_IC23(mm)']] = y_pos
                                    record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == submap, ['X_WID_IC23(mm)']] = x_wid
                                    record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == submap, ['Y_WID_IC23(mm)']] = y_wid

                                prev_submap = submap
                                prev_x, prev_y = x_pos, y_pos

                            record_file_df.drop_duplicates(subset=['SUBMAP_NUMBER'], keep='last', inplace=True)  # keep only last entries for each spot (most accurate)
                            if record_file_df.empty:
                                continue

                            record_file_df.index = record_file_df['TIME']   
                            record_file_df.drop(columns=['TIME'], inplace=True)
                            
                            for specif_file in record_specifs:  # draw machine parameters from *map_specif*.csv
                                if int(specif_file.split('_')[2].split('.')[0].split('_')[0]) == layer_id:
                                    with open(specif_file, 'r') as specif_file:
                                        lines = specif_file.readlines()
                                        ic_offsets = lines[3]
                                        ic_offset_x, ic_offset_y = float(ic_offsets.split(',')[3]), float(ic_offsets.split(',')[2])  # coordinate system switch x <--> y
                                        range_at_degrader = float(lines[1].split(',')[1])
                            
                            # calculate energy at isocenter from range
                            nozzle_wet = np.polyval(nozzle_wet_poly, range_at_degrader)  # [mm]
                            range_at_iso = range_at_degrader - nozzle_wet
                            layer_energy = np.exp(np.polyval(iba_gtr2_poly, np.log(range_at_iso)))  # [MeV]
                            record_file_df['LAYER_ENERGY(MeV)'] = layer_energy

                            if not np.round(layer_energy, 1) in qa_energies:
                                continue
                            
                            # coordinate system transform iba <-> raystation (x <-> y)
                            # further change x -> -y, y -> x : Lynx rotated 90° ccw when attached to nozzle holder
                            record_file_df.drop_duplicates(subset=['X_POSITION(mm)', 'Y_POSITION(mm)'], inplace=True)
                            record_file_df['Y_POSITION(mm)'] = record_file_df[['Y_POS_IC23(mm)']].apply(map_spot_pos, args=(ic_offset_x, sad_x, ictoiso_x))
                            record_file_df['X_POSITION(mm)'] = - record_file_df[['X_POS_IC23(mm)']].apply(map_spot_pos, args=(ic_offset_y, sad_y, ictoiso_y))

                            valid_spots = True
                            for (x, y) in zip(np.round(record_file_df['X_POSITION(mm)'], -1), np.round(record_file_df['Y_POSITION(mm)'], -1)):
                                if not x in qa_xy or y not in qa_xy:
                                    valid_spots = False
                                    break

                            if not valid_spots: continue
                            
                            # change x_width <-> y_width : Lynx rotated 90° ccw when attached to nozzle holder
                            record_file_df['Y_WIDTH(mm)'] = record_file_df[['Y_WID_IC23(mm)']].apply(map_spot_width, args=(sad_x, ictoiso_x))
                            record_file_df['X_WIDTH(mm)'] = record_file_df[['X_WID_IC23(mm)']].apply(map_spot_width, args=(sad_y, ictoiso_y))
                            record_file_df.drop(columns=['X_POS_IC23(mm)', 'Y_POS_IC23(mm)', 'X_WID_IC23(mm)', 'Y_WID_IC23(mm)'], inplace=True)
                            record_file_df['DIST_TO_ISO(mm)'] = np.sqrt(np.square(record_file_df['X_POSITION(mm)']) + np.square(record_file_df['Y_POSITION(mm)']))
                            record_file_df.reindex()
                            
                            if len(to_do_layers) > 0:
                                previous_record_file_df = to_do_layers[-1]
                                prev_x_last, prev_y_last = previous_record_file_df['X_POSITION(mm)'].iloc[-1], previous_record_file_df['Y_POSITION(mm)'].iloc[-1]
                                this_x_first, this_y_first = record_file_df['X_POSITION(mm)'].iloc[0], record_file_df['Y_POSITION(mm)'].iloc[0]

                                # check proximity to last spot entry of previous file (if same coordinates, spot is split over 2 files)
                                if abs(this_x_first - prev_x_last) < 1 and abs(this_y_first - prev_y_last) < 1:  
                                    print(f'''  /!\ Last spot of layer-ID {layer_id} split over 2 files, merging..''')
                                    record_file_df = record_file_df.iloc[1:]  # drop first entry of current file

                            if not record_file_df.empty: to_do_layers.append(record_file_df)  # if only spot was deleted in previous step
                    
                    if len(to_do_layers) > 0:  # can be zero, if only one spot in layer and omitted by high-weighted tuning
                        layer_df = pd.concat(to_do_layers)  # concatenate layers, assign additional columns
                        layer_df['LAYER_ID'] = layer_id
                        layer_df['TOTAL_LAYERS'] = num_layers
                        layer_df['FRACTION_ID'] = fraction_id
                        layer_df['BEAM_ID'] = beam_id
                        layer_df['GANTRY_ANGLE'] = gantry_angle
                        # layer_df['TEMPERATURE(K)'] = temperature
                        # layer_df['PRESSURE(hPa)'] = pressure
                        layer_df.drop(columns=['SUBMAP_NUMBER'], inplace=True)
                        layer_df = layer_df[~layer_df.index.duplicated(keep='first')]
                    else:
                        # print(f'  /!\ No QA record found for layer-ID {layer_id} in beam {beam_id}, continuing..')
                        continue
                    
                    del to_do_layers

                    # filter only relevant qa data
                    layer_df['X_ROUND'] = np.round(layer_df['X_POSITION(mm)'], -1)
                    layer_df['Y_ROUND'] = np.round(layer_df['Y_POSITION(mm)'], -1)  
                    layer_df['E_ROUND'] = np.round(layer_df['LAYER_ENERGY(MeV)'], 1)                
                    layer_df['DELTA_X(mm)'] = layer_df['X_POSITION(mm)'] - layer_df['X_ROUND']
                    layer_df['DELTA_Y(mm)'] = layer_df['Y_POSITION(mm)'] - layer_df['Y_ROUND']
                    layer_df['LAYER_ENERGY(MeV)'] = layer_df['E_ROUND']
                    layer_df.drop(columns=['X_ROUND', 'Y_ROUND', 'E_ROUND', 'TOTAL_LAYERS', 'LAYER_ID'], inplace=True)
                    finalized_layers.append(layer_df)

                if no_exceptions:
                    char = '#'
                else:
                    char = '_'

                if beam_no == (num_beams - 1):  # progress visualization
                    print('  ', '[' + (beam_no + 1) * char + (num_beams - beam_no - 1) * '-' + ']', end=f' Fraction {fraction_id} complete\n')
                else:
                    print('  ', '[' + (beam_no + 1) * char + (num_beams - beam_no - 1) * '-' + ']', end=f' Beam {str(beam_no + 1).zfill(2)}/{str(num_beams).zfill(2)}\r')
                                    
                no_exceptions = True

                # remove temporary files
                for file in os.listdir('.'):
                    if file.__contains__('temp'):
                        os.remove(file)

                self.qa_record_df = pd.concat(finalized_layers, sort=True)
                os.chdir(self.logfile_dir)

        # write out as .csv
        os.chdir(self.df_destination)
        self.record_df_name = 'QA_2017-2023_records_data.csv'
        
        print(f'''  ..Writing dataframe to '{self.df_destination}' as .CSV.. ''')
        while True:   
            try:
                self.qa_record_df.to_csv(self.record_df_name)
                break
            except PermissionError:
                input('  Permission denied, close target file and press ENTER.. ')

        self.qa_record_df = pd.read_csv('QA_2017-2023_records_data.csv', dtype={'BEAM_ID':str, 'FRACTION_ID':str})
        print('Complete')
    

    '''
    Shared searching function invoked by other methods, identify treatment plan matching given log-file data
    Requires:   Written .csv dataframe in self.dataframe_destination
    Arguments:  fraction_id/beam_id - Desired treatment fraction/beam present in dataframe
                verbose - if True, print output and warning, else supress them
    Operation:  Perform os.walk() combined with filter criteria to find plan containing desired fraction beam
    Returns:    plan_dcm - DICOM filename
                beam_ds - beam dataset, read-in with pydicom
    '''
    def dicom_finder(self, fraction_id, beam_id, verbose=False):
        fraction_df = self.patient_record_df.loc[self.patient_record_df['FRACTION_ID'] == fraction_id]  # slice patient dataframe
        
        # beam naming by planners might be inconsistent (BeamName vs. BeamDescription tag in dicom)
        for bid in fraction_df['BEAM_ID'].drop_duplicates():
            if str(bid).__contains__(str(beam_id)) or str(beam_id).__contains__(str(bid)):
                beam_df = fraction_df.loc[fraction_df['BEAM_ID'] == bid]
                break

        gtr_angle = beam_df['GANTRY_ANGLE'].iloc[0]
        n_layers = beam_df['TOTAL_LAYERS'].iloc[0]
        log_energies = np.array([np.round(e, 1) for e in sorted(beam_df['LAYER_ENERGY(MeV)'].drop_duplicates().to_list())])
        found = False
        
        if verbose: print(f'Looking for beam-ID {beam_id} with GTR{gtr_angle}, n_layers = {n_layers}, E = [{np.round(min(log_energies), 1)} .. {np.round(max(log_energies), 1)}]')

        # auto-location of plan DICOM    
        delivered = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\Logfiles\DeliveredPlans'
        plan_dir = os.path.join(delivered, self.patient_id)
        if not os.path.isdir(delivered):
            plan_dir = r'/home/luke/Scripts/Logfile_Extraction/1676348/DeliveredPlans'        
        
        if verbose: print('''... will only accept RayStation plan exports ('RP.*.dcm') ...''')
        for file in os.listdir(plan_dir):
            if file.__contains__('RP') and file.endswith('.dcm') and not file.__contains__('log'):
                ds = pydicom.read_file(os.path.join(plan_dir, file))
                try:
                    beam_seq = ds.IonBeamSequence
                except AttributeError:
                    if verbose: print(f'  /x\ DICOM {file} is not ion RTPLAN')
                    continue

                for beam in beam_seq:
                    plan_energies = np.array(pd.Series(sorted([layer.NominalBeamEnergy for layer in beam.IonControlPointSequence])).drop_duplicates().to_list())
                    
                    # check gantry angle, total layers, beam energies. Do not check names, they are non-standardized
                    if float(beam.IonControlPointSequence[0].GantryAngle) == gtr_angle and len(beam.IonControlPointSequence) == n_layers * 2:  # layer sequence in DICOM has double entries, BUT first/last layer might be missing due to single-spot layer tuning
                        try:
                            max_energy_diff = np.max(np.abs(plan_energies - log_energies))
                            if max_energy_diff < 0.1:  # tolerance on energy diff to plan, logs normally below dE = 0.1MeV
                                plan_dcm = os.path.join(plan_dir, file)
                                beam_ds = beam
                                found = True
                                break
                    
                        except:
                            for plan_e in plan_energies:
                                if plan_e not in log_energies:
                                    if plan_e == max(plan_energies):
                                        plan_dcm = os.path.join(plan_dir, file)
                                        beam_ds = beam
                                        found = True
                                        if verbose: 
                                            n_spots = beam.IonControlPointSequence[0].NumberOfScanSpotPositions
                                            if n_spots <= 1: print(f'  /!\ First layer with energy {plan_e} contains {n_spots} spot(s) - omitted by tuning?')
                                            else:
                                                print(f'EXIT CODE 1 - {plan_e}, {n_spots} > 1') 
                                                return None

                                    elif plan_e == min(plan_energies):
                                        plan_dcm = os.path.join(plan_dir, file)
                                        beam_ds = beam
                                        found = True
                                        if verbose: 
                                            n_spots = beam.IonControlPointSequence[-2].NumberOfScanSpotPositions
                                            if n_spots <= 1: print(f'  /!\ Last layer with energy {plan_e} contains {n_spots} spot(s) - omitted by tuning?')
                                            else: 
                                                print(f'EXIT CODE 2 - {plan_e}, {n_spots} > 1')
                                                return None
                                        
                                    else:
                                        if verbose: print(f'  /x\ Log-file is missing energy {plan_e}, skipping beam {beam_id}..')
                                        return None

                            if found: break

                        if verbose:
                            print('\n\t\tLOG\tPLAN')
                            print(f'Gantry\t\t{gtr_angle}\t{beam.IonControlPointSequence[0].GantryAngle}')
                            print(f'Layers\t\t{n_layers}\t{len(beam.IonControlPointSequence) / 2}')
                            print(f'Energy max.\t{max(log_energies)}\t{max(plan_energies)}')
                            print(f'Energy min.\t{min(log_energies)}\t{min(plan_energies)}')

        # except:   
        #     print('Pre-defined directory search failed')
        #     return None
        #     for path, dirnames, filenames in os.walk(os.path.join(self.logfile_dir, '..')):
        #         for fname in filenames:
        #             if fname.endswith('.dcm') and not fname.__contains__('log'):
        #                 ds = pydicom.read_file(os.path.join(path, fname))
        #                 for beam in ds.IonBeamSequence:
        #                     plan_energies = np.array(pd.Series(sorted([layer.NominalBeamEnergy for layer in beam.IonControlPointSequence])).drop_duplicates().to_list())
        #                     log_energies = np.array(sorted(beam_df['LAYER_ENERGY(MeV)'].drop_duplicates().to_list()))
                            
        #                     # check gantry angle, total layers, beam energies. Do not check names, they are non-standardized
        #                     if float(beam.IonControlPointSequence[0].GantryAngle) == gtr_angle and len(beam.IonControlPointSequence) == n_layers * 2:  # layer sequence in DICOM has double entries
        #                         # print(beam.IonControlPointSequence[0].GantryAngle, gtr_angle, len(beam.IonControlPointSequence), n_layers * 2)
        #                         max_energy_diff = np.max(np.abs(plan_energies - log_energies))
        #                         if max_energy_diff < 0.1:  # tolerance on energy diff to plan, logs normally below dE = 0.1MeV
        #                             plan_dcm = os.path.join(path, fname)
        #                             beam_ds = beam
        #                             if verbose:
        #                                 if not found:
        #                                     print('    Found matching DICOM -', fname)
        #                                 else:
        #                                     print('    /!\ Further DICOM match -', fname)
        #                             found = True
                            
        while not found:  # fallback: open dicom file manually if failed
            if not verbose: return None

            root = Tk()
            print(f'    /!\ Auto-location failed for beam-ID {beam_id} in fraction-ID {fraction_id}, select plan DICOM..')
            plan_dcm = filedialog.askopenfilename(initialdir=os.path.join(self.logfile_dir, '..'))
            root.destroy()
            if plan_dcm == '':
                print('    /x\ Process cancelled by user')
                return None
                
            ds = pydicom.read_file(plan_dcm)
            for beam in ds.IonBeamSequence:  # check only gtr-angle and number of layers, omit energy check
                if float(beam.IonControlPointSequence[0].GantryAngle) == gtr_angle:
                    beam_ds = beam
                    found = True
        
        if verbose: print('EXIT CODE 0')
        return plan_dcm, beam_ds


    '''
    Shared function invoked by other methods, changes log-file spot sequence order to match plan
    Requires:   Written .csv dataframe in self.dataframe_destination, existent plan DICOM within parent dir
    Arguments:  fraction_id/beam_id - Desired treatment fraction/beam present in dataframe
                mode:
                    record - sort irradiated spots from record file to plan, neglect tuning spots
                    tuning - identify spots omitted by ScanAlgo in case of high-weighted tuning (important for dose recon)
    Operation:  Map spots in each log layer to plan layer sequence by minimum euclidian distance (nearest neighbour search)
    Returns:    sorting_dict - nested dictionary working on spot_id in respective layer_id, i.e. sorting_dict[layer_id][spot_id] = plan_id
    '''
    def spot_sorter(self, fraction_id, beam_id, mode='record'):
        to_be_sorted = self.patient_record_df.loc[(self.patient_record_df['FRACTION_ID'] == fraction_id) & (self.patient_record_df['BEAM_ID'] == beam_id)]  # slice dataframe
        n_layers = to_be_sorted['TOTAL_LAYERS'].iloc[0]  # get total beam layers
        sorting_dict = {lid:{} for lid in range(n_layers)}  # initialize empty dictionary

        # use dicom_finder() to draw plan and beam
        plan_dcm, beam_ds = self.dicom_finder(fraction_id, beam_id, verbose=True)
        
        # start sorting procedure
        print(f'  Sorting spots for beam-ID {beam_id}..')
        for layer_id in to_be_sorted['LAYER_ID'].drop_duplicates():
            try:
                plan_layer = beam_ds.IonControlPointSequence[layer_id * 2]
            except IndexError:  # warning, possibly wrong plan
                print(f'''  /!\ Layer-ID mismatch, skipping layer #{layer_id + 1} in beam {beam_id}, fraction {fraction_id}''')
                continue

            plan_spotmap = plan_layer.ScanSpotPositionMap
            plan_x, plan_y = [], []

            # alternating read-out of x- and y-values in plan spotmap
            for i, spot in enumerate(plan_spotmap):
                if i % 2 == 0:
                    plan_x.append(spot)
                else:
                    plan_y.append(spot)
            plan_xy = [tup for tup in zip(plan_x, plan_y)]

            log_layer = to_be_sorted.loc[(to_be_sorted['LAYER_ID'] == layer_id) & (to_be_sorted['FRACTION_ID'] == fraction_id)]
            log_xy = [tup for tup in zip(log_layer['X_POSITION(mm)'].to_list(), log_layer['Y_POSITION(mm)'].to_list())]
            log_xy_sorted = [(np.nan, np.nan) for _ in range(max(len(log_xy), len(plan_xy)))]  # initialize empty sorted list with sufficient length

            # match (x,y)-positions to plan
            for i, log_spot in enumerate(log_xy):
                shifts = [np.array(plan_spot) - np.array(log_spot) for plan_spot in plan_xy]
                dists = [abs((shift).dot(shift)) for shift in shifts]
                index = dists.index(min(dists))
                log_xy_sorted[index] = log_xy[i]

            # drop NaN's from sorted list, these represent planned spots omitted by ScanAlgo in case of high-weighted tuning spot
            dropped = 0
            for drop, xy in enumerate(log_xy_sorted):
                if xy == (np.nan, np.nan):
                    log_xy_sorted.remove(xy)
                    drop_id = drop
                    dropped += 1
            
            # only the tuning spot might be dropped, if more than one, show corrupted spot map to user and return None
            if dropped > 1:
                print(f'Dropped {dropped} spots | fx-ID {fraction_id} | beam-ID {beam_id} | layer-ID {layer_id}')
                plt.plot(*zip(*plan_xy), marker='x', ls='-', color='tab:blue', label='plan')
                # plt.plot(*plan_xy[drop_id], marker='x', ls='-', color='black', label='omitted')
                plt.plot(*zip(*log_xy), marker='o', ls='--', color='tab:grey', label='log')
                plt.plot(*zip(*log_xy_sorted), marker='o', ls='-', color='black', label='sorted')
                plt.legend()
                plt.show()
                print(f'  /!\ Log beam {beam_id} spots do not match plan beam {beam_ds.BeamName}')
                return None

            # if sorting_dict is used for tuning spot sorting, use omitted spot information
            if mode == 'tuning':
                try:
                    sorting_dict[layer_id] = drop_id
                except:
                    sorting_dict[layer_id] = None
            else:
                sorting_dict[layer_id] = {log_xy.index(log_xy[i]) : log_xy.index(log_xy_sorted[i]) for i in range(len(log_xy))}  # do not iterate over elements directly, index() fails in this case

        return sorting_dict       
        

    '''
    Plot all spotmaps (layers) of last fraction of selected beam
    '''
    def plot_beam_layers(self):
        # For all layers and last fraction of selected beam:
        # derive scatter plot of (x,y)-positions from logfile, 
                                    # compare vs. planned positions extracted from RP*.dcm
        beam_list = self.patient_record_df['BEAM_ID'].drop_duplicates()
        indices = beam_list.index.to_list()

        print('\nSelect beam key for layer-wise spot plotting:\n')
        print('Key\tBeam-ID\t\tGantry Angle\tLayers\tFrom Fraction\n')  # choose beam
        for choice, (index, beam) in enumerate(zip(indices, beam_list)):
            num_layers = int(self.patient_record_df['TOTAL_LAYERS'][index].mean())
            fraction = self.patient_record_df['FRACTION_ID'][index]
            try:
                gantry_angle = self.patient_record_df['GANTRY_ANGLE'][index].iloc[0]
            except:
                gantry_angle = self.patient_record_df['GANTRY_ANGLE'][index]
            print(f'({choice + 1})\t{beam}\t\t{gantry_angle}\t\t{num_layers}\t{fraction}')

        while True:
            try:
                key = int(input('\n Select beam key: '))
                if key > len(beam_list) or key <= 0:
                    print('Key out of bounds, select another..')
                    continue
                else:
                    break
            except:
                print('Invalid input, try again..')

        beam_id = str(beam_list[key - 1])
        
        scope_record_df = self.patient_record_df.loc[self.patient_record_df['BEAM_ID'] == beam_id]  # slice currently needed portion from patient df
        scope_tuning_df = self.patient_tuning_df.loc[self.patient_tuning_df['BEAM_ID'] == beam_id]
        print('Selected beam is:', beam_id)

        fx_list = scope_record_df['FRACTION_ID'].drop_duplicates()
        indices = fx_list.index.to_list()

        print('\nSelect fraction to plot:\n')
        print('Key\tFraction\n')  # choose beam
        for choice, (index, fx) in enumerate(zip(indices, fx_list)):
            print(f'({choice + 1})\t{fx}')

        while True:
            try:
                key = int(input('\n Select fraction key: '))
                if key > len(fx_list) or key <= 0:
                    print('Key out of bounds, select another..')
                    continue
                else:
                    break
            except:
                print('Invalid input, try again..')

        fx_id = str(fx_list[key - 1])
        print('Selected fraction is:', fx_id)

        scope_record_df = scope_record_df.loc[scope_record_df['FRACTION_ID'] == fx_id]
        scope_tuning_df = scope_tuning_df.loc[scope_tuning_df['FRACTION_ID'] == fx_id]

        print('Tuning violations:')
        # for i, tune_spot in enumerate(scope_tuning_df):
        #     if abs(scope_tuning_df['X_POSITION(mm)'].iloc[i]) > 150 or abs(scope_tuning_df['Y_POSITION(mm)'].iloc[i]) > 193:
        #         print('  ', '150mm < x =', scope_tuning_df['X_POSITION(mm)'].iloc[i])
        #         print('  ', '193mm < y =', scope_tuning_df['Y_POSITION(mm)'].iloc[i])
        
        # print('Record violations:')
        # for i, rec_spot in enumerate(scope_record_df):
        #     if abs(scope_record_df['X_POSITION(mm)'].iloc[i]) > 150 or abs(scope_record_df['Y_POSITION(mm)'].iloc[i]) > 193:
        #         print('  ', '150mm < x =', scope_record_df['X_POSITION(mm)'].iloc[i])
        #         print('  ', '193mm < y =', scope_record_df['Y_POSITION(mm)'].iloc[i])

        print('\nGenerating layer plot..')
        plan_dcm, beam_ds = self.dicom_finder(fraction_id=fx_id, beam_id=beam_id, verbose=True)
        beam_sorting_dict = self.spot_sorter(fraction_id=fx_id, beam_id=beam_id)

        fig, axs = plt.subplots(6, 8, sharex=True, sharey=True, figsize=(24, 24 * 6/8), dpi=150)  # initiate matrix-like layer plot
        ax0 = fig.add_subplot(111, frameon=False)
        fig.subplots_adjust(hspace=0.0, wspace=0.0)
        axs = axs.ravel()

        for layer_id in scope_record_df['LAYER_ID'].drop_duplicates():
            layer_spot_df = scope_record_df.loc[scope_record_df['LAYER_ID'] == layer_id]
            layer_tuning_df = scope_tuning_df.loc[scope_tuning_df['LAYER_ID'] == layer_id]
            dcm_layer = beam_ds.IonControlPointSequence[layer_id * 2]
            plan_spotmap = dcm_layer.ScanSpotPositionMap
            plan_mu = dcm_layer.ScanSpotMetersetWeights
            if type(plan_mu) == float:
                plan_mu = [plan_mu]
            plan_x_positions, plan_y_positions = [], []
            for i, coord in enumerate(plan_spotmap):
                if i % 2 == 0:
                    plan_x_positions.append(coord)
                else:
                    plan_y_positions.append(coord)
            
            spot_points_log = [tuple for tuple in zip(layer_spot_df['X_POSITION(mm)'].to_list(), layer_spot_df['Y_POSITION(mm)'].to_list())]
            spot_mu_log = layer_spot_df['MU'].to_list()
            tuning_points_log = [tuple for tuple in zip(layer_tuning_df['X_POSITION(mm)'].to_list(), layer_tuning_df['Y_POSITION(mm)'].to_list())]
            spot_points_sorted = [spot_points_log[beam_sorting_dict[layer_id][i]] for i in range(len(spot_points_log))]
            spot_mu_sorted = [spot_mu_log[beam_sorting_dict[layer_id][i]] for i in range(len(spot_mu_log))]
            spot_mu_delta = [spot_mu_sorted[i] - plan_mu[i] for i in range(len(spot_mu_sorted))]
            
            if self.patient_id == '1669130' and beam_id == '3' and layer_id == 13:
                sns.set(style='ticks', context='paper', font_scale=1.5)
                sns.set_style({"xtick.direction": "in","ytick.direction": "in", "ax.axesbelow":True, "xtick.top":True, "ytick.right":True})
                fig2, ax2 = plt.subplots(figsize=(4, 4))
                ax2.plot(plan_x_positions, plan_y_positions, marker='x', linestyle='-', ms=7, label='Planned')
                ax2.plot(*zip(*spot_points_sorted), marker='o', markerfacecolor='None', linestyle='-', ms=7, color='black', label='Log file')
                ax2.plot(*zip(*tuning_points_log), marker='o', linestyle='None', color='limegreen', ms=7, markeredgecolor='black', alpha=0.8, label='Tuning spot(s)')
                ax2.set_xlabel('$x$ [mm]')
                ax2.set_ylabel('$y$ [mm]')
                ax2.set_xlim(-50, 0)
                ax2.set_ylim(-30, 20)
                ax2.legend(loc='upper left')
                ax2.grid(axis='both')
                ax2.text(0.05, 0.05, 'B', color='black', fontweight='bold', fontsize=24., transform=ax2.transAxes, ha='left', va='bottom', bbox=dict(facecolor='white', edgecolor='black', pad=10.0))
                # ax2.grid()

                fig2.tight_layout()
                fig2.savefig(f'{output_dir}/{self.patient_id}_{beam_id}_{fx_id}_spotmap_show.png', dpi=300)
                # plt.show()
               
                # return None
            if self.patient_id == '1635796' and beam_id == '1' and layer_id == 24:
                sns.set(style='ticks', context='paper', font_scale=1.5)
                sns.set_style({"xtick.direction": "in","ytick.direction": "in", "ax.axesbelow":True, "xtick.top":True, "ytick.right":True})
                fig2, ax2 = plt.subplots(figsize=(4, 4))
                ax2.plot(plan_x_positions, plan_y_positions, marker='x', linestyle='-', ms=7, label='Planned')
                ax2.plot(*zip(*spot_points_sorted), marker='o', markerfacecolor='None', linestyle='-', ms=7, color='black', label='Log file')
                ax2.plot(*zip(*tuning_points_log), marker='o', linestyle='None', color='limegreen', ms=7, markeredgecolor='black', alpha=0.8)
                ax2.set_xlabel('$x$ [mm]')
                ax2.set_ylabel('$y$ [mm]')
                ax2.set_xlim(-100, -50)
                ax2.set_ylim(45, 95)
                ax2.legend(loc='upper left')
                ax2.grid(axis='both')
                # ax2.text(0.05, 0.05, 'B', color='black', fontweight='bold', fontsize=24., transform=ax2.transAxes, ha='left', va='bottom', bbox=dict(facecolor='white', edgecolor='black', pad=10.0))
                ax2.text(0.05, 0.05, 'A', color='black', fontweight='bold', fontsize=24., transform=ax2.transAxes, ha='left', va='bottom', bbox=dict(facecolor='white', edgecolor='black', pad=10.0))

                
                # ax2.grid()

                fig2.tight_layout()
                fig2.savefig(f'{output_dir}/{self.patient_id}_{beam_id}_{fx_id}_spotmap_show.png', dpi=300)

            axs[layer_id].plot(plan_x_positions, plan_y_positions, marker='x', linestyle='-', markersize=2.0, markeredgewidth=0.2, linewidth=0.2, label='Planned')
            # MU colour coding still bugged, needs deltaframe for functionality
            axs[layer_id].scatter(*zip(*spot_points_sorted), c=spot_mu_delta, cmap='bwr', vmin=-0.1, vmax=0.1, edgecolors='black', linewidths=0.2, s=7, label='Log-file sorted', zorder=5)
            axs[layer_id].plot(*zip(*spot_points_sorted), marker='None', linestyle='-', lw=0.2, color='black')
            # axs[layer_id].plot(*zip(*spot_points_log), marker='o', markerfacecolor='None', linestyle='--', color='black', markersize=2.0, markeredgewidth=0.2, linewidth=0.2, label='Log-file original')
            axs[layer_id].plot(*zip(*tuning_points_log), marker='o', markerfacecolor='None', linestyle='None', markersize=2.0, markeredgewidth=0.2, color='limegreen', label='Tuning spot(s)')
            axs[layer_id].annotate(f'Layer #{str(layer_id + 1).zfill(2)} | $\Delta$ = {abs(len(plan_x_positions) - len(spot_points_log))}', xy=(1.0, 1.0), xycoords='axes points', fontsize=8)
            
            if len(tuning_points_log) > 1:
                axs[layer_id].annotate(f'Layer #{str(layer_id + 1).zfill(2)} | $\Delta$ = {abs(len(plan_x_positions) - len(spot_points_log))} | $n_t$ > 1', xy=(1.0, 1.0), xycoords='axes points', fontsize=8) 
            axs[layer_id].legend(loc='upper right', fontsize=8)

        plt.suptitle(f'Spot Positions for Patient-ID {self.patient_id}, Beam-ID: {beam_id}', fontweight='bold', y=0.9)
        ax0.set_xlabel('X [mm]', fontweight='bold', labelpad=10)
        ax0.set_ylabel('Y [mm]', fontweight='bold', labelpad=10)
        plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)

        print(f'  Saving figure to {output_dir}..')
        while True:
            try:
                fig.savefig(f'{output_dir}/{self.patient_id}_{beam_id}_spots.pdf')
                break
            except PermissionError:
                input('  /!\ Permission denied, close target file and press ENTER.. ')
        plt.close()
        plt.clf()


    def plot_spot_statistics(self):     # For one beam over all fractions:
                                        # derive histograms and MAE-evolution by comparing log vs. plan spot positions and MU's

        beam_list = [str(i) for i in self.patient_record_df['BEAM_ID'].drop_duplicates()]

        print('\nTrying to auto-locate patient plan dicoms..')  # read RT plan dicom via filedialog
        patient_dicoms = []
        for path, dirnames, filenames in os.walk(os.path.join(self.logfile_dir, '..')):
            for fname in filenames:
                if fname.__contains__('RP') and fname.endswith('.dcm'):
                    ds = pydicom.read_file(os.path.join(path, fname))
                    for i, beam in enumerate(ds.IonBeamSequence):
                        if (beam.BeamName in beam_list or beam.BeamDescription in beam_list) and float(beam.IonControlPointSequence[0].GantryAngle) in self.patient_record_df['GANTRY_ANGLE'].drop_duplicates().to_list():
                            patient_dicoms.append(os.path.join(path, fname))

        print('Collecting data..')
        # fig, axs = plt.subplots(3, len(self.fraction_list), sharey='row', figsize=(100, 15), dpi=100)
        fig, axs = plt.subplots(3, 2, sharey='row', figsize=(10, 10), dpi=100)
        ax0 = fig.add_subplot(111, frameon=False)
        ax0.set_xticks([])
        ax0.set_yticks([])
        axs[0, 0].set_ylabel('MAE [mm]')
        axs[1, 0].set_ylabel('MAE [mm]')
        axs[1, 0].set_xticks([])
        axs[2, 0].set_ylabel('Counts')
        fig.subplots_adjust(wspace=0.0)
        axs.ravel()
        for f, fraction_id in enumerate(self.fraction_list):
            if fraction_id != '20200731':
                continue
            f = 0
            print('!')

            fraction_x_maes, fraction_y_maes, layer_axis = [], [], []
            total_x_diffs, total_y_diffs = [], []
            beams_in_frac = self.patient_record_df.loc[self.patient_record_df['FRACTION_ID'] == fraction_id]['BEAM_ID'].drop_duplicates()
            for b, beam_id in enumerate(beams_in_frac):
                scope_record_df = self.patient_record_df.loc[(self.patient_record_df['BEAM_ID'] == beam_id) & (self.patient_record_df['FRACTION_ID'] == fraction_id)]  # slice currently needed portion from patient df
                scope_tuning_df = self.patient_tuning_df.loc[(self.patient_tuning_df['BEAM_ID'] == beam_id) & (self.patient_tuning_df['FRACTION_ID'] == fraction_id)]

                for dcm in patient_dicoms:
                    ds = pydicom.read_file(dcm)
                    for i, beam in enumerate(ds.IonBeamSequence):
                        if (beam.BeamName == beam_id or beam.BeamDescription == beam_id) and float(beam.IonControlPointSequence[0].GantryAngle) == scope_record_df.loc[scope_record_df['BEAM_ID'] == beam_id]['GANTRY_ANGLE'].mean():
                            dcm_beam = ds.IonBeamSequence[i]

                layer_x_maes, layer_y_maes = [], []
                for layer_id in scope_record_df['LAYER_ID'].drop_duplicates():
                    x_positions, y_positions = [], []
                    x_tunings, y_tunings, = [], []
                    dcm_layer = dcm_beam.IonControlPointSequence[layer_id * 2]
                    plan_spotmap = dcm_layer.ScanSpotPositionMap
                    plan_x_positions, plan_y_positions = [], []
                    layer_x_diffs, layer_y_diffs = [], []
                    for i, coord in enumerate(plan_spotmap):
                        if i % 2 == 0:
                            plan_x_positions.append(coord)
                        else:
                            plan_y_positions.append(coord)
                    
                    for spot_id in scope_record_df.loc[scope_record_df['LAYER_ID'] == layer_id]['SPOT_ID'].drop_duplicates():
                        x_positions.append(scope_record_df.loc[(scope_record_df['LAYER_ID'] == layer_id) & (scope_record_df['SPOT_ID'] == spot_id)]['X_POSITION(mm)'].mean())
                        y_positions.append(scope_record_df.loc[(scope_record_df['LAYER_ID'] == layer_id) & (scope_record_df['SPOT_ID'] == spot_id)]['Y_POSITION(mm)'].mean())
                    for tuning_id in scope_tuning_df.loc[scope_tuning_df['LAYER_ID'] == layer_id]['SPOT_ID'].drop_duplicates():   
                        x_tunings.append(scope_tuning_df.loc[(scope_tuning_df['LAYER_ID'] == layer_id) & (scope_tuning_df['SPOT_ID'] == tuning_id)]['X_POSITION(mm)'].mean())
                        y_tunings.append(scope_tuning_df.loc[(scope_tuning_df['LAYER_ID'] == layer_id) & (scope_tuning_df['SPOT_ID'] == tuning_id)]['Y_POSITION(mm)'].mean())

                    spot_points_log = [tuple for tuple in zip(x_positions, y_positions)]
                    spot_points_plan = [tuple for tuple in zip(plan_x_positions, plan_y_positions)]
                    for xy_log in spot_points_log:
                        x_differences, y_differences = [], []
                        distances = []
                        for xy_plan in spot_points_plan:  # calculate position difference and euclidian distance for every spot pair, find minimum to reindex (return sorted spots)
                            diff = np.array(xy_plan) - np.array(xy_log)
                            sqdist = abs((diff).dot(diff))
                            distances.append(sqdist)
                            x_differences.append(diff[0]), y_differences.append(diff[1])

                        min_dist = min(distances)
                        index = distances.index(min_dist)
                        total_x_diffs.append(x_differences[index]), total_y_diffs.append(y_differences[index])
                        layer_x_diffs.append(abs(x_differences[index])), layer_y_diffs.append(abs(y_differences[index]))
                    
                    layer_axis.append(str(layer_id))
                    layer_x_mae, layer_y_mae = np.mean(layer_x_diffs), np.mean(layer_y_diffs)
                    layer_x_maes.append(layer_x_mae), layer_y_maes.append(layer_y_mae)
                
                fraction_x_maes.append(layer_x_maes), fraction_y_maes.append(layer_y_maes)

            fraction_x_beam_maes = [np.mean(xerr_array) for xerr_array in fraction_x_maes]
            fraction_y_beam_maes = [np.mean(yerr_array) for yerr_array in fraction_y_maes]

            print(f'Generating MAE evolution plot {str(f + 1).zfill(2)}/{len(self.fraction_list)}', end='\r')
            colors = ['black', 'tab:blue', 'tab:orange', 'tab:red']
            axs[1, f].set_xticklabels([])
            for b, beam_id in enumerate(beams_in_frac):
                axs[0, f].plot(fraction_x_maes[b], label=f'{beam_id} - $x$', linewidth=1.0, linestyle='-', color=colors[b])
                axs[0, f].plot(fraction_y_maes[b], label=f'{beam_id} - $y$', linewidth=1.0, linestyle='--', color=colors[b])
                bp1 = axs[1, f].boxplot(fraction_x_maes[b], positions=[b + 1], widths=0.7)
                plt.setp(bp1['medians'], color=colors[b], linestyle='-')
                axs[1, f].axvline(len(beams_in_frac) + 1, color='black', lw=1.0, ls='--')
                bp2 = axs[1, f].boxplot(fraction_y_maes[b], positions=[b + 2 + len(beams_in_frac)], widths=0.7)
                plt.setp(bp2['medians'], color=colors[b], linestyle='--')
            axs[0, f].grid(axis='y')
            axs[0, f].set_xlabel('Layer-ID')
            axs[0, f].legend(loc='upper right')
            axs[0, f].set_title(f'{fraction_id}')
            axs[1, f].set_xlabel('Beams')
            axs[1, f].grid(axis='y')
            axs[1, f].set_xticks([])
            bin_width = 0.05  # mm
            x_bins = np.arange(min(total_x_diffs), max(total_x_diffs) + bin_width, bin_width)
            y_bins = np.arange(min(total_y_diffs), max(total_y_diffs) + bin_width, bin_width)
            axs[2, f].hist(total_x_diffs, bins=x_bins, color='tab:green', alpha=0.5, edgecolor='black', linewidth=1.0, label='$\Delta x$')
            axs[2, f].hist(total_y_diffs, bins=y_bins, color='tab:purple', alpha=0.5, edgecolor='black', linewidth=1.0, label='$\Delta y$')
            axs[2, f].axvline(np.mean(total_x_diffs), color='tab:green', alpha=0.5, lw=1.0, ls='--', label='mean')
            axs[2, f].axvline(np.mean(total_y_diffs), color='tab:purple', alpha=0.5, lw=1.0, ls='--')
            axs[2, f].set_xlim(-1.75, 1.75)
            axs[2, f].set_xlabel('$\Delta$(Plan, Log) [mm]')
            axs[2, f].annotate(f'$\mu_x =$ {np.mean(total_x_diffs):.3f}\n$\mu_y =$ {np.mean(total_y_diffs):.3f}', xy=(10.0, 120.0), xycoords='axes points')
            axs[2, f].legend(loc='upper right')
        
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        while True:
            try:
                plt.savefig(f'{output_dir}/{self.patient_id}_position_MAE_evo.pdf')
                break
            except PermissionError:
                input('  Permission denied, close target file and press ENTER.. ')
        plt.clf()


    def prepare_deltaframe(self):
        for df_file in sorted(os.listdir(self.df_destination)):
            if df_file.__contains__('delta') and df_file.__contains__(str(self.patient_id)) and df_file.endswith('.csv'):
                re_init = input(f'''Found existing patient deltaframe '{df_file}', re-initialize [y/n]? ''')
                if re_init == 'y':
                    break
                else:
                    self.patient_delta_df = pd.read_csv(os.path.join(self.df_destination, df_file), index_col='UNIQUE_INDEX', dtype={'BEAM_ID':str, 'FRACTION_ID':str})
                    return None

        # initialize deltaframe
        print(f'Initializing deltaframe for patient-ID {self.patient_id}..')
        self.patient_delta_df = pd.DataFrame(columns=self.patient_record_df.columns, index=self.patient_record_df.index)
        self.patient_delta_df.rename(columns={  'X_POSITION(mm)':'DELTA_X(mm)',
                                                'Y_POSITION(mm)':'DELTA_Y(mm)',
                                                'MU':'DELTA_MU',
                                                'LAYER_ENERGY(MeV)':'DELTA_E(MeV)'
                                             }, inplace=True)
        self.patient_delta_df['UNIQUE_INDEX'] = np.nan
        self.patient_delta_df.index = self.patient_delta_df['UNIQUE_INDEX']
        self.patient_delta_df = self.patient_delta_df.iloc[0:0]  # empty deltaframe
    	
        to_concat, unique_index = [], 0
        for f, fx_id in enumerate(self.patient_record_df['FRACTION_ID'].drop_duplicates()):
            fraction_df = self.patient_record_df.loc[self.patient_record_df['FRACTION_ID'] == fx_id]
            for beam_id in fraction_df['BEAM_ID'].drop_duplicates():
                has_plan = False
                try:    
                    plan_dcm, beam_ds = self.dicom_finder(fraction_id=fx_id, beam_id=beam_id, verbose=True)
                    has_plan = True
                except:
                    print(f'''  /!\ No plan dicom found for beam-ID {beam_id} in fraction-ID {fx_id}, setting NaN's..''')
                    
                # for bid in beam_df['BEAM_ID'].drop_duplicates():
                #     if str(bid).__contains__(str(beam_id)) or str(beam_id).__contains__(str(bid)):
                #         beam_df = beam_df.loc[beam_df['BEAM_ID'] == bid]
                #         break
                beam_df = fraction_df.loc[fraction_df['BEAM_ID'] == beam_id]
                num_layers = beam_df['TOTAL_LAYERS'].iloc[0]
                gtr_angle = float(beam_df['GANTRY_ANGLE'].iloc[0])

                for layer_id in beam_df['LAYER_ID'].drop_duplicates():
                    log_layer = beam_df.loc[beam_df['LAYER_ID'] == layer_id]
                    layer_missing = False
                    if has_plan:
                        try:
                            plan_layer = beam_ds.IonControlPointSequence[layer_id * 2]
                        except IndexError:
                            print(f'''  /!\ Layer-ID mismatch, skipping layer #{layer_id + 1} in beam {beam_id}, fraction {fx_id}''')
                            layer_missing = True
                        
                        if not layer_missing:
                            plan_spotmap = plan_layer.ScanSpotPositionMap
                            plan_x, plan_y = [], []
                            for i, spot in enumerate(plan_spotmap):
                                if i % 2 == 0:
                                    plan_x.append(spot)
                                else:
                                    plan_y.append(spot)
                            plan_xy = [tup for tup in zip(plan_x, plan_y)]
                            plan_e = plan_layer.NominalBeamEnergy
                            try:
                                plan_mu = [mu for mu in plan_layer.ScanSpotMetersetWeights]
                            except TypeError:
                                plan_mu = [plan_layer.ScanSpotMetersetWeights]

                            log_xy = [tup for tup in zip(log_layer['X_POSITION(mm)'], log_layer['Y_POSITION(mm)'])]
                            log_mu, log_e = log_layer['MU'], log_layer['LAYER_ENERGY(MeV)'].iloc[0]
                            log_xy_sorted = [(np.nan, np.nan) for _ in range(max(len(log_xy), len(plan_xy)))]
                            log_mu_sorted = [np.nan for _ in log_xy_sorted]

                            delta_x, delta_y, delta_mu = [], [], []
                            delta_e = plan_e - log_e.mean()

                            # match (x,y)-positions to plan, transform MU list equally
                            for i, log_spot in enumerate(log_xy):
                                shifts = [np.array(log_spot) - np.array(plan_spot) for plan_spot in plan_xy]
                                dists = [np.abs((shift).dot(shift)) for shift in shifts]
                                index = dists.index(min(dists))
                                dx, dy = shifts[index]
                                delta_x.append(dx), delta_y.append(dy)
                                delta_mu.append(log_mu[i] - plan_mu[index])
                                log_xy_sorted[index] = log_xy[i]
                                log_mu_sorted[index] = log_mu[i]

                            dropped = 0
                            for i, xy in enumerate(log_xy_sorted):
                                if xy == (np.nan, np.nan) :
                                    log_xy_sorted.remove(xy)
                                    log_mu_sorted.remove(log_mu_sorted[i])
                                    dropped += 1
                            
                            if dropped > 1:
                                print(f'Dropped {dropped} spots | fx-ID {fx_id} | beam-ID {beam_id} | layer-ID {layer_id}')
                                # plt.plot(*zip(*plan_xy), marker='x', ls='-', color='tab:blue', label='plan')
                                # plt.plot(*zip(*log_xy), marker='o', ls='--', color='tab:grey', label='log')
                                # plt.plot(*zip(*log_xy_sorted), marker='o', ls='-', color='black', label='sorted')
                                # plt.legend()
                                # plt.show()
                                print(f'  /!\ Log beam {beam_id} spots do not match plan beam {beam_ds.BeamName}')

                    # generate new dataframe 
                    fx_delta_df = pd.DataFrame(columns=self.patient_delta_df.columns)
                    fx_delta_df['UNIQUE_INDEX'] = np.array(range(len(log_layer))) + unique_index
                    fx_delta_df.index = fx_delta_df['UNIQUE_INDEX']
                    fx_delta_df.sort_index(inplace=True)
                    fx_delta_df.drop(columns=['UNIQUE_INDEX'], inplace=True)
                    fx_delta_df['LAYER_ID'] = layer_id
                    fx_delta_df['SPOT_ID'] = [id for id in range(len(fx_delta_df))]
                    fx_delta_df['FRACTION_ID'] = fx_id
                    fx_delta_df['BEAM_ID'] = beam_id
                    fx_delta_df['GANTRY_ANGLE'] = gtr_angle
                    fx_delta_df['TOTAL_LAYERS'] = num_layers

                    if has_plan and not layer_missing:
                        fx_delta_df['DELTA_X(mm)'] = delta_x
                        fx_delta_df['DELTA_Y(mm)'] = delta_y
                        fx_delta_df['DELTA_MU'] = delta_mu
                        fx_delta_df['DELTA_E(MeV)'] = delta_e

                    to_concat.append(fx_delta_df)
                    unique_index += len(fx_delta_df)

                    if layer_id == (num_layers - 1):  # progress visualization
                        print('  ', '[' + (layer_id + 1) * '#' + (num_layers - layer_id - 1) * '-' + ']', end=f' Beam {beam_id} complete\n')
                    else:
                        print('  ', '[' + (layer_id + 1) * '#' + (num_layers - layer_id - 1) * '-' + ']', end=f' Layer {str(layer_id + 1).zfill(2)}/{str(num_layers).zfill(2)}\r')

            print(f'  ..Fraction {f + 1}/{self.num_fractions} complete..\n')

        # concatenate layer deltaframes
        print('  ..Concatenating..')
        self.patient_delta_df = pd.concat(to_concat, sort=True)
        self.patient_delta_df.dropna(axis=1, how='all', inplace=True)
        if len(self.patient_delta_df) != len(self.patient_record_df):
            print(f'  /!\ Deltaframe length does not match dataframe length ({len(self.patient_delta_df)} vs. {len(self.patient_record_df)})')

        # write deltaframe
        os.chdir(self.df_destination)
        print(f'''  ..Writing deltaframe to '{self.df_destination}' as .CSV.. ''')
        while True:   
            try:
                self.patient_delta_df.to_csv(f'patient_{self.patient_id}_log-plan_delta.csv')
                break
            except PermissionError:
                input('  Permission denied, close target file and press ENTER.. ')
        plt.show()
        print('Complete')


    def corrupted_maps(self, gtr_only=False):
        other_records, other_deltas, other_tunings = [], [], []
        for file in sorted(os.listdir(self.df_destination)):
            if file.__contains__(str(self.patient_id)) and file.__contains__('delta') and file.endswith('.csv'):
                print(f'''Found patient deltaframe '{file}', reading in..''')
                self.patient_delta_df = pd.read_csv(os.path.join(self.df_destination, file), index_col='UNIQUE_INDEX', dtype={'BEAM_ID':str, 'FRACTION_ID':str})
            if file.__contains__('records') and file.endswith('.csv'):
                other_records.append(os.path.join(self.df_destination, file))           
            if file.__contains__('delta') and file.endswith('.csv'):
                other_deltas.append(os.path.join(self.df_destination, file))   
            if file.__contains__('tuning') and file.endswith('.csv'):
                other_tunings.append(os.path.join(self.df_destination, file))   
        
        to_drop, to_concat_rec, to_concat_tun = ['TOTAL_LAYERS'], [], []

        print('Gathering data from patient database..')
        for this_record, this_tuning in zip(other_records, other_tunings):
            has_delta = False
            this_patient_id = this_record.split('\\')[-1].split('_')[1]
            for this_delta in other_deltas:
                this_delta_id = this_delta.split('\\')[-1].split('_')[1]
                if this_patient_id == this_delta_id:
                    this_record_df = pd.read_csv(this_record, index_col='TIME', dtype={'BEAM_ID':str, 'FRACTION_ID':str})
                    this_tuning_df = pd.read_csv(this_tuning, index_col='TIME', dtype={'BEAM_ID':str, 'FRACTION_ID':str})
                    this_joint_df = pd.read_csv(this_delta, index_col='UNIQUE_INDEX', dtype={'BEAM_ID':str, 'FRACTION_ID':str})
                    this_joint_df['PATIENT_ID'] = int(this_patient_id)
                    this_tuning_df['PATIENT_ID'] = int(this_patient_id)
                    has_delta = True
                    break

            if not has_delta:
                print(f'  /!\ No deltaframe found for patient-ID {this_patient_id}, skipping..')
                continue

            if not gtr_only:
                if this_record_df.shape[0] == this_joint_df.shape[0]:
                    this_joint_df.drop(columns=to_drop, inplace=True)
                    this_joint_df['LAYER_ENERGY(MeV)'] = this_record_df['LAYER_ENERGY(MeV)'].to_list()
                    this_joint_df['MU'] = this_record_df['MU'].to_list()
                    this_joint_df['X_POSITION(mm)'] = this_record_df['X_POSITION(mm)'].to_list()
                    this_joint_df['Y_POSITION(mm)'] = this_record_df['Y_POSITION(mm)'].to_list()
                    this_joint_df['X_WIDTH(mm)'] = this_record_df['X_WIDTH(mm)'].to_list()
                    this_joint_df['Y_WIDTH(mm)'] = this_record_df['Y_WIDTH(mm)'].to_list()
                    this_joint_df['PRESSURE(hPa)'] = this_record_df['PRESSURE(hPa)'].to_list()
                    this_joint_df['TEMPERATURE(K)'] = this_record_df['TEMPERATURE(K)'].to_list()
                    to_concat_rec.append(this_joint_df)
                    to_concat_tun.append(this_tuning_df)
                else:
                    print(f'  /!\ Dataframe shapes do not match [{this_record_df.shape} vs. {this_joint_df.shape}]')
                    # continue
            
            else:
                this_joint_df.drop(columns=to_drop, inplace=True)
                to_concat_rec.append(this_joint_df)
                to_concat_tun.append(this_tuning_df)
        
        joint_df = pd.concat(to_concat_rec, ignore_index=True)
        joint_tu = pd.concat(to_concat_tun, ignore_index=True)
        joint_df['DIST(mm)'] = np.sqrt(joint_df['DELTA_X(mm)']**2 + joint_df['DELTA_Y(mm)']**2)
        joint_df['PLAN_X(mm)'] = joint_df['X_POSITION(mm)'] - joint_df['DELTA_X(mm)']
        joint_df['PLAN_Y(mm)'] = joint_df['Y_POSITION(mm)'] - joint_df['DELTA_Y(mm)']
        joint_df['PLAN_MU'] = joint_df['MU'] - joint_df['DELTA_MU']

        # print(joint_df.loc[abs(joint_df['DELTA_MU']) > 0.2][['PATIENT_ID', 'FRACTION_ID', 'BEAM_ID', 'LAYER_ID', 'SPOT_ID', 'DELTA_MU']])
        # plt.hist(joint_df['DELTA_MU'], bins=2000, edgecolor='black')
        # plt.yscale('log')
        # plt.xlabel('MU')
        # plt.ylabel('Count')
        # plt.show()
        # return None

        # Check all layer spotmaps for outliers
        failed, x_fail, y_fail, total, crit = 0, 0, 0, 0, 1.
        for id in joint_df.PATIENT_ID.drop_duplicates():
            # if id != 1230180: continue
            pat_df = joint_df.loc[joint_df.PATIENT_ID == id]
            pat_tu = joint_tu.loc[joint_tu.PATIENT_ID == id]
            for fx in pat_df.FRACTION_ID.drop_duplicates():
                fx_df = pat_df.loc[pat_df.FRACTION_ID == fx]
                fx_tu = pat_tu.loc[pat_tu.FRACTION_ID == fx]
                for beam in fx_df.BEAM_ID.drop_duplicates():
                    beam_df = fx_df.loc[fx_df.BEAM_ID == beam]
                    beam_tu = fx_tu.loc[fx_tu.BEAM_ID == beam]
                    for layer in beam_df.LAYER_ID.drop_duplicates():
                        layer_df = beam_df.loc[beam_df.LAYER_ID == layer]
                        layer_tu = beam_tu.loc[beam_tu.LAYER_ID == layer]
                        mean_dist = layer_df['DIST(mm)'].mean()
                        mean_x = layer_df['DELTA_X(mm)'].mean()
                        mean_y = layer_df['DELTA_Y(mm)'].mean()
                        max_delta = max(abs(layer_df['DELTA_X(mm)'].max()), abs(layer_df['DELTA_Y(mm)'].max()))
                        max_mu = max(abs(layer_df['DELTA_MU'].max()), abs(layer_df['DELTA_MU'].min()))
                        if (abs(mean_x) > crit or abs(mean_y) > crit):  # layer isotropically shifted
                            if abs(mean_x) > abs(mean_y):
                                x_fail += 1
                            else:
                                y_fail += 1
                                # plt.plot(layer_df['PLAN_X(mm)'], layer_df['PLAN_Y(mm)'], 'x-', c='tab:blue', label='Planned')
                                # plt.plot(layer_df['X_POSITION(mm)'], layer_df['Y_POSITION(mm)'], 'o-', c='black', markerfacecolor='None', label='Log file')
                                # plt.scatter(layer_tu['X_POSITION(mm)'], layer_tu['Y_POSITION(mm)'], c='lightgreen', label=f'Tuning pulses: {len(layer_tu)}')
                                # plt.legend(title=f'Mean dist: {mean_dist:.3f} mm')
                                # plt.title(f'{id} | Date: {fx} | Beam: {beam} ({layer_df.GANTRY_ANGLE.iloc[0]}) | Layer: {layer + 1}')
                                # # plt.xlim(x_min, x_max)
                                # # plt.xlim(y_min, y_max)
                                # plt.tight_layout()
                                # plt.show()
                                # print(len(layer_tu) > 1)
                                failed += 1
                        if max_mu > 0.065:
                            # plt.plot(layer_df['PLAN_X(mm)'], layer_df['PLAN_Y(mm)'], 'x-', c='tab:blue', label='Planned', zorder=1)
                            # plt.scatter(layer_df['X_POSITION(mm)'], layer_df['Y_POSITION(mm)'], c=layer_df['DELTA_MU'], cmap='bwr', vmin=-0.1, vmax=0.1, zorder=3)
                            # plt.plot(layer_df['X_POSITION(mm)'], layer_df['Y_POSITION(mm)'], 'o-',  c='black', markerfacecolor='white', label='Log file', zorder=2)
                            # plt.colorbar(label='$\Delta$MU to plan')
                            # plt.legend(title=f'Max MU diff: {max_mu:.4f} MU')
                            # plt.title(f'{id} | Date: {fx} | Beam: {beam} | Layer-ID: {layer}')
                            # # plt.xlim(x_min, x_max)
                            # # plt.xlim(y_min, y_max)
                            # plt.tight_layout()
                            # print(layer_df[['MU', 'PLAN_MU']].iloc[43:45])
                            # plt.show()
                            pass
                        total += 1
        
        # print(f'{failed}/{total} ({failed/total * 100:.1f}%) layers failed {crit} mm mean distance criterion')
        print(f'{x_fail + y_fail}/{total} ({failed/total * 100:.3f}%) of layers failed {crit} mm delta x/y criterion')
        print(f'From this, {x_fail} ({x_fail/(x_fail + y_fail) * 100:.1f}%) are caused by X error')
        print(f'From this, {y_fail} ({y_fail/(x_fail + y_fail) * 100:.1f}%) are caused by Y error')


    def delta_dependencies(self):
        for file in sorted(os.listdir(self.df_destination)):
            if file.__contains__(str(self.patient_id)) and file.__contains__('delta') and file.endswith('.csv'):
                print(f'''Found patient deltaframe '{file}', reading in..''')
                self.patient_delta_df = pd.read_csv(os.path.join(self.df_destination, file), index_col='UNIQUE_INDEX', dtype={'BEAM_ID':str, 'FRACTION_ID':str})
                break            

        try:
            if len(self.patient_record_df) != len(self.patient_delta_df):
                print(f'  /!\ Dataframe length does not match deltaframe length [{len(self.patient_record_df)} vs. {len(self.patient_record_df)}], proceed with caution..')
        except AttributeError:
            print(f'''\nUnable to locate patient deltaframe for patient-ID {self.patient_id}, calling prepare_deltaframe()..''')
            self.prepare_deltaframe()
        
        options = ['Histograms', 'Gantry angle vs. delta(x,y)', 'MU vs. delta(x,y)', '(x,y) vs. delta(x,y)', 'Energy vs. delta(x,y)', 'MU(plan) vs. MU(log)']
        print('Choose dependency to derive:')
        for i, option in enumerate(options):
            print(f'  ({i + 1})\t{option}')
        
        while True:
            try:
                choice = abs(int(input('\n Selected option: ')))
                if choice <= len(options):
                    break
                else:
                    print('  /x\ Selection not in list, retry..')

            except TypeError and ValueError:
                print('  /x\ Invalid input, retry..')

        print(f'''\nSelected option is ({choice}) {options[choice - 1]}. Starting..''')

        if choice == 1:
            beams = self.patient_delta_df.loc[~self.patient_delta_df['DELTA_X(mm)'].isna(), 'BEAM_ID'].drop_duplicates().to_list()  # only beams with existing delta
            fig, axs = plt.subplots(len(beams), 4, figsize=(12, 3 * len(beams)), sharex='col')
            ax0 = fig.add_subplot(111, frameon=False)
            ax0.set_xticks([])
            ax0.set_yticks([])
            # axs = axs.flatten()
            bins = 100
            bin_xy = np.histogram(np.hstack((self.patient_delta_df.loc[~self.patient_delta_df['DELTA_X(mm)'].isna(), 'DELTA_X(mm)'], self.patient_delta_df.loc[~self.patient_delta_df['DELTA_X(mm)'].isna(), 'DELTA_Y(mm)'])), bins=bins)[1]
            bin_mu = np.histogram(self.patient_delta_df.loc[~self.patient_delta_df['DELTA_X(mm)'].isna(), 'DELTA_MU'], bins=bins)[1]
            bin_e = np.histogram(self.patient_delta_df.loc[~self.patient_delta_df['DELTA_X(mm)'].isna(),'DELTA_E(MeV)'], bins=bins)[1]

            for row, beam in enumerate(beams):
                beam_df = self.patient_delta_df.loc[self.patient_delta_df['BEAM_ID'] == beam]

                axs[row, 0].hist(beam_df['DELTA_X(mm)'], bins=bin_xy, alpha=0.7, label=f'''$\mu_x =$ {beam_df['DELTA_X(mm)'].mean():.3f} mm\n$\sigma_x =$ {beam_df['DELTA_X(mm)'].std():.3f} mm''')
                axs[row, 0].axvline(beam_df['DELTA_X(mm)'].mean(), ls='-', color='black', lw=0.5)
                axs[row, 0].axvline(0.0, ls='--', color='black', lw=0.5)
                axs[row, 0].set_xlabel('$\Delta x$ to plan [mm]')
                axs[row, 0].set_ylabel(f'Beam {beam} - {beam_df.GANTRY_ANGLE.iloc[0]}°')
                axs[row, 0].set_yscale('log')
                # axs[row, 0].set_xlim(-2, 2)
                axs[row, 1].hist(beam_df['DELTA_Y(mm)'], bins=bin_xy, alpha=0.7, label=f'''$\mu_y =$ {beam_df['DELTA_Y(mm)'].mean():.3f} mm\n$\sigma_y =$ {beam_df['DELTA_Y(mm)'].std():.3f} mm''')
                axs[row, 1].axvline(beam_df['DELTA_Y(mm)'].mean(), ls='-', color='black', lw=0.5)
                axs[row, 1].axvline(0.0, ls='--', color='black', lw=0.5)
                axs[row, 1].set_xlabel('$\Delta y$ to plan [mm]')
                axs[row, 1].set_yscale('log')
                # axs[row, 1].set_xlim(-2, 2)
                axs[row, 2].hist(beam_df['DELTA_MU'], bins=bin_mu, color='tab:green', alpha=0.7, label=f'''$\mu_D =$ {beam_df['DELTA_MU'].mean():.4f} MU\n$\sigma_D =$ {beam_df['DELTA_MU'].std():.4f} MU''')
                axs[row, 2].axvline(beam_df['DELTA_MU'].mean(), ls='-', color='black', lw=0.5)
                axs[row, 2].axvline(0.0, ls='--', color='black', lw=0.5)
                axs[row, 2].set_xlabel('$\Delta D$ to plan [MU]')
                axs[row, 2].set_yscale('log')
                # axs[row, 2].set_xlim(-0.005, 0.005)
                # axs[2].set_ylim(0, 30000)
                axs[row, 3].hist(beam_df['DELTA_E(MeV)'].drop_duplicates(), bins=bin_e, color='tab:red', alpha=0.7, label=f'''$\mu_E =$ {beam_df['DELTA_E(MeV)'].mean():.4f} MeV\n$\sigma_E =$ {beam_df['DELTA_E(MeV)'].std():.4f} MeV''')
                axs[row, 3].axvline(beam_df['DELTA_E(MeV)'].mean(), ls='-', color='black', lw=0.5)
                axs[row, 3].axvline(0.0, ls='--', color='black', lw=0.5)
                axs[row, 3].set_xlabel('$\Delta E$ to plan [MeV]')
                axs[row, 3].set_yscale('log')
                # axs[row, 3].set_xlim(-0.04, 0.04)

                print(f'>> BEAM {beam} <<')
                print(f'''dx[mm]:\t\tmin={round(beam_df['DELTA_X(mm)'].min(), 3)}\tmax={round(beam_df['DELTA_X(mm)'].max(), 3)}\tmean={round(beam_df['DELTA_X(mm)'].mean(), 3)}\tstd={round(beam_df['DELTA_X(mm)'].std(), 3)}''')
                print(f'''dy[mm]:\t\tmin={round(beam_df['DELTA_Y(mm)'].min(), 3)}\tmax={round(beam_df['DELTA_Y(mm)'].max(), 3)}\tmean={round(beam_df['DELTA_Y(mm)'].mean(), 3)}\tstd={round(beam_df['DELTA_Y(mm)'].std(), 3)}''')
                print(f'''dD[MU]:\t\tmin={round(beam_df['DELTA_MU'].min(), 4)}\tmax={round(beam_df['DELTA_MU'].max(), 4)}\tmean={round(beam_df['DELTA_MU'].mean(), 4)}\tstd={round(beam_df['DELTA_MU'].std(), 4)}''')
                print(f'''dE[MeV]:\tmin={round(beam_df['DELTA_E(MeV)'].min(), 4)}\tmax={round(beam_df['DELTA_E(MeV)'].max(), 4)}\tmean={round(beam_df['DELTA_E(MeV)'].mean(), 4)}\tstd={round(beam_df['DELTA_E(MeV)'].std(), 4)}''')

                for j in range(4):
                    if j == 3: which = 'both'
                    else: which = 'major'
                    axs[row, j].grid(which=which, axis='y', zorder=-1)
                    axs[row, j].set_axisbelow(True)
                    axs[row, j].legend()

            # ax0.set_title(f'Delta Histograms for Patient-ID {self.patient_id}', fontweight='bold')
            plt.title(f'Log-file report for patient-ID {self.patient_id}', fontweight='bold')
            plt.tight_layout()
            plt.savefig(f'{output_dir}/{self.patient_id}_histograms.png', dpi=300)
            # plt.show()

        if choice == 2:  # gantry angle vs. delta(x,y)
            other_dfs = [pd.read_csv(os.path.join(self.df_destination, file), index_col='UNIQUE_INDEX', dtype={'BEAM_ID':str, 'FRACTION_ID':str}) for file in sorted(os.listdir(self.df_destination)) if file.__contains__('delta') and file.endswith('.csv')]
            print('  Concatenating existent patient dataframes..')
            concat_df = pd.concat([df for df in other_dfs])
            
            gtr_angles = sorted(concat_df['GANTRY_ANGLE'].drop_duplicates().to_list())

            # fig, axs = plt.subplots(2, 1, sharey=True, figsize=(10, 10))
            # for alpha in gtr_angles:
            #     delta_x, delta_y = concat_df.loc[concat_df['GANTRY_ANGLE'] == alpha, 'DELTA_X(mm)'], concat_df.loc[concat_df['GANTRY_ANGLE'] == alpha ,'DELTA_Y(mm)']
            #     dx_mean, dy_mean = delta_x.mean(), delta_y.mean()
            #     axs[0].boxplot(delta_x, positions=[alpha], widths = 20)
            #     axs[0].set_ylabel('$\mathbf{\Delta x}$ [mm]', fontweight='bold')
            #     axs[1].boxplot(delta_y, positions=[alpha], widths = 20)
            #     axs[1].set_ylabel('$\mathbf{\Delta (x,y)}$ [mm]', fontweight='bold')
            #     axs[1].set_xlabel('Gantry angle [°]', fontweight='bold')
            dx_means, dy_means = [], []
            for alpha in gtr_angles:
                delta_x, delta_y = concat_df.loc[concat_df['GANTRY_ANGLE'] == alpha, 'DELTA_X(mm)'], concat_df.loc[concat_df['GANTRY_ANGLE'] == alpha ,'DELTA_Y(mm)']
                dx_mean, dy_mean = delta_x.mean(), delta_y.mean()
                dx_means.append(dx_mean), dy_means.append(dy_mean)
            plt.axhline(0.0, color='grey', lw=0.5)
            plt.plot(gtr_angles, dx_means, marker='o', c='black', ls='-', label='$\Delta x$ (mean)')
            plt.plot(gtr_angles, dy_means, marker='o', c='black', ls='--', label='$\Delta y$ (mean)', markerfacecolor='white')
            plt.xlabel('Gantry angle [°]')
            plt.ylabel('$\Delta (x,y)$ [mm]')
            plt.legend()
            plt.savefig(f'{output_dir}/delta_vs_gantry_angle.png', dpi=150)
            plt.show()
        
        if choice == 4:  # (x,y) vs. delta(x,y)
            x, y = self.patient_record_df['X_POSITION(mm)'], self.patient_record_df['Y_POSITION(mm)']
            delta_x, delta_y = self.patient_delta_df['DELTA_X(mm)'], self.patient_delta_df['DELTA_Y(mm)']
            fig, axs = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(10, 10))
            axs[0].plot(x, delta_x, marker='o', ls='-', label='$\Delta x$')
            axs[0].plot(x, delta_y, marker='o', ls='-', label='$\Delta y$')
            axs[0].set_xlabel('$\mathbf{x}$ position [mm]', fontweight='bold')
            axs[0].set_ylabel('$\mathbf{\Delta x}$ [mm]', fontweight='bold')
            axs[0].legend()
            axs[1].plot(y, delta_x, marker='o', ls='-', label='$\Delta x$')
            axs[1].plot(y, delta_y, marker='o', ls='-', label='$\Delta y$')
            axs[1].set_ylabel('$\mathbf{\Delta (x,y)}$ [mm]', fontweight='bold')
            axs[1].set_xlabel('$\mathbf{y}$ position [mm]', fontweight='bold')
            axs[1].legend()

            plt.show()
        
        if choice == 6:  # MU(plan) vs. MU(log)
            log_mu, delta_mu = self.patient_record_df['MU'].to_numpy(), self.patient_delta_df['DELTA_MU'].to_numpy()
            plan_mu = log_mu + delta_mu
            lim = max(max(plan_mu), max(log_mu))
            lim = lim + 0.05 * lim
            fig, axs = plt.subplots(3, 1, sharex=True, figsize=(7, 10))
            ax0 = fig.add_subplot(111, frameon=False)
            ax0.set_xticks([]), ax0.set_yticks([])     
            axs[0].plot([0, lim], [0, lim], ls='--', lw=0.2, c='black')    
            axs[0].plot(plan_mu, log_mu, c='tab:red', marker='x', linestyle='none', markersize=1, markeredgewidth=0.5, label='Spot meterset weights')
            axs[0].set_ylabel('MU log-file')
            axs[0].set_xlim(0, lim)
            axs[0].set_ylim(0, lim)
            axs[0].legend()
            axs[1].plot(plan_mu, log_mu - plan_mu, c='tab:blue', marker='x', linestyle='none', markersize=1, markeredgewidth=0.5, label='Spot meterset difference')
            axs[1].set_ylabel('MU difference (Log - Plan)')
            axs[1].axhline(0.0, ls='--', lw=0.2, color='black')
            axs[1].set_ylim(-0.05, 0.02)
            axs[1].legend()
            axs[2].plot(plan_mu, np.divide(log_mu, plan_mu), c='tab:orange', marker='x', linestyle='none', markersize=1, markeredgewidth=0.5, label='Spot meterset ratio')
            axs[2].axhline(1.0, ls='--', lw=0.2, color='black')
            axs[2].set_xlabel('MU planned')
            # axs[1].set_ylabel('$\Delta$MU (plan - log)')
            axs[2].set_ylabel('MU ratio (Log/Plan)')
            axs[2].set_ylim(0.65, 1.1)
            axs[2].legend()
            ax0.set_title(f'Meterset weight comparison for patient-ID {self.patient_id}', fontweight='bold', fontsize=10)
            plt.savefig(f'{output_dir}/{self.patient_id}_spot_metersets.png', dpi=600)
            plt.show()

    
    def beam_timings(self):
        logfile_root = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\Logfiles\converted'
        beam_list = self.patient_record_df['BEAM_ID'].drop_duplicates()
        indices = beam_list.index.to_list()

        # get body site of patient
        bs_file = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\cct_entities.csv'
        bs_data = pd.read_csv(bs_file, sep=';', index_col='PATIENT_ID')
        bs_dict = bs_data.to_dict()
        print('_______')
        print(self.patient_id, bs_dict['BODY_SITE'][int(self.patient_id)])

        print('\nSelect beam key for beam timing plot:\n')
        print('Key\tBeam-ID\t\tGantry Angle\tLayers\tFrom Fraction\n')  # choose beam
        for choice, (index, beam) in enumerate(zip(indices, beam_list)):
            num_layers = int(self.patient_record_df['TOTAL_LAYERS'][index].mean())
            fraction = self.patient_record_df['FRACTION_ID'][index]
            try:
                gantry_angle = self.patient_record_df['GANTRY_ANGLE'][index].iloc[0]
            except:
                gantry_angle = self.patient_record_df['GANTRY_ANGLE'][index]
            print(f'({choice + 1})\t{beam}\t\t{gantry_angle}\t\t{num_layers}\t{fraction}')

        while True:
            try:
                key = input('\n Select beam key(s) separated by comma: ')
                key = [int(k) for k in key.split(',')]
                if max(key) > len(beam_list) or min(key) <= 0:
                    print('Key out of bounds, select another..')
                    continue
                else:
                    break
            except:
                print('Invalid input, try again..')

        beam_ids = [str(beam_list[k - 1]) for k in key]
        print(beam_ids)

        cmap = sns.color_palette('inferno_r', 6)
        sns.set(context='paper', style='ticks', font_scale=1.5)
        sns.set_style({"xtick.direction": "in","ytick.direction": "in", "axes.axisbelow":"line"})
        global_plan_delivery_times = [0 for fx in self.fraction_list]
        for beam_id in beam_ids:
            fig, axs = plt.subplots(4, 1, sharex=True, figsize=(7, 7), gridspec_kw={'height_ratios': [3, 1, 1, 1]})
            # fig.subplots_adjust(hspace=0.0, wspace=0.0)
            ax0 = fig.add_subplot(111, frameon=False)
            ax0.set_xticks([]), ax0.set_yticks([])
            ax0.set_ylabel('Time [s]', labelpad=50)
            ax0.set_xlabel('Fraction', labelpad=30)
            # ax0.set_title(f'Beam timings for patient-ID {self.patient_id}', fontweight='bold')
            # axs.flatten()
            x_axis = self.fraction_list
            # for ax in axs:
            #     ax.set_xticks(range(len(x_axis)))
            #     ax.set_xticklabels(range(1, len(x_axis) + 1))

            print('\nGenerating beam timing plot..')
            global_drills, global_spot_switches, global_energy_switches, global_interlocks, totals = [], [], [], [], []
            for x, fx_id in enumerate(x_axis):  # remember type(fx_id) = <str>
                if self.patient_id == '1180747' and x == 20: continue
                if self.patient_id == '1635796' and x == 10: continue
                beam_df = self.patient_record_df.loc[(self.patient_record_df['BEAM_ID'] == beam_id) & (self.patient_record_df['FRACTION_ID'] == fx_id)]
                beam_tuning_df = self.patient_tuning_df.loc[(self.patient_tuning_df['BEAM_ID'] == beam_id) & (self.patient_tuning_df['FRACTION_ID'] == fx_id)]
                if beam_df.empty or beam_tuning_df.empty:
                    continue

                total_drill_time = (beam_df['DRILL_TIME(ms)'].sum() + beam_tuning_df['DRILL_TIME(ms)'].sum()) / 1000
                total_layer_time, total_energy_switch = 0.0, 0.0
                layer_dfs = [beam_df.loc[beam_df['LAYER_ID'] == lid] for lid in beam_tuning_df['LAYER_ID'].drop_duplicates()]
                layer_tuning_dfs = [beam_tuning_df.loc[beam_tuning_df['LAYER_ID'] == lid] for lid in beam_tuning_df['LAYER_ID'].drop_duplicates()]
                for layer_id, layer_df in enumerate(layer_dfs):
                    if layer_df.empty:
                        layer_df = layer_tuning_dfs[layer_id]  # all layer spots covered by tuning

                    start_irr = pd.to_datetime(layer_df.first_valid_index(), format='%Y-%m-%d %H:%M:%S.%f')
                    start_tun = pd.to_datetime(layer_tuning_dfs[layer_id].first_valid_index(), format='%Y-%m-%d %H:%M:%S.%f')
                    end_irr = pd.to_datetime(layer_df.last_valid_index(), format='%Y-%m-%d %H:%M:%S.%f')
                    layer_time = end_irr - start_irr
                    total_layer_time += layer_time.total_seconds()
                    
                    if layer_id > 0:
                        if layer_dfs[layer_id - 1].empty:
                            end_prev = pd.to_datetime(layer_tuning_dfs[layer_id - 1].last_valid_index(), format='%Y-%m-%d %H:%M:%S.%f')
                        else:
                            end_prev = pd.to_datetime(layer_dfs[layer_id - 1].last_valid_index(), format='%Y-%m-%d %H:%M:%S.%f')

                        energy_switch = start_tun - end_prev
                        total_energy_switch += energy_switch.total_seconds()
                
                total_spot_switch = total_layer_time - total_drill_time
                
                fx_dir = os.path.join(logfile_root, self.patient_id, fx_id, beam_id)
                has_interlock, total_interlock = False, 0.0
                for file in os.listdir(fx_dir):
                    if file.__contains__('events') and file.endswith('.csv'):
                        events = open(os.path.join(fx_dir, file)).readlines()
                        for i, line in enumerate(events):
                            if line.__contains__('Layer'):
                                layer_id = int(line.split('Layer ')[-1].split(':')[0])
                            if (line.__contains__('RESUME_REQUESTED') or line.__contains__('PAUSE_REQUESTED')) and not has_interlock:  # dont count resume statements where nothing happened in between
                                print(f'  /!\ Interlock detected in fraction {x + 1} ({fx_id})')
                                stop_time = pd.to_datetime(events[i - 1].split(',')[0], format='%d/%m/%Y %H:%M:%S.%f')
                                has_interlock = True

                            if has_interlock:
                                if line.__contains__('START_MAP_IRRADIATION'):
                                    resume_time = pd.to_datetime(line.split(',')[0], format='%d/%m/%Y %H:%M:%S.%f')
                                    interlock_layer_df = beam_df.loc[beam_df.LAYER_ID == layer_id]
                                    time_index = interlock_layer_df.index.to_series()
                                    time_diff = interlock_layer_df.index.to_series().diff()  # Calculate the time difference between consecutive rows
                                    time_df = pd.DataFrame({'t':time_index.to_list(), 'dt':time_diff.to_list()})
                                    time_df['post_stop'] = time_df.t > stop_time
                                    time_df = time_df.loc[time_df.post_stop == True]

                                    # Check if any time differences are greater than 1 second
                                    if (time_diff > pd.Timedelta(seconds=1)).any():  # Interlock during spot switch
                                        print('      (Spot switch)')
                                        resume_idx = time_diff.index.get_loc(time_df.loc[(time_df.dt > pd.Timedelta(seconds=1)), 't'].iloc[0])
                                        resume_time = time_index.iloc[resume_idx]
                                        print(f'      {resume_time}')
                                        stop_time = time_index.iloc[resume_idx - 1]
                                        # resume_time = pd.to_datetime(time_df.loc[time_df.dt > pd.Timedelta(seconds=1)].t.iloc[0],format='%Y-%m-%d %H:%M:%S.%f')  # overwrite resume time if spot switch (more precise)
                                        interlock_time = (resume_time - stop_time).total_seconds()
                                        if global_spot_switches:
                                            if abs(total_spot_switch - global_spot_switches[-1]) > 0.1:
                                                total_spot_switch -= interlock_time
                                        else:
                                            total_spot_switch -= interlock_time
                                    else:  # Interlock during energy switch
                                        print('      (Energy switch)')
                                        interlock_time = (resume_time - stop_time).total_seconds()
                                        if global_energy_switches:
                                            if total_energy_switch - global_energy_switches[-1] > 3:
                                                total_energy_switch -= interlock_time
                                        else:
                                            total_energy_switch -= interlock_time
                                                                        
                                    total_interlock += interlock_time
                                    has_interlock = False  # interlock handled
                                   
                global_drills.append(total_drill_time), global_spot_switches.append(total_spot_switch), global_energy_switches.append(total_energy_switch)
                if total_interlock > 0:
                    global_interlocks.append(total_interlock)
                else:
                    global_interlocks.append(np.nan)
                totals.append(total_drill_time + total_spot_switch + total_energy_switch)
                global_plan_delivery_times[x] += (total_drill_time + total_spot_switch + total_energy_switch)
            

            global_drills = np.array(global_drills)
            global_spot_switches = np.array(global_spot_switches)
            global_energy_switches = np.array(global_energy_switches)
            global_interlocks = np.array(global_interlocks)

            # print(f'\n\t\tMean\tMin\tMax\tSigma')
            # print(f'Drill\t\t{np.mean(global_drills):.3f}\t{min(global_drills):.3f}\t{max(global_drills):.3f}\t{np.std(global_drills):.3f}')
            # print(f'Spot switch\t{np.mean(global_spot_switches):.3f}\t{min(global_spot_switches):.3f}\t{max(global_spot_switches):.3f}\t{np.std(global_spot_switches):.3f}')
            # print(f'Energy switch\t{np.mean(global_energy_switches):.3f}\t{min(global_energy_switches):.3f}\t{max(global_energy_switches):.3f}\t{np.std(global_energy_switches):.3f}')
            # print(f'Total beamtime\t{np.mean(totals):.3f}\t{min(totals):.3f}\t{max(totals):.3f}\t{np.std(totals):.3f}\n')
            print(f'Mean field delivery time:\t{round(np.mean(totals), 2)} ({round(np.std(totals), 2)})\n')

            ec = 'none'
            axs[0].bar(range(1, len(global_drills) + 1), global_drills, label='Spot drill', color=cmap[0], edgecolor=ec, zorder=0)
            axs[0].bar(range(1, len(global_spot_switches) + 1), global_spot_switches, bottom=global_drills, label='Spot switching', color=cmap[1], edgecolor=ec, zorder=2)
            axs[0].bar(range(1, len(global_energy_switches) + 1), global_energy_switches, bottom=global_spot_switches + global_drills, label='Energy switching', color=cmap[2], edgecolor=ec, zorder=2)
            axs[0].bar(range(1, len(global_interlocks) + 1), global_interlocks, bottom=global_spot_switches + global_drills + global_energy_switches, label='Pause (interlock)', color=cmap[4], edgecolor=ec, zorder=2)
            axs[3].plot(range(1, len(global_drills) + 1), global_drills, marker='o', color='black', markerfacecolor=cmap[0], label='Spot drill time')
            axs[3].axhline(np.mean(global_drills), ls='--', color='black', lw=1.0, zorder=-1, label='Mean')
            axs[3].set_xticks([1, 5, 10, 15, 20, 25, 30, 35, 40])
            axs[2].plot(range(1, len(global_spot_switches) + 1), global_spot_switches, marker='s', color='black', markerfacecolor=cmap[1], label='Spot switching time')  
            axs[2].axhline(np.mean(global_spot_switches), ls='--', color='black', lw=1.0, zorder=-1, label='Mean')
            axs[1].plot(range(1, len(global_energy_switches) + 1), global_energy_switches, marker='^', color='black', markerfacecolor=cmap[2], label='Energy switching time')  
            axs[1].axhline(np.mean(global_energy_switches), ls='--', color='black', lw=1.0, zorder=-1, label='Mean')

            axs[0].set_ylim(0, 125)
            axs[1].set_ylim(29, 34)
            axs[1].set_yticks([29, 31, 33])
            axs[2].set_ylim(4.46, 4.51)
            axs[2].set_yticks([4.46, 4.48, 4.5])
            axs[3].set_ylim(6.2, 7.6)
            axs[3].set_yticks([6.2, 6.7, 7.2])

            for i, ax in enumerate(axs):
                if i == 0:
                    ax.legend(loc='upper left')
                    ax.grid(axis='y')
                else:
                    ax.legend(loc='upper left', ncol=2)
                ax.set_axisbelow(True)
                # ax.grid(axis='both')        
            plt.tight_layout()
            plt.savefig(f'{output_dir}/{self.patient_id}_{beam_id}_beamtime.pdf')        
            # plt.show()
        global_plan_delivery_times = [t for t in global_plan_delivery_times if t != 0]
        print('_____________________________________________________')
        print(f'Total beamtime\t{round(np.mean(global_plan_delivery_times), 2)} ({round(np.std(global_plan_delivery_times), 2)})\n')


    def plan_creator(self, fraction, mode):
        log_xlim = 150  # mm half maximum field extent
        log_ylim = 193
        if fraction == 'last':
            fx_list = [self.fraction_list[-1]]
        elif fraction == 'all':
            fx_list = [fx for fx in self.fraction_list]
        else: 
            root = Tk()
            fx_list = [filedialog.askdirectory(initialdir=self.logfile_dir).split('/')[-1]]
            root.destroy()
        
        for fx_no, fx_id in enumerate(fx_list):
            target_record = self.patient_record_df.loc[self.patient_record_df['FRACTION_ID'] == fx_id]
            target_tuning = self.patient_tuning_df.loc[self.patient_tuning_df['FRACTION_ID'] == fx_id]
            beam_list = target_record['BEAM_ID'].drop_duplicates().to_list()
            try:
                plan_dcm, beam_ds = self.dicom_finder(fx_id, beam_list[0])
            except TypeError:
                print(f'/!\ Fraction {fx_id} with beams {self.beam_list[fx_no]} will be skipped (NoPlanFoundError)')
                continue

            dcm_path = os.path.dirname(plan_dcm)
            
            print(f'\nWriting RT plan for fraction-ID {fx_id}..')
            ds = pydicom.read_file(plan_dcm)
            for i, (beam_id, plan_beam) in enumerate(zip(beam_list, ds.IonBeamSequence)):
                beam_df = target_record.loc[target_record['BEAM_ID'] == beam_id]
                if not mode == 'all':
                    sorting_dict = self.spot_sorter(fx_id, beam_id)
                    tuning_dict = self.spot_sorter(fx_id, beam_id, mode='tuning')

                total_layers = beam_df['TOTAL_LAYERS'][0]
                cumulative_mu = 0
                for layer_id in range(target_record.loc[target_record['BEAM_ID'] == beam_id]['TOTAL_LAYERS'].iloc[0]):
                    layer_df = target_record.loc[(target_record['BEAM_ID'] == beam_id) & (target_record['LAYER_ID'] == layer_id)]
                    if layer_df.empty: print(f'  /!\ Layer-ID {layer_id} has no record, probably one-spot layer')
                    tuning_df = target_tuning.loc[(target_tuning['BEAM_ID'] == beam_id) & (target_tuning['LAYER_ID'] == layer_id)]
                    layer_xy, layer_mu, tuning_xy, tuning_mu = [], [], [], []
                    # layer_energy = tuning_df['LAYER_ENERGY(MeV)'].drop_duplicates().iloc[0]
                    plan_xy = plan_beam.IonControlPointSequence[layer_id * 2].ScanSpotPositionMap

                    # generate spot position list equal to DICOM tag --> [x1, y1, x2, y2, ..]
                    for log_spot in range(len(layer_df['X_POSITION(mm)'].to_list())):
                        x = layer_df['X_POSITION(mm)'][log_spot]
                        y = layer_df['Y_POSITION(mm)'][log_spot]
                        if abs(x) > log_xlim:
                            print(f'  /!\ X-limit exceeded by {abs(x) - log_xlim}mm - Spot-ID in Fx-ID {fx_id} beam {beam_id} layer {layer_id}, setting spot to border of field..')
                            x = np.sign(x) * log_xlim
                        if abs(y) > log_ylim:
                            print(f'  /!\ Y-limit exceeded by {abs(y) - log_ylim}mm - Spot-ID in Fx-ID {fx_id} beam {beam_id} layer {layer_id}, setting spot to border of field..')
                            y = np.sign(y) * log_ylim
                        layer_xy.append(x)
                        layer_xy.append(y)
                        layer_mu.append(layer_df['MU'][log_spot])
                    for tune_spot in range(len(tuning_df['X_POSITION(mm)'].to_list())):
                        tuning_xy.append(tuning_df['X_POSITION(mm)'][tune_spot])
                        tuning_xy.append(tuning_df['Y_POSITION(mm)'][tune_spot])
                        tuning_mu.append(tuning_df['MU'][tune_spot])
                                      
                    n_log_spots = len(layer_mu)
                    n_plan_spots = plan_beam.IonControlPointSequence[layer_id * 2].NumberOfScanSpotPositions

                    if mode == 'all':  # except energy
                        if layer_id == 0:
                            print(f'  Will overwrite planned spot positions and metersets for beam-iD {beam_id}..')
                        
                        # combine record and tuning spots
                        layer_xy += tuning_xy
                        layer_mu += tuning_mu
                        n_log_spots = len(layer_mu)
                        spot_diff = n_log_spots - n_plan_spots
                        if layer_id != 0 and spot_diff > 1:  # allow 1 spot difference for all layers except the first (which might have multiple tuning spots)
                            print(f'  /!\ Warning: {spot_diff} more spots recorded than planned in layer-ID {layer_id}, fraction-ID {fx_id}\n      > Check for additional tuning pulses')
                            if spot_diff > 2:
                                plt.plot([x for x in layer_xy if layer_xy.index(x) % 2 == 0], [y for y in layer_xy if layer_xy.index(y) % 2 != 0], 'o')
                                plt.plot([x for x in plan_xy if plan_xy.index(x) % 2 == 0], [y for y in plan_xy if plan_xy.index(y) % 2 != 0], 'x')
                                plt.show()
                        elif layer_id != 0 and spot_diff < 0:
                            print(f'  /x\ Critical: {abs(spot_diff)} less spots recorded than planned in layer-ID {layer_id}, fraction-ID {fx_id}\n      > Review log-files')
                            plt.plot([x for x in layer_xy if layer_xy.index(x) % 2 == 0], [y for y in layer_xy if layer_xy.index(y) % 2 != 0], 'o')
                            plt.plot([x for x in plan_xy if plan_xy.index(x) % 2 == 0], [y for y in plan_xy if plan_xy.index(y) % 2 != 0], 'x')
                            plt.show()

                        plan_beam.IonControlPointSequence[layer_id * 2].NumberOfScanSpotPositions = plan_beam.IonControlPointSequence[layer_id * 2 + 1].NumberOfScanSpotPositions = n_log_spots
                        plan_beam.IonControlPointSequence[layer_id * 2].ScanSpotPositionMap = plan_beam.IonControlPointSequence[layer_id * 2 + 1].ScanSpotPositionMap = layer_xy
                        plan_beam.IonControlPointSequence[layer_id * 2].CumulativeMetersetWeight = cumulative_mu
                        cumulative_mu += sum(layer_mu)
                        plan_beam.IonControlPointSequence[layer_id * 2 + 1].CumulativeMetersetWeight = cumulative_mu
                        plan_beam.IonControlPointSequence[layer_id * 2].ScanSpotMetersetWeights = layer_mu
                        plan_beam.IonControlPointSequence[layer_id * 2 + 1].ScanSpotMetersetWeights = [0.0 for _ in range(len(layer_mu))]
                        # plan_beam.IonControlPointSequence[layer_id * 2].NominalBeamEnergy = plan_beam.IonControlPointSequence[layer_id * 2 + 1].NominalBeamEnergy = layer_energy
                    
                    elif mode == 'pos':
                        if layer_id == 0:
                            print(f'  Will overwrite only planned spot positions for beam-iD {beam_id}..')

                        record_xy_tuples = [(layer_xy[i], layer_xy[i + 1]) for i in range(len(layer_xy)) if i % 2 == 0]
                        xy_sorted = [record_xy_tuples[sorting_dict[layer_id][i]] for i in range(len(record_xy_tuples))]
                        
                        if n_log_spots < n_plan_spots:  # tuning spot replaces one planned spot, map it to closest plan spot
                            tuning_xy_tuples = [np.array((tuning_xy[i], tuning_xy[i + 1])) for i in range(len(tuning_xy)) if i % 2 == 0]
                            avg_tuning_pos = sum(tuning_xy_tuples) / len(tuning_xy_tuples)
                            xy_sorted.insert(tuning_dict[layer_id], avg_tuning_pos)                            
                        
                        if n_log_spots > n_plan_spots:
                            print(f'  /!\ Critical: More spots recorded than planned in layer-ID {layer_id} ({n_log_spots} vs. {n_plan_spots}), skipping this..')
                            continue

                        layer_xy = []
                        for tup in xy_sorted:
                            layer_xy.append(tup[0])
                            layer_xy.append(tup[1])
                        
                        plan_beam.IonControlPointSequence[layer_id * 2].ScanSpotPositionMap = layer_xy
                        plan_beam.IonControlPointSequence[layer_id * 2 + 1].ScanSpotPositionMap = layer_xy

                    elif mode == 'mu':
                        if layer_id == 0:
                            print(f'  Will overwrite only planned spot metersets for beam-iD {beam_id}..')
                        
                        plan_xy_tuples = [np.array((plan_xy[i], plan_xy[i + 1])) for i in range(len(plan_xy)) if i % 2 == 0]
                        mu_sorted = [layer_mu[sorting_dict[layer_id][i]] for i in range(len(layer_mu))]
                        cumulative_tuning_mu = sum(tuning_mu)

                        if n_log_spots == n_plan_spots:
                            tuning_xy_tuples = [np.array((tuning_xy[i], tuning_xy[i + 1])) for i in range(len(tuning_xy)) if i % 2 == 0]
                            for t, tuning_spot in enumerate(tuning_xy_tuples):
                                shifts = [tuning_spot - plan_spot for plan_spot in plan_xy_tuples]
                                dists = [shift.dot(shift) for shift in shifts]
                                closest = dists.index(min(dists))
                                mu_sorted[closest] += tuning_mu[t]

                        if n_log_spots < n_plan_spots:  # tuning spot replaces one planned spot, map it to closest plan spot
                            mu_sorted.insert(tuning_dict[layer_id], cumulative_tuning_mu)
                        
                        if n_log_spots > n_plan_spots:
                            print(f'  /!\ Critical: More spots recorded than planned in layer-ID {layer_id} ({n_log_spots} vs. {n_plan_spots}), skipping this..')
                            continue

                        plan_beam.IonControlPointSequence[layer_id * 2].CumulativeMetersetWeight = cumulative_mu
                        cumulative_mu += sum(mu_sorted)
                        plan_beam.IonControlPointSequence[layer_id * 2 + 1].CumulativeMetersetWeight = cumulative_mu
                        plan_beam.IonControlPointSequence[layer_id * 2].ScanSpotMetersetWeights = mu_sorted
                        plan_beam.IonControlPointSequence[layer_id * 2 + 1].ScanSpotMetersetWeights = [0.0 for _ in range(len(mu_sorted))]
                    
                    else:
                        print(f'/x\ Mode {mode} not supported, exiting..')
                        return None
                
                if mode == 'all' or mode == 'mu':
                    ds.FractionGroupSequence[0].ReferencedBeamSequence[i].BeamMeterset = cumulative_mu
                    plan_beam.FinalCumulativeMetersetWeight = cumulative_mu
                
                plan_beam.NumberOfControlPoints = total_layers * 2
                
                if len(plan_beam.IonControlPointSequence) > 2 * total_layers:
                    diff = len(plan_beam.IonControlPointSequence) - 2 * total_layers
                    print(f'  /!\ Deleting {diff} entries from beam {beam_id} control point sequence..')
                    plan_beam.IonControlPointSequence = plan_beam.IonControlPointSequence[: - diff]
                
                ds.IonBeamSequence[i] = plan_beam
            
            ds.FractionGroupSequence[0].NumberOfFractionsPlanned = 1
            sop_instance = str(ds.SOPInstanceUID).split('.')
            sop_instance[-1] = '99999'
            new_sop_instance = '.'.join(sop_instance)
            ds.SOPInstanceUID = new_sop_instance
            
            if not ds.RTPlanName.__contains__('log'):
                ds.RTPlanLabel = f'(RPR)-{fx_id}-' + ds.RTPlanName.strip('_unapprove')
                ds.RTPlanName = ds.RTPlanName.strip('_unapprove') + f'_LOG_{fx_id}_{mode}'
            ds.StudyDate = pd.to_datetime(fx_id).date()

            # ds.ReferringPhysicianName = 'Wolter^Lukas'
            ds.ApprovalStatus = 'UNAPPROVED'

            print('\nWriting dicom..')
            destination = os.path.join(dcm_path, f'RP{ds.SOPInstanceUID}_fx_{fx_id}_log_{mode}.dcm')
            pydicom.write_file(destination, ds)
            print(f'Wrote log-based plan to RP{ds.SOPInstanceUID}_fx_{fx_id}_log_{mode}.dcm')
                

    def fractional_evolution(self, all=False):
        # bs_file = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\cct_entities.csv'
        # bs_data = pd.read_csv(bs_file, index_col='PATIENT_ID')
        # bs_dict = bs_data.to_dict()
        # bodysites = bs_data['BODY_SITE'].drop_duplicates().to_list()  # to be continued
        gcf_update = pd.to_datetime('20200510').date()
        ic1_change = pd.to_datetime('20200713').date()

        other_records, mobil_dfs = [], r'N:\fs4-HPRT\HPRT-Docs\Lukas\Logfile_Extraction\dataframes'
        for file in sorted(os.listdir(self.df_destination)):
            if file.__contains__('records') and file.endswith('.csv') and not file.__contains__('QA'):
                other_records.append(os.path.join(self.df_destination, file))
        for file in sorted(os.listdir(mobil_dfs)):
            if file.__contains__('records') and file.endswith('.csv') and not file.__contains__('QA'):
                other_records.append(os.path.join(mobil_dfs, file))

        fig = plt.figure(1, figsize=(10, 4))
        global_dist_means = []
        for n, record_file in enumerate(other_records):
            if not all:
                print(f'Starting patient {self.patient_id}..')
                this_record_df = self.patient_record_df
            else:
                print(f'Starting record dataframe ({n + 1}/{len(other_records)}) {record_file}..')
                this_record_df = pd.read_csv(record_file, index_col='TIME', dtype={'BEAM_ID':str, 'FRACTION_ID':str})         
            
            beam_list = this_record_df['BEAM_ID'].drop_duplicates()
            patient_id = this_record_df['PATIENT_ID'].iloc[0]
            # if patient_id != 1617814: continue
            # bodysite = bs_dict[patient_id]
            for beam_id in beam_list:
                # print(f'  Processing beam-ID {beam_id}..')
                beam_df = this_record_df.loc[this_record_df['BEAM_ID'] == beam_id]
                beam_fxs = beam_df['FRACTION_ID'].drop_duplicates()
                ref_df = beam_df.loc[beam_df['FRACTION_ID'] == beam_fxs.iloc[0]]
                date_axis, fx_dist_means = [], []
                for fx_id in beam_fxs[1:]:
                    ref_df_copy = ref_df.copy()
                    fx_df = beam_df.loc[beam_df['FRACTION_ID'] == fx_id]
                    while len(fx_df) != len(ref_df_copy):
                        # print(f'    Correcting fx-ID {fx_id}', end='\r')
                        dropped = False
                        for index, pair in enumerate(zip(fx_df['X_POSITION(mm)'], ref_df_copy['X_POSITION(mm)'])):
                            if abs(pair[0] - pair[1]) > 2:  # mm difference
                                if len(fx_df) > len(ref_df_copy):
                                    fx_df = fx_df.drop(fx_df.index[index])
                                else:
                                    ref_df_copy = ref_df_copy.drop(ref_df_copy.index[index])
                                
                                dropped = True
                                break
                        
                        if not dropped:
                            print(f'  /!\ Slicing dataframe due to shape mismatch (beam-ID {beam_id}, fx-ID {fx_id}, delta {abs(len(fx_df) - len(ref_df_copy))} spot)..')
                            if len(fx_df) > len(ref_df_copy):
                                fx_df = fx_df.loc[:fx_df.index[len(ref_df_copy) - 1]]
                            else:
                                ref_df_copy = ref_df_copy.loc[:ref_df_copy.index[len(fx_df) - 1]]

                    date = pd.to_datetime(fx_df.index[0], dayfirst=False).date()
                    dx = fx_df['X_POSITION(mm)'].to_numpy() - ref_df_copy['X_POSITION(mm)'].to_numpy()
                    dy = fx_df['Y_POSITION(mm)'].to_numpy() - ref_df_copy['Y_POSITION(mm)'].to_numpy()
                    shift_vectors = [np.array(tup) for tup in zip(dx, dy)]
                    mean_dist = np.linalg.norm(sum(shift_vectors) / len(shift_vectors))
                    if mean_dist > 2:
                        print('/x\\', mean_dist, patient_id, fx_id, date, beam_id)
                        # shifts = [np.linalg.norm(shift) for shift in shift_vectors]
                        # max_shift = max(shifts)
                        # max_index = shifts.index(max_shift)
                        # display = ref_df['LAYER_ID'].iloc[max_index]
                        # plt.scatter(ref_df.loc[ref_df['LAYER_ID'] == display, ['X_POSITION(mm)']], ref_df.loc[ref_df['LAYER_ID'] == display, ['Y_POSITION(mm)']], marker='o', edgecolors='black', facecolors='white', label=f'{beam_fxs.iloc[0]} (first)', zorder=5)
                        # plt.scatter(fx_df.loc[fx_df['LAYER_ID'] == display, ['X_POSITION(mm)']], fx_df.loc[fx_df['LAYER_ID'] == display, ['Y_POSITION(mm)']], marker='o', label=f'{fx_id} (this)', c='red')
                        # plt.scatter(beam_df.loc[(beam_df['LAYER_ID'] == display) & (beam_df['FRACTION_ID'] != fx_id), ['X_POSITION(mm)']], beam_df.loc[(beam_df['LAYER_ID'] == display) & (beam_df['FRACTION_ID'] != fx_id), ['Y_POSITION(mm)']], marker='o', alpha=0.3, label=f'remaining', c='tab:blue', zorder=-5)
                        # plt.title(f'Pat. {patient_id}, Beam {beam_id}, (showing Layer-iD {display}/{fx_df["TOTAL_LAYERS"].iloc[0]})\nEntire field mean pos. $\Delta$ = {mean_dist:.3f}')
                        # plt.legend()
                        # plt.savefig(f'{output_dir}/{patient_id}_{fx_id}_{beam_id}_failure.png', dpi=300)
                        # plt.show()
                        # plt.clf()
                        # continue

                    fx_dist_means.append(mean_dist), date_axis.append(date)
                    global_dist_means.append(mean_dist)

                if not all:
                    plt.errorbar(x=sorted(date_axis), y=fx_dist_means, yerr=None, fmt='o-', capsize=3, label=beam_id)
                    pass
                else:
                    plt.errorbar(x=date_axis, y=fx_dist_means, yerr=None, fmt='o', capsize=3,  color='black', alpha=0.3, markersize=4)
                    pass

            if not all:
                break    
            else:
                # print(f'\n>> FINISHED', n + 1, 'of', len(other_records), 'datasets <<\n')  
                pass  
        
        passed, passmark = 0, 0.5
        for dist in global_dist_means:
            if dist < passmark:
                passed += 1
        
        percentage = passed / len(global_dist_means) * 100
        
        plt.xlabel('Date [YYYY-MM-DD]')
        # plt.ylabel('Mean field $\Delta$ to Fx-01 [mm]')
        plt.ylabel('Mean spot position offset [mm]')
        plt.ylim(0.0, 5.0)
        plt.grid(axis='y')
        if not all:
            plt.title(f'Delivery variance in recorded spot position (pat.-ID {self.patient_id})', fontweight='bold')
            plt.legend(title='Beam-ID')
            plt.tight_layout()
            plt.savefig(f'{output_dir}/{self.patient_id}_fractional_fluctuation.png', dpi=1000)
        else:
            # plt.title(f'Delivery variance in recorded spot position | N={len(other_records)} patients, 04/2020 - present', fontweight='bold')
            plt.axvline(gcf_update, color='tab:blue', zorder=-10)
            plt.axvline(ic1_change, color='tab:blue', zorder=-10)
            plt.legend(title=f'{round(percentage, 2)}% of fields within $\pm${passmark}mm')
            plt.tight_layout()
            plt.savefig(f'{output_dir}/full_fractional_fluctuation.png', dpi=1000)
        plt.show()    
    

    def beam_histos(self):
        for file in sorted(os.listdir(self.df_destination)):
            if file.__contains__(str(self.patient_id)) and file.__contains__('delta') and file.endswith('.csv'):
                print(f'''Found patient deltaframe '{file}', reading in..''')
                self.patient_delta_df = pd.read_csv(os.path.join(self.df_destination, file), index_col='UNIQUE_INDEX', dtype={'BEAM_ID':str, 'FRACTION_ID':str})
                break 

        beam_list = self.patient_record_df['BEAM_ID'].drop_duplicates()
        fig, axs = plt.subplots(2, len(beam_list), figsize=(50, 10), dpi=80)
        ax0 = fig.add_subplot(111, frameon=False)
        ax0.set_xticks([]), ax0.set_yticks([])
        for i, beam_id in enumerate(beam_list):
            beam_df = self.patient_record_df.loc[self.patient_record_df['BEAM_ID'] == beam_id]
            delta_df = self.patient_delta_df.loc[self.patient_delta_df['BEAM_ID'] == beam_id]
            delta_df.index = beam_df.index
            joint_df = pd.concat([beam_df, delta_df], axis=1)
            joint_df.dropna(inplace=True)
            plan_x = joint_df['DELTA_X(mm)'].to_numpy() + joint_df['X_POSITION(mm)'].to_numpy()
            plan_y = joint_df['DELTA_Y(mm)'].to_numpy() + joint_df['Y_POSITION(mm)'].to_numpy()
            
            if len(beam_df) != len(delta_df):
                print(len(beam_df), len(delta_df))
                continue
            
            bins = 80
            axs[0, i].hist(joint_df['DELTA_X(mm)'], bins=bins, label='log-file X')
            axs[0, i].hist(plan_x - joint_df['X_POSITION_CORR(mm)'].to_numpy(), bins=bins, label='log-file X\n(corrected for iso shift)', alpha=0.5)
            axs[0, i].set_xlim(-2, 2)
            axs[0, i].annotate(f'{beam_df.GANTRY_ANGLE.mean():.2f}°', xy=(0.05, 0.90), xycoords='axes fraction')
            axs[1, i].hist(joint_df['DELTA_Y(mm)'], bins=bins, label='log-file Y')
            axs[1, i].hist(plan_y - joint_df['Y_POSITION_CORR(mm)'].to_numpy(), bins=bins, label='log-file Y\n(corrected for iso shift)', alpha=0.5)
            axs[1, i].set_xlim(-2, 2)
            axs[1, i].annotate(f'{beam_df.GANTRY_ANGLE.mean():.2f}°', xy=(0.05, 0.90), xycoords='axes fraction')
        
        axs[0, -1].legend()
        axs[1, -1].legend()
        plt.title(f'Positional difference to plan (patient-ID {self.patient_id})', fontweight='bold')
        plt.savefig(f'{output_dir}/{self.patient_id}_gtr_corrected_hist.png', dpi=600)
        plt.show()
        plt.clf()

        # for i, beam_id in enumerate(beam_list):
        #     beam_df = self.patient_record_df.loc[self.patient_record_df['BEAM_ID'] == beam_id]
        #     delta_df = self.patient_delta_df.loc[self.patient_delta_df['BEAM_ID'] == beam_id]
        #     delta_df.index = beam_df.index
        #     diff_cols = delta_df.columns.difference(beam_df.columns)
        #     joint_df = pd.concat([beam_df, delta_df[diff_cols]], axis=1)
        #     joint_df.dropna(inplace=True)

        #     if joint_df.empty or i != 3:
        #         continue
            
        #     plan_x = joint_df['DELTA_X(mm)'].to_numpy() + joint_df['X_POSITION(mm)'].to_numpy()
        #     plan_y = joint_df['DELTA_Y(mm)'].to_numpy() + joint_df['Y_POSITION(mm)'].to_numpy()
        #     joint_df['CORR_DELTA_X'] = plan_x - beam_df['X_POSITION_CORR(mm)'].to_numpy()
        #     joint_df['CORR_DELTA_Y'] = plan_y - beam_df['Y_POSITION_CORR(mm)'].to_numpy()

            # x_deltas = {}
            # for i, fx in enumerate(joint_df.FRACTION_ID.drop_duplicates().to_list()):
            #     print(fx, type(fx))
            #     fx_df = joint_df.loc[joint_df['FRACTION_ID'] == fx]
            #     x_deltas[fx] = fx_df['DELTA_X(mm)'].to_list()
            #     temp_df = pd.DataFrame(x_deltas)
            #     sns.histplot(temp_df, kde=True, cbar_kws={'color':'None'})
            #     plt.show()


            # sns.histplot(joint_df[['DELTA_X(mm)', 'CORR_DELTA_X']])
            # print(joint_df.columns)
            # hue = 'Y_POSITION(mm)'
            # sns.pairplot(joint_df.loc[joint_df['FRACTION_ID'] == self.fraction_list[1]][['DELTA_X(mm)', 'DELTA_Y(mm)', hue]], vars=['DELTA_X(mm)', 'DELTA_Y(mm)'], hue=hue, corner=True)
            # plt.show()

    def prepare_sss_dataframe(self):
        # locate existing dataframes
        delta_exists, single_stats_exists = False, False
        for file in sorted(os.listdir(self.df_destination)):
            if file.__contains__(str(self.patient_id)) and file.__contains__('delta') and file.endswith('.csv'):
                print(f'''Found patient deltaframe '{file}', reading in..''')
                self.patient_delta_df = pd.read_csv(os.path.join(self.df_destination, file), index_col='UNIQUE_INDEX', dtype={'BEAM_ID':str, 'FRACTION_ID':str})
                delta_exists = True
            elif file.__contains__(str(self.patient_id)) and file.__contains__('sss_data') and file.endswith('.csv'):
                print(f'''Found patient spot statistics data '{file}', reading in..''')
                self.patient_sss_df = pd.read_csv(os.path.join(self.df_destination, file), dtype={'BEAM_ID':str})
                single_stats_exists = True
            if delta_exists and single_stats_exists:
                re_init = input(f'''Re-initialize [y/N]? ''')
                if re_init == 'y':
                    break
                else:
                    return None
        
        if not delta_exists:
            print(f'  /x\ No deltaframe detected for patient_id {self.patient_id}, please run prepare_deltaframe()')
        
        # get spot counts per beam and fraction
        print(f'\nInitializing spot statistics dataframe for patient-ID {self.patient_id}..')
        self.patient_delta_df['DIST(mm)'] = np.sqrt(self.patient_delta_df['DELTA_X(mm)'] ** 2 + self.patient_delta_df['DELTA_Y(mm)'] ** 2)
        beams = self.patient_delta_df.BEAM_ID.drop_duplicates()
        spots_data = {beam_id:[] for beam_id in beams}
        operation_columns = ['DELTA_X(mm)', 'DELTA_Y(mm)', 'DELTA_MU', 'DELTA_E(MeV)']

        for fx_id in self.fraction_list:
            for beam_id in beams:
                beam_df = self.patient_delta_df.loc[(self.patient_delta_df['FRACTION_ID'] == fx_id) & (self.patient_delta_df['BEAM_ID'] == beam_id)]
                if beam_df[operation_columns].isnull().values.all():
                    spots_data[beam_id].append(0)
                else:
                    spots_data[beam_id].append(len(beam_df))
        
        operation_columns += ['X_WIDTH(mm)', 'Y_WIDTH(mm)', 'DIST(mm)']
        to_concat = []
        for beam_id in spots_data.keys():
            beam_df = self.patient_delta_df.loc[self.patient_delta_df['BEAM_ID'] == beam_id].reset_index()
            beam_rec = self.patient_record_df.loc[self.patient_record_df['BEAM_ID'] == beam_id].reset_index()
            if len(beam_df) == len(beam_rec):
                try:
                    beam_df = pd.concat([beam_df, beam_rec[['X_WIDTH(mm)', 'Y_WIDTH(mm)', 'X_POSITION(mm)', 'Y_POSITION(mm)', 'MU', 'BEAM_CURRENT(V)', 'DRILL_TIME(ms)', 'LAYER_ENERGY(MeV)']]], axis=1)
                except KeyError:
                    beam_df = pd.concat([beam_df, beam_rec[['X_WIDTH(mm)', 'Y_WIDTH(mm)', 'X_POSITION(mm)', 'Y_POSITION(mm)', 'MU', 'BEAM_CURRENT(A)', 'DRILL_TIME(ms)', 'LAYER_ENERGY(MeV)']]], axis=1)
            else:
                print('\n  /!\ Record DF doesnt match Delta DF! \n')
                return None
            
            # accept only beams with >= 15 fx
            valid_counts = [count for count in spots_data[beam_id] if count != 0]
            if len(valid_counts) < 15:
                continue

            # set reference fraction (first fx holding the most common spot count)
            most_common = max(set(valid_counts), key=valid_counts.count)
            reference_fx = self.fraction_list[spots_data[beam_id].index(most_common)]
            print(f'  Referencing all spots for beam-ID {beam_id} to fx-ID {reference_fx}..')
            reference_df = beam_df.loc[beam_df['FRACTION_ID'] == reference_fx].reset_index(drop=True)

            # fractional evolution by layer
            beam_sss_df = pd.DataFrame()
            for layer_id in reference_df.LAYER_ID.drop_duplicates():
                print(f'    Starting layer-ID {layer_id}')
                layer_record = self.patient_record_df.loc[(self.patient_record_df.BEAM_ID == beam_id) & (self.patient_record_df.LAYER_ID == layer_id)] 
                fx_dfs = []
                for fx_id in self.fraction_list:
                    fx_df = beam_df.loc[(beam_df.FRACTION_ID == fx_id) & (beam_df.LAYER_ID == layer_id)]
                    fx_df = fx_df.set_index('SPOT_ID')[operation_columns]
                    layer_reference = reference_df.loc[reference_df.LAYER_ID == layer_id].copy()  # keep this here, might be modified

                    if fx_df.empty: 
                        continue
                    if fx_id == reference_fx:
                        fx_dfs.append(fx_df)
                        continue  
                    
                    record_reference, record_current = layer_record.loc[layer_record.FRACTION_ID == reference_fx].set_index('SPOT_ID').copy(), layer_record.loc[layer_record.FRACTION_ID == fx_id].set_index('SPOT_ID').copy()
                    diff = abs(record_reference[['X_POSITION(mm)', 'Y_POSITION(mm)']] - record_current[['X_POSITION(mm)', 'Y_POSITION(mm)']])
                    if len(fx_df) == len(layer_reference):
                        if (diff > 3).any(axis=None):
                            print(f'    /!\ X/Y fluctuation max of {np.round(diff.max().max(), 2)}mm in fraction {fx_id}, handling..')
                        
                        # verify correct sorting, cases exist where two different spots are missing in each df, but len() is equal
                        while (diff > 3).any(axis=None):
                            xy_ref, xy_cur = pd.DataFrame(), pd.DataFrame()
                            xy_ref[['X_REF', 'Y_REF']] = record_reference[['X_POSITION(mm)', 'Y_POSITION(mm)']]
                            xy_cur[['X_CUR', 'Y_CUR']] = record_current[['X_POSITION(mm)', 'Y_POSITION(mm)']]
                            comp = pd.concat([xy_ref, xy_cur], axis=1)
                            err = comp.loc[(abs(comp.X_REF - comp.X_CUR) > 3) | (abs(comp.Y_REF - comp.Y_CUR) > 3)].index[0]
                            if abs(xy_cur.X_CUR.iloc[err] - xy_ref.X_REF[err + 1]) < abs(xy_cur.X_CUR.iloc[err + 1] - xy_ref.X_REF[err]):  # current is missing spot
                                print('CUR')
                                record_current.index = record_current.index * 2 + 1
                                record_current.loc[err * 2] = np.nan
                                record_current = record_current.sort_index().reset_index(drop=True)
                                record_current.index.name = 'SPOT_ID'
                                fx_df.index = fx_df.index * 2 + 1
                                fx_df.loc[err * 2] = np.nan
                                fx_df = fx_df.sort_index().reset_index(drop=True)
                                fx_df.index.name = 'SPOT_ID'
                            else:
                                print('REF')
                                record_reference.index = record_reference.index * 2 + 1
                                record_reference.loc[err * 2] = np.nan
                                record_reference = record_reference.sort_index().reset_index(drop=True)
                                record_reference.index.name = 'SPOT_ID'
                                layer_reference.index = layer_reference.index * 2 + 1
                                layer_reference.loc[err * 2] = np.nan
                                layer_reference = layer_reference.sort_index().reset_index(drop=True)
                                layer_reference.index.name = 'SPOT_ID'

                            diff = abs(record_reference[['X_POSITION(mm)', 'Y_POSITION(mm)']] - record_current[['X_POSITION(mm)', 'Y_POSITION(mm)']])

                        if (diff > 1).any(axis=None):
                            print(f'    /!\ X/Y fluctuation max of {np.round(diff.max().max(), 2)}mm in fraction {fx_id}')

                    # need to drop certain values
                    else:
                        print(f'    /!\ Length mismatch in fraction {fx_id}, processing..')
                        record_reference, record_current = layer_record.loc[layer_record.FRACTION_ID == reference_fx].copy(), layer_record.loc[layer_record.FRACTION_ID == fx_id].copy()
                        record_reference['X_REF'], record_reference['Y_REF'] = record_reference['X_POSITION(mm)'].round(0), record_reference['Y_POSITION(mm)'].round(0)
                        record_current['X_CURRENT'], record_current['Y_CURRENT'] = record_current['X_POSITION(mm)'].round(0), record_current['Y_POSITION(mm)'].round(0)
                        record_reference.set_index('SPOT_ID', inplace=True), record_current.set_index('SPOT_ID', inplace=True)
                        while len(record_current) != len(record_reference):
                            comp = pd.concat([record_reference, record_current[['X_CURRENT', 'Y_CURRENT']]], axis=1)
                            indices = comp.loc[(abs(comp['X_REF'] - comp['X_CURRENT']) > 1) | (abs(comp['Y_REF'] - comp['Y_CURRENT']) > 1) | (comp['X_CURRENT'].apply(np.isnan)) | (comp['X_REF'].apply(np.isnan))].index
                            # if current fx has to many spots, drop them
                            if len(record_current) > len(record_reference):
                                record_current.drop(record_current.index[indices[0]], inplace=True)
                                fx_df.drop(fx_df.index[indices[0]], inplace=True)
                            # if current fx has less spots, fill nan
                            else:
                                record_reference.drop(record_reference.index[indices[0]], inplace=True)
                                fx_df.index = fx_df.index * 2 + 2
                                fx_df.loc[indices[0] * 2 + 1] = np.nan

                            fx_df = fx_df.sort_index().reset_index(drop=True)
                            record_reference.reset_index(drop=True, inplace=True), record_current.reset_index(drop=True, inplace=True)
                        
                        fx_df.index.name = 'SPOT_ID'
                    
                    fx_dfs.append(fx_df)
                
                layer_stacked = pd.concat(fx_dfs)
                layer_means, layer_stds = layer_stacked.groupby('SPOT_ID').mean(), layer_stacked.groupby('SPOT_ID').std()
                layer_means.rename(columns={col:col + '_MEAN' for col in operation_columns}, inplace=True)
                layer_stds.rename(columns={col:col + '_STD' for col in operation_columns}, inplace=True)
                sss_data = pd.concat([layer_means, layer_stds], axis=1)
                beam_sss_df = pd.concat([beam_sss_df, sss_data])

            beam_sss_df.reset_index(drop=True, inplace=True)
            reference_df.drop(columns=operation_columns, inplace=True)
            finalized = pd.concat([reference_df, beam_sss_df], axis=1)
            to_concat.append(finalized)
        
        self.patient_sss_df = pd.concat(to_concat, ignore_index=True)
        self.patient_sss_df = self.patient_sss_df.loc[:,~self.patient_sss_df.columns.duplicated()]

        # write deltaframe
        os.chdir(self.df_destination)
        print(f'''  ..Writing spot statistics data to '{self.df_destination}' as .CSV.. ''')
        while True:   
            try:
                self.patient_sss_df.to_csv(f'patient_{self.patient_id}_sss_data.csv')
                break
            except PermissionError:
                input('  Permission denied, close target file and press ENTER.. ')
        print('Complete')
                

    def sss_histograms(self, mode='pos'):
        # locate existing dataframe
        single_stats_exists = False
        for file in sorted(os.listdir(self.df_destination)):
            if file.__contains__(str(self.patient_id)) and file.__contains__('sss_data') and file.endswith('.csv'):
                print(f'''Found patient spot statistics data '{file}', reading in..''')
                self.patient_sss_df = pd.read_csv(os.path.join(self.df_destination, file), dtype={'BEAM_ID':str})
                single_stats_exists = True

        if not single_stats_exists:
            print(f'/x\ No single spot statistics data found for patient-ID {self.patient_id}, exiting..')
            return None
        
        beams = self.patient_sss_df.dropna().BEAM_ID.drop_duplicates().to_list()

        # histograms (one patient)
        sns.set()
        cmap = sns.color_palette()
        
        if mode == 'pos':
            fig, axs = plt.subplots(len(beams), 3, figsize=(20, 10 * len(beams) / 3), sharex='col', sharey='row')
            fig.suptitle(f'Patient {self.patient_id} Delivery Report', fontweight='bold')
            axs[0, 0].set_title('X-position', fontstyle='italic')
            axs[0, 1].set_title('Y-position', fontstyle='italic')
            axs[0, 2].set_title('Dose', fontstyle='italic')
            target_cols = [['DELTA_X(mm)_MEAN', 'DELTA_X(mm)_STD'], ['DELTA_Y(mm)_MEAN', 'DELTA_Y(mm)_STD'], ['DELTA_MU_MEAN', 'DELTA_MU_STD']]
            ax_labels = ['$\Delta X$ [mm]', '$\Delta Y$ [mm]', '$\Delta D$ [MU]']
            ax_units = ['mm', 'mm', 'MU']
            label_text = ['Accuracy ($\mu$):', 'Precision ($\sigma$):']
            for nr, bid in enumerate(beams):
                df = self.patient_sss_df.loc[self.patient_sss_df.BEAM_ID == bid]
                for i, (ax, cols) in enumerate(zip(axs[nr], target_cols)):
                    if i != 2:  # positions
                        sns.histplot(df[cols], kde=True, ax=ax, binwidth=0.02)
                        ax.set_xlim(-1.5, 1.5)
                    else:  # dose
                        sns.histplot(df[cols], kde=True, ax=ax)
                        ax.set_xlim(-0.02, 0.02)
                    for j, col in enumerate(cols):
                        mean, std = df[col].mean(), df[col].std()
                        ax.axvline(mean, ls='--', lw=1, color=cmap[j], label=f'{label_text[j]} {np.round(mean, 3)} $\pm$ {np.round(std, 3)} {ax_units[i]}')
                    
                    ax.set_ylabel('')
                    ax.set_xlabel(ax_labels[i])
                    ax.legend(title='Single spot statistics', loc='upper left')
                    ax.text(0.98, .92, f'{df.GANTRY_ANGLE.iloc[0]}°', horizontalalignment='right', verticalalignment='center', transform = ax.transAxes, fontweight='bold', fontsize=14)

            out = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\Logfiles\Reports'
            plt.tight_layout()
            plt.savefig(os.path.join(out, f'{id}_sss_hist.png'), dpi=300)
            plt.show()
        
        elif mode == 'width':
            import matplotlib.colors as mcolors
            from matplotlib.patches import Patch
            from matplotlib.lines import Line2D
            from scipy import interpolate
    
            beam_model = pd.read_csv(r'N:\fs4-HPRT\HPRT-Docs\Lukas\BA_Anna_Leimbach\Beam_Modell\beam_sigma_RS.csv')
            bm_energies = beam_model.iloc[:, 0]
            bm_fwhm = beam_model.iloc[:, 1] * 10 * 2.355
            bm_interp = interpolate.interp1d(bm_energies, bm_fwhm, kind='cubic')
            cmaps, colors = [plt.get_cmap('Reds_r'), plt.get_cmap('Blues_r')], []
            values = np.linspace(0, 1, 10)
            for cmap in cmaps:
                colors.append([mcolors.to_hex(cmap(value)) for value in values])

            markers = ['s', 'o', 'v', '^']
            legend_elements = [Patch(facecolor=colors[0][1], edgecolor='None', label='$x$-FWHM'), 
                               Patch(facecolor=colors[1][3], edgecolor='None', label='$y$-FWHM')]
            
            fig1, axs1 = plt.subplots(1, 1, figsize=(6, 6))
            fig2, axs2 = plt.subplots(1, 1, figsize=(6, 6))
            for nr, bid in enumerate(beams):
                df = self.patient_record_df.loc[self.patient_record_df.BEAM_ID == bid]
                energies = sorted(df['LAYER_ENERGY(MeV)'].drop_duplicates().to_list())
                gantry = df.GANTRY_ANGLE.iloc[0]
                x_fwhm_mean, x_fwhm_std = [], []
                y_fwhm_mean, y_fwhm_std = [], []
                delta_x_rel, delta_y_rel = [], []
                for energy in energies:
                    e_df = df.loc[df['LAYER_ENERGY(MeV)'] == energy]
                    x_fwhm_mean.append(e_df['X_WIDTH(mm)'].mean()), x_fwhm_std.append(e_df['X_WIDTH(mm)'].std())
                    y_fwhm_mean.append(e_df['Y_WIDTH(mm)'].mean()), y_fwhm_std.append(e_df['Y_WIDTH(mm)'].std())
                    delta_x_rel.append((e_df['X_WIDTH(mm)'].mean() - bm_interp(energy)) / bm_interp(energy) * 100)
                    delta_y_rel.append((e_df['Y_WIDTH(mm)'].mean() - bm_interp(energy)) / bm_interp(energy) * 100)

                # axs1.errorbar(x=energies, y=y_fwhm_mean, marker=markers[nr], yerr=y_fwhm_std, color=colors[0][nr + 1])  # x <-> y
                # axs1.errorbar(x=energies, y=x_fwhm_mean, marker=markers[nr], yerr=x_fwhm_std, color=colors[1][nr + 3])
                axs1.plot(energies, delta_x_rel, color='tab:orange', label='X')  # relative diff to bm
                axs1.plot(energies, delta_y_rel, color='tab:blue', label='Y')

                legend_elements.append(Line2D([0], [0], color='black', marker=markers[nr], linestyle='None', label=f'Beam {nr + 1}: {gantry}°'))

            legend_elements.append(Line2D([0], [0], color='black', linestyle='--', label='Beam-Modell'))
            axs1.plot(bm_energies, bm_fwhm, '--', color='black', zorder=-10)
            # axs1.set_ylabel('FWHM [mm]')
            axs1.set_xlabel('Energie [MeV]')
            axs1.set_ylabel('Rel. Diff. zum BM [%]')
            axs2.set_ylabel('Anzahl')
            axs2.set_xlabel('$\sigma_{FWHM}$ [mm]')
            # axs[1].set_xlim(0, 0.25)
            # axs1.legend(handles=legend_elements, title='Energie-/Winkelabhängigkeit')
            # plt.show()

            # fig.suptitle(f'Patient {log.patient_id} Delivery Report', fontweight='bold')
            target_cols = ['Y_WIDTH(mm)_STD', 'X_WIDTH(mm)_STD']  # x <-> y
            labels = ['$x$-FWHM', '$y$-FWHM']
            for i, cols in enumerate(target_cols):
                sns.histplot(self.patient_sss_df[cols], color=colors[i][1], kde=True, ax=axs2, label=labels[i], binwidth=0.0025, element='step')

            out = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\Logfiles\Reports'
            # plt.tight_layout()
            axs2.legend(title='Einzelspot-Reprod.')
            axs2.set_xlim(0, 0.3)
            fig1.tight_layout(), fig2.tight_layout()
            fig1.savefig(os.path.join(out, f'{id}_fwhm_vs_energy.png'), dpi=300)
            fig2.savefig(os.path.join(out, f'{id}_histograms.png'), dpi=300)
            # plt.show()


    def sss_boxplot(self):
        from matplotlib.patches import Patch, Rectangle
        from matplotlib.lines import Line2D
        from matplotlib.ticker import FormatStrFormatter
        out = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\Logfiles\Reports'
        ponaqua_qualified = [id.strip('\n') for id in open(r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\qualified_IDs.txt', 'r').readlines()]

        # boxplots (all patient with sss data)
        sss_files = [file for file in os.listdir(self.df_destination) if file.__contains__('_sss_') and file.endswith('.csv')]
        rec_files = [file for file in os.listdir(self.df_destination) if file.__contains__('_delta') and file.endswith('.csv')]
        patient_dfs = []
        for sss_file in sss_files:
            print(f'  Reading {sss_file}..')
            id = int(sss_file.split('_')[1])
            sss_df = pd.read_csv(os.path.join(self.df_destination, sss_file))
            sss_df['PATIENT_ID'] = id
            patient_dfs.append(sss_df)
        
        print('Calculating record worst cases:')
        worst_x, worst_y, worst_mu = 0, 0, 0
        pat_x, pat_y, pat_mu = None, None, None
        for rec_file in rec_files:
            id = rec_file.split('_')[1]
            rec_df = pd.read_csv(os.path.join(self.df_destination, rec_file))
            try: rec_df.drop(rec_df.loc[rec_df.BEAM_ID.astype(int) > 4].index, inplace=True)
            except: pass
            min_x, min_y = rec_df['DELTA_X(mm)'].min() , rec_df['DELTA_Y(mm)'].min()
            max_x, max_y = rec_df['DELTA_X(mm)'].max(), rec_df['DELTA_Y(mm)'].max()
            min_mu, max_mu = rec_df['DELTA_MU'].min() , rec_df['DELTA_MU'].max()

            if abs(worst_x) < max(np.abs([min_x, max_x])):
                if abs(min_x) > abs(max_x):
                    worst_x = min_x
                else:
                    worst_x = max_x
                pat_x = id

            if abs(worst_y) < max(np.abs([min_y, max_y])):
                if abs(min_y) > abs(max_y):
                    worst_y = min_y
                else:
                    worst_y = max_y
                pat_y = id
            
            if abs(worst_mu) < max(np.abs([min_mu, max_mu])):
                if abs(min_mu) > abs(max_mu):
                    worst_mu = min_mu
                else:
                    worst_mu = max_mu
                pat_mu = id
        
        # print('  X: ', worst_x, 'for patient', pat_x)
        # print('  Y: ', worst_y, 'for patient', pat_y)
        # print('  MU:', worst_mu, 'for patient', pat_mu)
        
        df, n_fields = pd.concat(patient_dfs), 0
        df.reset_index(inplace=True)
        # over_twomm = df.loc[df['DIST(mm)_MEAN'] > 2.0]
        # print(len(df.loc[df['DIST(mm)_MEAN'] > 2.0]), len(df), len(df.loc[df['DIST(mm)_MEAN'] > 2.0])/len(df) * 100)
        # for layer_id in over_twomm.LAYER_ID.drop_duplicates():
        #     plot_df = over_twomm.loc[over_twomm.LAYER_ID == layer_id]
        #     previous_fx = df.loc[(df.PATIENT_ID == plot_df.PATIENT_ID.iloc[0]) & (df.FRACTION_ID.astype(int) == plot_df.FRACTION_ID.iloc[0] - 1) & (df.BEAM_ID == plot_df.BEAM_ID.iloc[0]) & (df.LAYER_ID == plot_df.LAYER_ID.iloc[0])]
        #     plt.plot(previous_fx['X_POSITION(mm)'], previous_fx['Y_POSITION(mm)'], 'x-')
        #     plt.plot(plot_df['X_POSITION(mm)'], plot_df['Y_POSITION(mm)'], 'o-')
        #     plt.show()
        # return None
        # print(len(df.GANTRY_ANGLE.dropna().drop_duplicates()), df.GANTRY_ANGLE.dropna().drop_duplicates())
        for id in df.PATIENT_ID.drop_duplicates():
            n_fields += len(df.loc[df.PATIENT_ID == id, 'BEAM_ID'].drop_duplicates())
        
        for gtr in df.GANTRY_ANGLE.drop_duplicates():
            gtr_df = df.loc[df.GANTRY_ANGLE == gtr]
            df.loc[df.GANTRY_ANGLE == gtr, 'MEDIAN_X_STD'] = gtr_df['DELTA_X(mm)_STD'].median()
            df.loc[df.GANTRY_ANGLE == gtr, 'MEDIAN_Y_STD'] = gtr_df['DELTA_Y(mm)_STD'].median()
            df.loc[df.GANTRY_ANGLE == gtr, 'MEDIAN_MU_STD'] = gtr_df['DELTA_MU_STD'].median()
        # df.reset_index(inplace=True)

        worst, spots = [], 0
        for i, patient in enumerate(ponaqua_qualified):
            pdf = df.loc[df.PATIENT_ID == int(patient)]
            # print(i + 1, patient, len(pdf), 'spots')
            patient_mean_dist = pdf['DIST(mm)_MEAN'].mean()
            if not pdf.loc[pdf.GANTRY_ANGLE == 270.].empty:
                print(patient, pdf.loc[pdf.GANTRY_ANGLE == 270., 'BEAM_ID'].iloc[0], pdf.loc[pdf.GANTRY_ANGLE == 270., 'GANTRY_ANGLE'].iloc[0])

            worst.append(np.round(pdf['DIST(mm)_MEAN'].max(), 1))
            spots += len(pdf)
            print('\tGenauigkeit X: ', pdf['DELTA_X(mm)_MEAN'].mean(), '+-', pdf['DELTA_X(mm)_MEAN'].std(), 'mm', 'MAX', pdf.loc[abs(pdf['DELTA_X(mm)_MEAN']).idxmax(), 'DELTA_X(mm)_MEAN'])
            print('\tGenauigkeit Y: ', pdf['DELTA_Y(mm)_MEAN'].mean(), '+-', pdf['DELTA_Y(mm)_MEAN'].std(), 'mm', 'MAX', pdf.loc[abs(pdf['DELTA_Y(mm)_MEAN']).idxmax(), 'DELTA_Y(mm)_MEAN'])
            print('\tGenauigkeit MU:', pdf['DELTA_MU_MEAN'].mean(), '+-', pdf['DELTA_MU_MEAN'].std(), 'mm', 'MAX', pdf.loc[abs(pdf['DELTA_MU_MEAN']).idxmax(), 'DELTA_MU_MEAN'])
            print('\tReproduzier X: ', pdf['DELTA_X(mm)_STD'].mean(), '+-', pdf['DELTA_X(mm)_STD'].std(), 'mm', 'MAX', pdf.loc[abs(pdf['DELTA_X(mm)_STD']).idxmax(), 'DELTA_X(mm)_STD'])
            print('\tReproduzier Y: ', pdf['DELTA_Y(mm)_STD'].mean(), '+-', pdf['DELTA_Y(mm)_STD'].std(), 'mm', 'MAX', pdf.loc[abs(pdf['DELTA_Y(mm)_STD']).idxmax(), 'DELTA_Y(mm)_STD'])
            print('\tReproduzier MU:', pdf['DELTA_MU_STD'].mean(), '+-', pdf['DELTA_MU_STD'].std(), 'mm', 'MAX', pdf.loc[abs(pdf['DELTA_MU_STD']).idxmax(), 'DELTA_MU_STD'])
        
        print('ALL:', spots, 'spots')
        print('Genauigkeit: X', df['DELTA_X(mm)_MEAN'].mean(), '+-', df['DELTA_X(mm)_MEAN'].std(), 'mm', 'MAX', df.loc[abs(df['DELTA_X(mm)_MEAN']).idxmax(), 'DELTA_X(mm)_MEAN'])
        print('             Y', df['DELTA_Y(mm)_MEAN'].mean(), '+-', df['DELTA_Y(mm)_MEAN'].std(), 'mm', 'MAX', df.loc[abs(df['DELTA_Y(mm)_MEAN']).idxmax(), 'DELTA_Y(mm)_MEAN'])
        print('             MU', df['DELTA_MU_MEAN'].mean(), '+-', df['DELTA_MU_MEAN'].std(), 'mm', 'MAX', df.loc[abs(df['DELTA_MU_MEAN']).idxmax(), 'DELTA_MU_MEAN'])
        print('Reproduzier: X', df['DELTA_X(mm)_STD'].mean(), '+-', df['DELTA_X(mm)_STD'].std(), 'mm', 'MAX', df.loc[abs(df['DELTA_X(mm)_STD']).idxmax(), 'DELTA_X(mm)_STD'])
        print('             Y', df['DELTA_Y(mm)_STD'].mean(), '+-', df['DELTA_Y(mm)_STD'].std(), 'mm', 'MAX', df.loc[abs(df['DELTA_Y(mm)_STD']).idxmax(), 'DELTA_Y(mm)_STD'])
        print('             MU', df['DELTA_MU_STD'].mean(), '+-', df['DELTA_MU_STD'].std(), 'mm', 'MAX', df.loc[abs(df['DELTA_MU_STD']).idxmax(), 'DELTA_MU_STD'])

        sns.set(style='ticks', context='paper', font_scale=2.)
        sns.set_style({"xtick.direction": "in","ytick.direction": "in", "ax.axesbelow":True})

        # # scatterplots for Kristina - delete
        # # sns.scatterplot(data=df, x='DRILL_TIME(ms)', y='BEAM_CURRENT(A)')
        # sns.histplot(df['DRILL_TIME(ms)'], label=f"Mean: {df['DRILL_TIME(ms)'].mean().round(3)} ms")
        # # plt.yscale('log')
        # plt.legend()
        # plt.show()
        # return None

        # plt.rcParams.update({'font.size': 10})
        mm = 1 / 25.4
        figh, axh = plt.subplots(3, 2, figsize=(16, 12), constrained_layout=True)
        axh = axh.flatten()
        cmap = sns.color_palette('Paired', 8)
        subfigures = ['(a)', '(d)', '(b)', '(e)', '(c)', '(f)']

        bw1 = 5 / 80
        bw2 = 0.6 / 80
        ec = 'black'
        violin_lw = 1.0

        # original
        h1 = sns.histplot(df['DELTA_X(mm)_MEAN'],color=cmap[0], ax=axh[0], label='$\mu_{x}$', binrange=(-2.5, 2.5), binwidth=bw1, edgecolor=ec, zorder=2)
        h2 = sns.histplot(df['DELTA_Y(mm)_MEAN'],color=cmap[1], ax=axh[0], label='$\mu_{y}$', binrange=(-2.5, 2.5), binwidth=bw1, edgecolor=ec, zorder=1)
        h3 = sns.histplot(df['DELTA_X(mm)_STD'],color=cmap[6], ax=axh[1], label='$\sigma_{x}$', binrange=(0., 0.6), binwidth=bw2, edgecolor=ec, zorder=2)
        h4 = sns.histplot(df['DELTA_Y(mm)_STD'],color=cmap[7], ax=axh[1], label='$\sigma_{y}$', binrange=(0., 0.6), binwidth=bw2, edgecolor=ec, zorder=1)

        # COV evaluation
        # h1 = sns.histplot(df['DELTA_X(mm)_STD'] / df['DELTA_X(mm)_MEAN'],color=cmap[0], ax=axh[0], label='$COV_{x}$', binrange=(-2.5, 2.5), binwidth=bw1, edgecolor=ec, zorder=2)
        # h2 = sns.histplot(df['DELTA_Y(mm)_STD'] / df['DELTA_Y(mm)_MEAN'],color=cmap[1], ax=axh[0], label='$COV_{y}$', binrange=(-2.5, 2.5), binwidth=bw1, edgecolor=ec, zorder=1)

        axh[0].set_xlabel('$\mu_{x,y}$ [mm]')
        # axh[0].set_xlabel('$COV_{x,y}$')
        axh[0].set_xlim(-2.5, 2.5)
        axh[0].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        axh[0].grid(axis='both')
        axh[1].set_xlabel('$\sigma_{x,y}$ [mm]')
        axh[1].set_ylabel(None)
        axh[1].set_xlim(0.01, 0.6)
        axh[1].grid(axis='both')
        for i, ax in enumerate(axh):
            ax.legend(loc='upper right')
            # rectangle = Rectangle((0.1, 0.9), 0.9, 0.9, edgecolor='black', facecolor='white', fill=True)
            # ax.add_patch(rectangle)
            # rx, ry = rectangle.get_xy()
            # cx = rx + rectangle.get_width()/2.0
            # cy = ry + rectangle.get_height()/2.0
            # ax.annotate(subfigures[i], (0, 1), color='black', weight='bold', xycoords='axes fraction', ha='left', va='top')
            ax.text(-0.15, 1.1, f'{subfigures[i]}', color='black', fontweight='bold', fontsize=24., transform=ax.transAxes, ha='left', va='top')
            for container in ax.patches:
                # container.set_linewidth(violin_lw)
                pass

        # figh.tight_layout()
        # figh.savefig(os.path.join(out, 'histograms_all.png'), dpi=600)
        # plt.show()
        
        # sns.set(style='darkgrid', context='paper', font_scale=1.5)
        # fig, axs = plt.subplots(2, 1, figsize=(7, 7))
        target_cols = ['DELTA_X(mm)_MEAN', 'DELTA_X(mm)_STD', 'DELTA_Y(mm)_MEAN', 'DELTA_Y(mm)_STD']  #, ['DELTA_MU_MEAN', 'DELTA_MU_STD', 'MEDIAN_MU_STD']] #, ['DELTA_MU_MEAN', 'DELTA_MU_STD', 'MEDIAN_MU_STD']]
        ms = 10
        legend = [[Patch(facecolor=cmap[0], edgecolor=ec, alpha=0.8, label=f'$x$-accuracy'), 
                  Line2D([0], [0], marker='|', linewidth=8, color='black', markeredgecolor='white', markersize=ms, label='IQR & median')],
                  [Patch(facecolor=cmap[6], edgecolor=ec, alpha=0.8, label=f'$x$-reproducibility'), 
                  Line2D([0], [0], marker='|', linewidth=8, color='black', markeredgecolor='white', markersize=ms, label='IQR & median')],
                  [Patch(facecolor=cmap[1], edgecolor=ec, alpha=0.8, label=f'$y$-accuracy'), 
                  Line2D([0], [0], marker='|', linewidth=8, color='black', markeredgecolor='white', markersize=ms, label='IQR & median')],
                  [Patch(facecolor=cmap[7], edgecolor=ec, alpha=0.8, label=f'$y$-reproducibility'), 
                  Line2D([0], [0], marker='|', linewidth=8, color='black', markeredgecolor='white', markersize=ms, label='IQR & median')]]

        for i, (ax, col) in enumerate(zip(axh[2:], target_cols)):
            if i < 2:
                if i % 2 == 0:
                    vp1 = sns.violinplot(data=df, x='GANTRY_ANGLE', y=col, color=cmap[0], ax=ax, width=10, density_norm='width', cut=0, order=[float(angle) for angle in np.arange(360)], linewidth=violin_lw)
                else:
                    vp1 = sns.violinplot(data=df, x='GANTRY_ANGLE', y=col, color=cmap[6], ax=ax, width=10, density_norm='width', cut=0, order=[float(angle) for angle in np.arange(360)], linewidth=violin_lw)

            else:
                if i % 2 == 0:
                    vp1 = sns.violinplot(data=df, x='GANTRY_ANGLE', y=col, color=cmap[1], ax=ax, width=10, density_norm='width', cut=0, order=[float(angle) for angle in np.arange(360)], linewidth=violin_lw)
                else:
                    vp1 = sns.violinplot(data=df, x='GANTRY_ANGLE', y=col, color=cmap[7], ax=ax, width=10, density_norm='width', cut=0, order=[float(angle) for angle in np.arange(360)], linewidth=violin_lw)
 
            # bp1 = sns.boxplot(data=df, x='GANTRY_ANGLE', y=cols[1], color=cmap[i], ax=ax, width=10, order=np.arange(360), showfliers=True)
            for violin in ax.collections:
                violin.set_edgecolor('black')
            for box in ax.artists:
                box.set_edgecolor('black')

            for whisker in ax.lines:
                whisker.set_color('black')

            ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            if i % 2 == 0:
                ax.axhline(0., ls='--', c='black', zorder=0.8)
                ax.set_ylim(-2., 3.)
            else:
                ax.set_ylim(0, 1.)
            
            ax.set_xlabel(None)
            ax.set_xlim(-10, )
            ax.set_xticks([0, 45, 90, 135, 180, 225, 270, 315])
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
            ax.legend(handles=legend[i], loc='upper right')
            ax.set_axisbelow(True)
            ax.grid(axis='y')
            tick_every = 45
            # [tick.set_visible(False) for (i, tick) in enumerate(ax.get_xticklabels()) if i % tick_every != 0]
        
        # axh[0].legend(handles=legend, loc='upper right')
        # axs[0].set_title('Treatment fields grouped by gantry position')
        # axs[0].set_ylim(0, .8)
        axh[2].set_xlabel('Gantry angle [°]')
        axh[2].set_ylabel('$\mu_{x}$ [mm]')
        axh[3].set_xlabel('Gantry angle [°]')
        axh[3].set_ylabel('$\sigma_{x}$ [mm]')
        axh[4].set_ylabel('$\mu_{y}$ [mm]')
        axh[4].set_xlabel('Gantry angle [°]')
        axh[5].set_ylabel('$\sigma_{y}$ [mm]')
        axh[5].set_xlabel('Gantry angle [°]')
        # axs[2].set_ylim(0, 0.005)
        # fig.suptitle(f'Log-file Reproducibility', fontweight='bold')
        figh.tight_layout()
        figh.subplots_adjust(hspace=0.3)
        figh.subplots_adjust(wspace=0.2)
        figh.savefig(os.path.join(out, 'Fig2_violinplots_print.pdf'))
        
        # outliers
        Q1 = df.quantile(0.25)
        Q3 = df.quantile(0.75)
        IQR = Q3 - Q1

        # plt.show()


    def split_sigma(self):
        out = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\Logfiles\Reports'

        # locate existing dataframe
        single_stats_exists = False
        for file in sorted(os.listdir(self.df_destination)):
            if file.__contains__(str(self.patient_id)) and file.__contains__('sss_data') and file.endswith('.csv'):
                print(f'''Found patient spot statistics data '{file}', reading in..''')
                self.patient_sss_df = pd.read_csv(os.path.join(self.df_destination, file), dtype={'BEAM_ID':str})
                single_stats_exists = True

        if not single_stats_exists:
            print(f'/x\ No single spot statistics data found for patient-ID {self.patient_id}, exiting..')
            return None
        
        sns.set(style='darkgrid', context='talk', font_scale=1.2)
        ponaqua_qualified = [id.strip('\n') for id in open(r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\qualified_IDs.txt', 'r').readlines()]
        bs_file = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\cct_entities.csv'
        bs_data = pd.read_csv(bs_file, delimiter=';')
        bs_data = bs_data.loc[bs_data.PATIENT_ID.astype(str).isin(ponaqua_qualified)]
        bs_data.loc[bs_data.BODY_SITE == 'HNC/Paranasal/Eye'] = 'Oral cavity'
        bs_data.loc[bs_data.BODY_SITE == 'Skull/-base'] = 'Skull base'
        bs_data.loc[bs_data.BODY_SITE == 'Liver/Kidney'] = 'Kidney'
        bs_data.loc[bs_data.BODY_SITE == 'Lung/BCA'] = 'Lung'
        bs_data.loc[bs_data.BODY_SITE == 'Head/Glio'] = 'Brain'
        values = bs_data.value_counts('BODY_SITE').values
        def make_autopct(values):
            def my_autopct(pct):
                total = sum(values)
                val = int(round(pct*total/100.0))
                return '{v:d}'.format(v=val)
            return my_autopct
        
        bs_data.value_counts('BODY_SITE').plot.pie(y='BODY_SITE', autopct=make_autopct(values))
        plt.tight_layout()
        plt.savefig(os.path.join(out, f'ponaqua_piechart.png'), dpi=300)
        # plt.show()
    
        bs_dict = bs_data.to_dict()
        
        sns.set()
        for beam_id in self.patient_sss_df.BEAM_ID.drop_duplicates():
            beam_df = self.patient_sss_df.loc[self.patient_sss_df.BEAM_ID == beam_id]

            dim = int(np.ceil(np.sqrt(len(beam_df.LAYER_ID.drop_duplicates()))))
            if dim == 0:
                print(beam_id)
                continue
            fig, axs = plt.subplots(dim, dim, figsize=(10, 10), sharex=True, sharey=True, dpi=100)
            ax0 = fig.add_subplot(111, frameon=False)
            ax0.set_xticks([]), ax0.set_yticks([])
            fig.subplots_adjust(hspace=0.0, wspace=0.0)
            axs = axs.flatten()

            for i, layer_id in enumerate(beam_df.LAYER_ID.drop_duplicates()):
                layer_df = beam_df.loc[beam_df.LAYER_ID == layer_id]
                sns.scatterplot(data=layer_df, x='X_POSITION(mm)', y='Y_POSITION(mm)', hue='DELTA_X(mm)_STD', size='MU', ax=axs[i], legend=None, palette='YlOrBr')
                axs[i].annotate(f"$E =$ {np.round(layer_df['LAYER_ENERGY(MeV)'].iloc[0], 1)}MeV", xy=(5, 5), xycoords='axes pixels')
                # axs[i].grid()
                axs[i].set_xlabel(None), axs[i].set_ylabel(None)
                # if i != 0:
                #     axs[i].get_legend().remove()
            
            # sns.move_legend(axs[0], "upper right", bbox_to_anchor=(0, 1))
            
            ax0.set_xlabel('X position (mm)', labelpad=30)
            ax0.set_ylabel('Y position (mm)', labelpad=30)
            # fig.suptitle(f'''PATIENT {self.patient_id} ({bs_dict['BODY_SITE'][int(self.patient_id)]}) - BEAM {beam_id} ({beam_df.GANTRY_ANGLE.iloc[0]})''', fontweight='bold')
            plt.tight_layout()
            plt.savefig(os.path.join(out, f'{self.patient_id}_beam_{beam_id}_sigma_spotmap.png'), dpi=300)
            # plt.show()             
                
            continue
            
            fig, axs = plt.subplots(2, 1, figsize=(8, 10), sharex=True, sharey=True)
            # cols = ['DELTA_X(mm)_MEAN', 'DELTA_X(mm)_STD']
            # cols = ['MU', 'LAYER_ENERGY(MeV)']
            cols = ['LAYER_ENERGY(MeV)']
            h1 = sns.histplot(below[cols], ax=axs[0], kde=True)
            h2 = sns.histplot(above[cols], ax=axs[1], kde=True)
            # for ax in axs:
            #     ax.axvline(split_at, color='tab:orange')

            fig.suptitle(f'Beam {beam_df.BEAM_ID.iloc[0]} - {beam_df.GANTRY_ANGLE.iloc[0]}°', fontweight='bold')
            axs[0].set_title(f'$\sigma \leq$ 0.18')
            axs[1].set_title(f'$\sigma >$ 0.18')
            plt.tight_layout()
            plt.show()

        if int(self.patient_id) == 1700535:
            df1 = self.patient_sss_df.loc[(np.round(self.patient_sss_df['LAYER_ENERGY(MeV)'], 1) == 165.1) & (self.patient_sss_df.BEAM_ID == '1')]
            hist1 = self.patient_sss_df.loc[(self.patient_sss_df.BEAM_ID == '1')]
            df2 = self.patient_sss_df.loc[(np.round(self.patient_sss_df['LAYER_ENERGY(MeV)'], 1) == 165.9) & (self.patient_sss_df.BEAM_ID == '2')]
            hist2 = self.patient_sss_df.loc[(self.patient_sss_df.BEAM_ID == '2')]

            sns.set(style='ticks', context='paper', font_scale=1.5)
            sns.set_style({"xtick.direction": "in","ytick.direction": "in", "ax.axesbelow":True, "xtick.top":True, "ytick.right":True})
            fig1, ax = plt.subplots(1, 1, figsize=(7, 4), constrained_layout=True)
            fig2, axs = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(7, 4), constrained_layout=True)
            # cax = fig2.add_axes([0.1, 0., 0.8, 0.05])
            fig2.subplots_adjust(wspace=0.1)
            norm = plt.Normalize(0, 0.3)
            cmap = plt.get_cmap('YlOrBr')
            ec = 'black'
            h1 = sns.histplot(data=hist1, x='DELTA_X(mm)_STD', ax=ax, color=cmap(0.4), ec=ec, binwidth=0.3/80, binrange=(0, 0.3), label='Beam 1 - 180°', zorder=1)
            h2 = sns.histplot(data=hist2, x='DELTA_X(mm)_STD', ax=ax, color=cmap(0.9), ec=ec, binwidth=0.3/80, binrange=(0, 0.3), label='Beam 2 - 90°', zorder=2)
            s1 = sns.scatterplot(data=df1, x='X_POSITION(mm)', y='Y_POSITION(mm)', hue='DELTA_X(mm)_STD', hue_norm=norm, ec=ec, size='MU', size_norm=(0.02, 0.8), sizes=(30, 140), ax=axs[0], palette='YlOrBr', legend=None)
            s2 = sns.scatterplot(data=df2, x='X_POSITION(mm)', y='Y_POSITION(mm)', hue='DELTA_X(mm)_STD', hue_norm=norm, ec=ec, size='MU', size_norm=(0.02, 0.8), sizes=(30, 140), ax=axs[1], palette='YlOrBr', legend='brief')
            h,l = s2.get_legend_handles_labels()
            ax.grid(axis='both')
            ax.set_axisbelow(True)  # grid behind data
            ax.legend()
            ax.set_ylim(0.001, 350)
            # ax.text(0.02, 0.042, 'A', color='black', fontweight='bold', fontsize=24., transform=ax.transAxes, ha='left', va='bottom', bbox=dict(facecolor='white', edgecolor='black', pad=10.0))
            ax.text(-0.15, .92, '(a)', color='black', fontweight='bold', fontsize=24., transform=ax.transAxes, ha='left', va='bottom')
            axs[1].get_legend().remove()
            axs[1].legend(h[6:], l[6:], bbox_to_anchor=(0., 1.), loc=2, markerscale=1.)
            # axs[1].text(-0.15, 1.1, 'C', color='black', fontweight='bold', fontsize=24., transform=axs[1].transAxes, ha='left', va='bottom')
            axs[1].grid(axis='both')
            # axs[1].set_title(f'Beam 2 ({df2.GANTRY_ANGLE.iloc[0]}°)')
            ax.set_xlabel('$\sigma_{x}$ [mm]')
            ax.set_xlim(0, 0.3)
            axs[0].set_xlabel('$x$ [mm]'), axs[1].set_xlabel('$x$ [mm]')
            # axs[0].set_title(f'Beam 1 ({df1.GANTRY_ANGLE.iloc[0]}°)')
            axs[0].set_ylabel('$y$ [mm]')
            axs[0].text(-0.3, .9, '(b)', color='black', fontweight='bold', fontsize=24., transform=axs[0].transAxes, ha='left', va='bottom')
            axs[0].grid(axis='both')
            axs[0].set_ylim(-50, 90)
            axs[0].set_xlim(-100, 80)
            axs[0].set_xticks([-90, -60, -30,  0, 30, 60])
            axs[0].set_yticks([-60, -30,  0, 30, 60, 90])
            # for ax in axs:
            #     ax.grid()
            sm = plt.cm.ScalarMappable(cmap='YlOrBr', norm=norm)
            # cbar_ax = fig2.add_axes((.75, .11, .2, .77))
            # cbar_ax.set_visible(False)
            # cbar_ax.figure.colorbar(sm, ax=cbar_ax, label='$\sigma_{\Delta x}$ [mm]')
            # fig2.colorbar(sm, orientation='horizontal', cax=cax, label='$\sigma_{\Delta x}$ [mm]')
            fig2.colorbar(sm, ax=axs.ravel().tolist(), orientation='horizontal', aspect=40, label='$\sigma_{x}$ [mm]')
            # fig2.suptitle(f'''PATIENT {self.patient_id} ({bs_dict['BODY_SITE'][int(self.patient_id)]})''', fontweight='bold')
            fig1.savefig(os.path.join(out, f'{self.patient_id}_beam_{beam_id}_split_paper.pdf'))
            fig2.savefig(os.path.join(out, f'{self.patient_id}_beam_{beam_id}_sigma_paper.pdf')) 


    def prepare_psqa(self) -> None:
        treatment_records = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\Logfiles\extracted'
        
        # get original treatment deltaframe
        for csv in os.listdir(treatment_records):
            if csv.__contains__(str(self.patient_id)) and csv.__contains__('delta') and csv.endswith('.csv'):
                delta = pd.read_csv(os.path.join(treatment_records, csv), index_col='UNIQUE_INDEX', dtype={'FRACTION_ID':str, 'BEAM_ID':str})
            if csv.__contains__(str(self.patient_id)) and csv.__contains__('records') and csv.endswith('.csv'):
                records = pd.read_csv(os.path.join(treatment_records, csv), index_col='TIME', parse_dates=['TIME'], dtype={'FRACTION_ID':str, 'BEAM_ID':str})

        df = pd.concat([records.reset_index(), delta.reset_index()], axis=1)
        df = df.iloc[:,~df.columns.duplicated()]
        df = df[['FRACTION_ID', 'BEAM_ID', 'GANTRY_ANGLE', 'LAYER_ID', 'SPOT_ID', 'X_POSITION(mm)', 'Y_POSITION(mm)', 'MU', 'DELTA_X(mm)', 'DELTA_Y(mm)', 'DELTA_MU']]
        df['PLAN_X'] = df['X_POSITION(mm)'] - df['DELTA_X(mm)']   
        df['PLAN_Y'] = df['Y_POSITION(mm)'] - df['DELTA_Y(mm)'] 
        df.drop(columns=['DELTA_X(mm)', 'DELTA_Y(mm)'], inplace=True)

        treated_beams = df.BEAM_ID.drop_duplicates()
        qa_beams = self.patient_record_df.BEAM_ID.drop_duplicates()
        beam_dict = {}
        for tr_beam in treated_beams:
            for qa_beam in qa_beams:
                if (qa_beam.__contains__('WP') or qa_beam.__contains__('QA')) and qa_beam.__contains__(tr_beam):  # get any of the QA beams irradiated to WP
                    beam_dict[tr_beam] = qa_beam

        psqa_df = pd.DataFrame()
        for fx_no, fx_id in enumerate(df.FRACTION_ID.drop_duplicates()):
            fx_df = df.loc[df.FRACTION_ID == fx_id]
            for beam_id in fx_df.BEAM_ID.drop_duplicates():
                print(f'Start fraction {fx_id} ({fx_no + 1}/{len(df.FRACTION_ID.drop_duplicates())}) .. Beam {beam_id}', end='\r')
                try:
                    last_qa_fx = self.patient_record_df.loc[self.patient_record_df.BEAM_ID == beam_dict[beam_id]].FRACTION_ID.drop_duplicates()[-1]
                    beam_df_treat = fx_df.loc[fx_df.BEAM_ID == beam_id]
                    beam_df_qa = self.patient_record_df.loc[(self.patient_record_df.BEAM_ID == beam_dict[beam_id]) & (self.patient_record_df.FRACTION_ID == last_qa_fx)]

                except KeyError:
                    print(f'  /!\ Beam {beam_id} not in records, will be skipped')
                    continue
                
                # OLD STYLE
                # iterations = 0
                # while len(beam_df_qa) != len(beam_df_treat):
                #     for i, (x_qa, x_tr) in enumerate(zip(beam_df_qa['X_POSITION(mm)'], beam_df_treat['X_POSITION(mm)'])):
                #         if abs(x_qa - x_tr) > 3.:
                #             if len(beam_df_qa) > len(beam_df_treat):
                #                 beam_df_qa.drop(beam_df_qa.index[i], inplace=True)
                #             else:
                #                 beam_df_treat.drop(beam_df_treat.index[i], inplace=True)
                            
                #             break
                    
                #     iterations += 1
                #     if iterations > 10:
                #         break
                # if len(beam_df_qa) == len(beam_df_treat):
                #         temp = beam_df_treat.copy()
                #         temp['X_PSQA'] = np.array(beam_df_qa['X_POSITION(mm)'])
                #         temp['Y_PSQA'] = np.array(beam_df_qa['Y_POSITION(mm)'])
                #         temp['MU_PSQA'] = np.array(beam_df_qa['MU'])
                #         psqa_df = pd.concat([psqa_df, temp])
                #     else:
                #         continue

                # NEW STYLE
                for layer_id in beam_df_treat.LAYER_ID.drop_duplicates():
                    layer_df_treat = beam_df_treat.loc[beam_df_treat.LAYER_ID == layer_id]
                    layer_df_qa = self.patient_record_df.loc[(self.patient_record_df.BEAM_ID == beam_dict[beam_id]) & (self.patient_record_df.FRACTION_ID == last_qa_fx) & (self.patient_record_df.LAYER_ID == layer_id)]

                    need_sorting = False
                    if len(layer_df_qa) == len(layer_df_treat):
                        max_x = max(abs(np.array(layer_df_qa['X_POSITION(mm)']) - np.array(layer_df_treat['X_POSITION(mm)'])))
                        max_y = max(abs(np.array(layer_df_qa['Y_POSITION(mm)']) - np.array(layer_df_treat['Y_POSITION(mm)'])))
                        if max_x > 3. or max_y > 3.:
                            need_sorting = True
                    else:
                        need_sorting = True

                    layer_copy = layer_df_treat.copy()
                    x_qa_sorted, y_qa_sorted = [], []
                    if need_sorting:
                        for i, (x_treat, y_treat) in enumerate(zip(layer_df_treat['X_POSITION(mm)'], layer_df_treat['Y_POSITION(mm)'])):
                            found = False
                            for j, (x_qa, y_qa) in enumerate(zip(layer_df_qa['X_POSITION(mm)'], layer_df_qa['Y_POSITION(mm)'])):
                                dist = (x_treat - x_qa) ** 2 + (y_treat - y_qa) ** 2
                                if dist < 9:  # < 3mm euclidean distance
                                    x_qa_sorted.append(x_qa)
                                    y_qa_sorted.append(y_qa)
                                    found = True
                                    break
                            
                            if not found:
                                x_qa_sorted.append(np.nan)
                                y_qa_sorted.append(np.nan)

                        # plt.plot(layer_df_treat['X_POSITION(mm)'],layer_df_treat['Y_POSITION(mm)'], 'o-', label='Treat')
                        # plt.plot(x_qa_sorted, y_qa_sorted, 'o-', label='QA')
                        # plt.legend()
                        # plt.show()
                        layer_copy['X_PSQA'] = np.array(x_qa_sorted)
                        layer_copy['Y_PSQA'] = np.array(y_qa_sorted)
                    
                    else:
                        layer_copy['X_PSQA'] = np.array(layer_df_qa['X_POSITION(mm)'])
                        layer_copy['Y_PSQA'] = np.array(layer_df_qa['Y_POSITION(mm)'])
                        
                    psqa_df = pd.concat([psqa_df, layer_copy])
                
            print(f'Start fraction {fx_id} ({fx_no + 1}/{len(df.FRACTION_ID.drop_duplicates())}) .. DONE            ', end='\n')
                    

        # psqa_df.dropna(axis=0, inplace=True)
        psqa_df['DELTA_TO_QA_X'] = psqa_df['X_POSITION(mm)'] - psqa_df['X_PSQA']
        psqa_df['DELTA_TO_QA_Y'] = psqa_df['Y_POSITION(mm)'] - psqa_df['Y_PSQA']

        print(f'SUCCESS - Writing data to "patient_{self.patient_id}_diff_to_PSQA.csv"')
        psqa_df.to_csv(os.path.join(self.df_destination, f'patient_{self.patient_id}_diff_to_PSQA.csv'))
        
        
        # indices = psqa_df.loc[abs(psqa_df.DELTA_TO_QA_X) > 5].index
        # already = []
        # for i in indices:
        #     fx, beam, layer = psqa_df.FRACTION_ID.loc[i], psqa_df.BEAM_ID.loc[i], psqa_df.LAYER_ID.loc[i]
        #     if (beam, layer) in already: continue
        #     layer_df = psqa_df.loc[(psqa_df.FRACTION_ID == fx) & (psqa_df.BEAM_ID == beam) & (psqa_df.LAYER_ID == layer)]
        #     plt.plot(layer_df['X_POSITION(mm)'],layer_df['Y_POSITION(mm)'], 'o-', label='treat')
        #     plt.plot(layer_df['X_PSQA'],layer_df['Y_PSQA'], 'o-', label='QA')
        #     plt.legend()
        #     print(layer_df.DELTA_TO_QA_X.max())
        #     plt.show()
        #     already.append((beam, layer))

        # sns.set(style='whitegrid')
        # sns.histplot(psqa_df['X_POSITION(mm)'] - psqa_df['X_PSQA'], label='x-position')
        # sns.histplot(psqa_df['Y_POSITION(mm)'] - psqa_df['Y_PSQA'], label='y-position')
        # plt.legend()
        # plt.xlabel('Treat - QA log file diff [mm]')
        # plt.ylabel('Count')
        # plt.title(f'{self.patient_id} - Treatment vs. Water Phantom')
        # plt.xlim(-3., 3.)
        # # plt.yscale('log')
        # plt.savefig(rf'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\Logfiles\Reports\{self.patient_id}_diff_to_QA.png')
        # plt.clf()


    def eval_psqa(self):
        psqa_records_dir = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\Logfiles\PSQA\extracted'
        other_dfs, loaded = [], 0
        for i, csv in enumerate(os.listdir(psqa_records_dir)):
            if csv.__contains__('diff') and csv.endswith('.csv'):
                id = csv.split('_')[1]
                df = pd.read_csv(os.path.join(psqa_records_dir, csv), dtype={'BEAM_ID':str, 'FRACTION_ID':str})
                df['PATIENT_ID'] = id
                for beam in df.BEAM_ID.drop_duplicates():
                    if len(df.loc[df.BEAM_ID == beam, 'FRACTION_ID'].drop_duplicates()) < 15:
                        df.drop(df.loc[df.BEAM_ID == beam].index, inplace=True)
                        print(id, 'dropped beam', beam)
                other_dfs.append(df)
                loaded += 1
                print(f'Loaded {loaded}/14 datasets ..', end='\r')
        
        df = pd.concat(other_dfs).reset_index()
        df.dropna(inplace=True)

        print(f'\tMin\tMax\tMean\tStd\tMed')
        print(f'X[mm]:\t{df.DELTA_TO_QA_X.min().round(3)}\t{df.DELTA_TO_QA_X.max().round(3)}\t{df.DELTA_TO_QA_X.mean().round(3)}\t{df.DELTA_TO_QA_X.std().round(3)}\t{df.DELTA_TO_QA_X.median().round(3)}')
        print(f'Y[mm]:\t{df.DELTA_TO_QA_Y.min().round(3)}\t{df.DELTA_TO_QA_Y.max().round(3)}\t{df.DELTA_TO_QA_Y.mean().round(3)}\t{df.DELTA_TO_QA_Y.std().round(3)}\t{df.DELTA_TO_QA_Y.median().round(3)}')

        # indices = df.loc[abs(df.DELTA_TO_QA_Y) > 10].index
        # for i in indices:
        #     fx, beam, layer = df.FRACTION_ID.loc[i], df.BEAM_ID.loc[i], df.LAYER_ID.loc[i]
        #     layer_df = df.loc[(df.FRACTION_ID == fx) & (df.BEAM_ID == beam) & (df.LAYER_ID == layer)]
        #     plt.plot(layer_df['X_POSITION(mm)'],layer_df['Y_POSITION(mm)'], 'o-', label='treat')
        #     plt.plot(layer_df['X_PSQA'],layer_df['Y_PSQA'], 'o-', label='QA')
        #     plt.legend()
        #     print(layer_df.DELTA_TO_QA_X.max())
        #     print(layer_df.loc[i])
        #     plt.show()

        # all patients
        # sns.set(style='whitegrid')
        # sns.histplot(df.DELTA_TO_QA_X, label='x-position')#, binrange=(-3, 3))        
        # sns.histplot(df.DELTA_TO_QA_Y, label='y-position')#, binrange=(-3, 3))  
        # plt.legend()
        # plt.show()

        delta_treat, delta_qa = [], []
        for pat in df.PATIENT_ID.drop_duplicates():
            pat_df = df.loc[df.PATIENT_ID == pat]
            for beam in pat_df.BEAM_ID.drop_duplicates():
                slice = pat_df.loc[pat_df.BEAM_ID == beam]
                delta_treat.append((slice.GANTRY_ANGLE.iloc[0], (slice['X_POSITION(mm)'] - slice.PLAN_X).mean(), (slice['X_POSITION(mm)'] - slice.PLAN_X).std()))
                delta_qa.append((slice.GANTRY_ANGLE.iloc[0], (slice.X_PSQA - slice.PLAN_X).mean(), (slice.X_PSQA - slice.PLAN_X).std()))

        sns.set(style='ticks')
        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        for i, (treat_point, qa_point) in enumerate(zip(delta_treat, delta_qa)):
            if i == len(delta_qa) - 1:
                ax.errorbar(x=treat_point[0], y=treat_point[1], yerr=treat_point[2], color='tab:blue', marker='o', capsize=3, ls='None',  label='Treatment: X° mean $\pm$ 1$\sigma$', zorder=2)
                ax.errorbar(x=qa_point[0], y=qa_point[1], yerr=qa_point[2], color='red', marker='s', capsize=3, ls='None', label='QA: 0° mean $\pm$ 1$\sigma$', zorder=3)
            else:
                ax.errorbar(x=treat_point[0], y=treat_point[1], yerr=treat_point[2], color='tab:blue', marker='o', capsize=3, ls='None', zorder=2)
                ax.errorbar(x=qa_point[0], y=qa_point[1], yerr=qa_point[2], color='red', marker='s', capsize=3, ls='None', zorder=3)
        
        ax.axhline(0., color='lightgray', label='Planned')
        ax.set_xlabel('Gantry angle [°]')
        ax.set_ylabel('Spot $\Delta x$ to plan [mm]')
        ax.set_ylim(-1.5, 1.5)
        ax.grid(False)
        ax.legend()
        fig.tight_layout()
        fig.savefig(r'N:\fs4-HPRT\HPRT-Docs\Lukas\Publications\01_LogfileUncertainty\figures\supplement04.png', dpi=600)
        fig.show()


if __name__ == '__main__':
    root_dir = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\Logfiles\converted'
    # root_dir = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\01_SpotShape\Logfiles_Spotshape_QA\converted'
    ponaqua_qualified = [id.strip('\n') for id in open(r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\qualified_IDs.txt', 'r').readlines()]
    for id in ponaqua_qualified:
        log = MachineLog(os.path.join(root_dir, id))
        # log.prepare_dataframe()
        # log.prepare_deltaframe()
        # log.prepare_sss_dataframe()
        # log.plan_creator(fraction='all', mode='all')
        # log.beam_timings()
        # log.sss_boxplot()
        # log.split_sigma()
    
