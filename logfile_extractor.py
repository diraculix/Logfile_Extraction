'''
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
from scipy import optimize
from tkinter import Tk, filedialog
from matplotlib import pyplot as plt

output_dir = r'N:\fs4-HPRT\HPRT-Docs\Lukas\Logfile_Extraction\output'  # set output dir, fallback alternative below (if code is used between multiple systems)
if not os.path.isdir(os.path.join(output_dir, '..')):
    output_dir = r'/home/luke/Logfile_Extraction/output'

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
                    print(f'''Chosen path '{self.logfile_dir}' may only contain directories (one per fraction). Please retry..''')
                    break
                elif index == len(os.listdir(self.logfile_dir)) - 1:
                    valid_dir = True
        
        # set dataframe storage directory, fallback alternative below
        # self.df_destination = r'N:\fs4-HPRT\HPRT-Docs\Lukas\Logfile_Extraction\dataframes'
        self.df_destination = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\Logfiles\extracted'
        if not os.path.isdir(os.path.join(self.df_destination, '..')):            
            self.df_destination = r'/home/luke/Logfile_Extraction/dataframes'
        
        if not os.path.isdir(self.df_destination):
            os.mkdir(self.df_destination)

        # initialize fraction and beam list by directory names
        self.fraction_list = sorted(os.listdir(self.logfile_dir))
        self.num_fractions = len(self.fraction_list)
        self.beam_list = []
        for f in self.fraction_list:
            beams_in_frac = []
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
                self.patient_record_df = pd.read_csv(os.path.join(self.df_destination, fname), index_col='TIME', dtype={'BEAM_ID':str, 'FRACTION_ID':str})
                record_df_exists = True
            elif fname.__contains__(f'{self.patient_id}_tuning') and fname.endswith('.csv'):
                self.tuning_df_name = fname
                print(f'''Found patient tuning dataframe '{self.tuning_df_name}', reading in..''')
                self.patient_tuning_df = pd.read_csv(os.path.join(self.df_destination, fname), index_col='TIME', dtype={'BEAM_ID':str, 'FRACTION_ID':str})
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

                            except:  # log-file may not be read if interlock caused restart of logging in same file 
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
                                record_file_df = record_file_df.loc[:, :'Y_POSITION(mm)']                       # slice dataframe, drop redundant columns
                                record_file_df['DOSE_PRIM(C)'] = charge_col
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
                                
                                # drop unusable rows AFTER charge extraction
                                record_file_df.drop(record_file_df.loc[(record_file_df['SUBMAP_NUMBER'] == current_spot_submap) & ((record_file_df['X_POSITION(mm)'] == -10000.0) | (record_file_df['Y_POSITION(mm)'] == -10000.0))].index, inplace=True)  # drop unusable rows for current submap
                                
                                # average over all spot entries for most accurate position/shape (recommended by IBA)                                
                                mean_xpos, mean_ypos = record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['X_POSITION(mm)']].mean().iloc[0], record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['Y_POSITION(mm)']].mean().iloc[0]
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

                                previous_xpos, previous_ypos, previous_spot_submap = mean_xpos, mean_ypos, current_spot_submap
                                current_spot_submap = record_file_df.loc[record_file_df['SUBMAP_NUMBER'] > current_spot_submap]['SUBMAP_NUMBER'].min()  # proceed to next submap
                                if not split_spot: current_spot_id += 1  # keep track of spot id
                            
                            record_file_df.drop(columns=['DOSE_PRIM(C)'], inplace=True)
                            record_file_df.drop_duplicates(subset=['SUBMAP_NUMBER'], keep='last', inplace=True)  # keep only last entry for each spot
                            record_file_df.index = record_file_df['TIME']  # change to datetime index AFTER filtering, timestamp is NOT unique!
                            record_file_df.drop(columns=['TIME', 'SUBMAP_NUMBER'], inplace=True)
                            record_file_df.reindex()

                            # draw machine parameters from *map_specif*.csv
                            for specif_file in record_specifs:  
                                if int(specif_file.split('_')[2].split('_')[0]) == layer_id:
                                    with open(specif_file, 'r') as specif_file:
                                        lines = specif_file.readlines()
                                        ic_offsets = lines[3]
                                        ic_offset_x, ic_offset_y = float(ic_offsets.split(',')[3]), float(ic_offsets.split(',')[2])  # coordinate system switch x <--> y
                                        range_at_degrader = float(lines[1].split(',')[1])  # beam range after energy selection system
                            
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
        qa_spots = [18, 22]
        qa_xy = [0., 30., 60., 120.]
        qa_xy += [- xy for xy in qa_xy]

        print(f'\nMining semi-annual spot-QA data..')
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
                
                if not gantry_angle in qa_angles:
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
                            for submap in record_file_df['SUBMAP_NUMBER'].drop_duplicates():
                                submap_df = record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == submap]
                                record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == submap, ['X_POS_IC23(mm)']] = submap_df['X_POSITION(mm)'].mean()
                                record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == submap, ['Y_POS_IC23(mm)']] = submap_df['Y_POSITION(mm)'].mean()
                                record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == submap, ['X_WID_IC23(mm)']] = submap_df['X_WIDTH(mm)'].mean()
                                record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == submap, ['Y_WID_IC23(mm)']] = submap_df['Y_WIDTH(mm)'].mean()

                            record_file_df.drop_duplicates(subset=['SUBMAP_NUMBER'], keep='last', inplace=True)  # keep only last entries for each spot (most accurate)
                            record_file_df.index = record_file_df['TIME']   
                            record_file_df.drop(columns=['TIME'], inplace=True)

                            if not len(record_file_df[['X_POSITION(mm)', 'Y_POSITION(mm)']].drop_duplicates()) in qa_spots:
                                continue
                            
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
                            record_file_df.drop_duplicates(subset=['X_POSITION(mm)', 'Y_POSITION(mm)'], inplace=True)
                            # record_file_df['X_POS_IC23(mm)'], record_file_df['Y_POS_IC23(mm)'] = record_file_df['X_POSITION(mm)'], record_file_df['Y_POSITION(mm)']
                            # record_file_df['X_WID_IC23(mm)'], record_file_df['Y_WID_IC23(mm)'] = record_file_df['X_WIDTH(mm)'], record_file_df['Y_WIDTH(mm)']
                            record_file_df['X_POSITION(mm)'] = record_file_df[['Y_POS_IC23(mm)']].apply(map_spot_pos, args=(ic_offset_x, sad_x, ictoiso_x))
                            record_file_df['Y_POSITION(mm)'] = record_file_df[['X_POS_IC23(mm)']].apply(map_spot_pos, args=(ic_offset_y, sad_y, ictoiso_y))
                            record_file_df['X_WIDTH(mm)'] = record_file_df[['Y_WID_IC23(mm)']].apply(map_spot_width, args=(sad_x, ictoiso_x))
                            record_file_df['Y_WIDTH(mm)'] = record_file_df[['X_WID_IC23(mm)']].apply(map_spot_width, args=(sad_y, ictoiso_y))
                            record_file_df.drop(columns=['X_POS_IC23(mm)', 'Y_POS_IC23(mm)', 'X_WID_IC23(mm)', 'Y_WID_IC23(mm)'], inplace=True)
                            record_file_df['DIST_TO_ISO(mm)'] = np.sqrt(np.square(record_file_df['X_POSITION(mm)']) + np.square(record_file_df['Y_POSITION(mm)']))
                            record_file_df.reindex()  # make sure modified layer df is consistent with indexing
                            to_do_layers.append(record_file_df)
                    
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
                    layer_df = layer_df.loc[layer_df['X_ROUND'].isin(qa_xy) & layer_df['Y_ROUND'].isin(qa_xy)]
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
        self.record_df_name = 'QA_2017-2022_records_data.csv'
        
        print(f'''  ..Writing dataframe to '{self.df_destination}' as .CSV.. ''')
        while True:   
            try:
                self.qa_record_df.to_csv(self.record_df_name)
                break
            except PermissionError:
                input('  Permission denied, close target file and press ENTER.. ')

        self.qa_record_df = pd.read_csv('QA_2017-2022_records_data.csv', dtype={'BEAM_ID':str, 'FRACTION_ID':str})
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
            plan_dir = r'/home/luke/Logfile_Extraction/1676348/DeliveredPlans'        
        
        for file in os.listdir(plan_dir):
            if file.endswith('.dcm'):
                ds = pydicom.read_file(os.path.join(plan_dir, file))
                try:
                    beam_seq = ds.IonBeamSequence
                except AttributeError:
                    if verbose: print(f'  /x\ DICOM {file} is SOP storage class, not actual RTPLAN')
                    continue

                for beam in ds.IonBeamSequence:
                    plan_energies = np.array(pd.Series(sorted([layer.NominalBeamEnergy for layer in beam.IonControlPointSequence])).drop_duplicates().to_list())
                    
                    # check gantry angle, total layers, beam energies. Do not check names, they are non-standardized
                    if float(beam.IonControlPointSequence[0].GantryAngle) == gtr_angle and len(beam.IonControlPointSequence) == n_layers * 2:  # layer sequence in DICOM has double entries, BUT first/last layer might be missing due to single-spot layer tuning
                        try:
                            max_energy_diff = np.max(np.abs(plan_energies - log_energies))
                            if max_energy_diff < 0.1:  # tolerance on energy diff to plan, logs normally below dE = 0.1MeV
                                plan_dcm = os.path.join(plan_dir, file)
                                beam_ds = beam
                                found = True
                    
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
                                            else: return None
                                    elif plan_e == min(plan_energies):
                                        plan_dcm = os.path.join(plan_dir, file)
                                        beam_ds = beam
                                        found = True
                                        if verbose: 
                                            n_spots = beam.IonControlPointSequence[-2].NumberOfScanSpotPositions
                                            if n_spots <= 1: print(f'  /!\ Last layer with energy {plan_e} contains {n_spots} spot(s) - omitted by tuning?')
                                            else: return None
                                    else:
                                        if verbose: print(f'  /x\ Log-file is missing energy {plan_e}, skipping beam {beam_id}..')
                                        return None

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
        scope_record_df = scope_record_df.loc[scope_record_df['FRACTION_ID'] == scope_record_df['FRACTION_ID'].max()]  # analyze only last irradiated selected beam
        scope_tuning_df = self.patient_tuning_df.loc[self.patient_tuning_df['BEAM_ID'] == beam_id]
        scope_tuning_df = scope_tuning_df.loc[scope_tuning_df['FRACTION_ID'] == scope_record_df['FRACTION_ID'].max()]
        print('Selected beam is:', beam_id)

        print('\nGenerating layer plot..')
        plan_dcm, beam_ds = self.dicom_finder(fraction_id=scope_record_df['FRACTION_ID'].iloc[0], beam_id=beam_id, verbose=True)
        beam_sorting_dict = self.spot_sorter(fraction_id=scope_record_df['FRACTION_ID'].iloc[0], beam_id=beam_id)

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
            
            # if beam_id == '2' and layer_id == 9:
            #     fig2, ax2 = plt.subplots(figsize=(6, 6))
            #     ax2.plot(plan_x_positions, plan_y_positions, marker='x', linestyle='-', label='Planned')
            #     # ax2.plot(*zip(*spot_points_log), marker='o', markerfacecolor='None', linestyle='--', color='grey', label='Log-file original')
            #     ax2.plot(*zip(*spot_points_sorted), marker='o', markerfacecolor='None', linestyle='-', color='black', label='Log-file')
            #     ax2.plot(*zip(*tuning_points_log), marker='o', markerfacecolor='None', linestyle='None', color='limegreen', label='Tuning spot(s)')
            #     # ax2.annotate(f'Beam {beam_id} | Layer #{str(layer_id + 1).zfill(2)} | $\Delta$ = {abs(len(plan_x_positions) - len(x_positions))}', xy=(1.0, 1.0), xycoords='axes points')
            #     ax2.set_xlabel('Spot $x$-position @ISO [mm]')
            #     ax2.set_ylabel('Spot $y$-position @ISO [mm]')
            #     ax2.legend()
            #     fig2.tight_layout()
            #     fig2.savefig(f'{output_dir}/{self.patient_id}_{beam_id}_spotmap_show.png', dpi=600)
            #     plt.show()
               
            #     return None

            axs[layer_id].plot(plan_x_positions, plan_y_positions, marker='x', linestyle='-', markersize=2.0, markeredgewidth=0.2, linewidth=0.2, label='Planned')
            # MU colour coding still bugged, needs deltaframe for functionality
            axs[layer_id].scatter(*zip(*spot_points_sorted), c=spot_mu_delta, cmap='bwr', vmin=-0.1, vmax=0.1, edgecolors='black', linewidths=0.2, s=7, label='Log-file sorted')
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
                    plan_dcm, beam_ds = self.dicom_finder(fraction_id=fx_id, beam_id=beam_id, verbose=False)
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

            print(f'  ..Fraction {f + 1}/{self.num_fractions} complete..')

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


    def delta_correlation_matrix(self, gtr_only=True):
        other_records, other_deltas = [], []
        for file in sorted(os.listdir(self.df_destination)):
            if file.__contains__(str(self.patient_id)) and file.__contains__('delta') and file.endswith('.csv'):
                print(f'''Found patient deltaframe '{file}', reading in..''')
                self.patient_delta_df = pd.read_csv(os.path.join(self.df_destination, file), index_col='UNIQUE_INDEX', dtype={'BEAM_ID':str, 'FRACTION_ID':str})
            if file.__contains__('records') and file.endswith('.csv'):
                other_records.append(os.path.join(self.df_destination, file))           
            if file.__contains__('delta') and file.endswith('.csv'):
                other_deltas.append(os.path.join(self.df_destination, file))   

        try:
            delta_shape = self.patient_delta_df.shape
            if self.patient_record_df.shape != delta_shape:
                print(f'  /x\ Dataframe shape does not match deltaframe shape [{self.patient_record_df.shape} vs. {self.patient_delta_df.shape}], proceed with caution..')

        except AttributeError:
            print(f'''\nUnable to locate patient deltaframe for patient-ID {self.patient_id}, calling prepare_deltaframe()..''')
            self.prepare_deltaframe()
        
        to_drop, to_concat = ['FRACTION_ID', 'LAYER_ID', 'TOTAL_LAYERS', 'SPOT_ID'], []

        print('Gathering correlation data from patient database..')
        for this_record in other_records:
            has_delta = False
            this_patient_id = this_record.split('\\')[-1].split('_')[1]
            for this_delta in other_deltas:
                this_delta_id = this_delta.split('\\')[-1].split('_')[1]
                if this_patient_id == this_delta_id:
                    this_record_df = pd.read_csv(this_record, index_col='TIME', dtype={'BEAM_ID':str, 'FRACTION_ID':str})
                    this_joint_df = pd.read_csv(this_delta, index_col='UNIQUE_INDEX', dtype={'BEAM_ID':str, 'FRACTION_ID':str})
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
                    to_concat.append(this_joint_df)
                else:
                    print(f'  /!\ Dataframe shapes do not match [{this_record_df.shape} vs. {this_joint_df.shape}]')
                    # continue
            
            else:
                this_joint_df.drop(columns=to_drop, inplace=True)
                to_concat.append(this_joint_df)
        
        joint_df = pd.concat(to_concat, ignore_index=True)
        corr_matrix = joint_df.corr(method='pearson')  # method='spearman'?
        fig, ax = plt.subplots(1, 1, figsize=(16, 10))
        sns.heatmap(corr_matrix, annot=True, cmap='icefire', vmin=-1.0, vmax=1.0)
        plt.tight_layout()
        plt.savefig(f'{output_dir}/N={len(other_records)}_correlation_matrix.png', dpi=150)
        plt.clf()
        # fig, ax = plt.subplots(1, 1, figsize=(15, 15))
        # sns.pairplot(joint_df, vars=['DELTA_X(mm)', 'DELTA_Y(mm)', 'DELTA_MU', 'MU', 'LAYER_ENERGY(MeV)'], hue='GANTRY_ANGLE')
        # try:
        #     g = sns.pairplot(joint_df, vars=['DELTA_X(mm)', 'DELTA_Y(mm)', 'LAYER_ENERGY(MeV)', 'MU', 'GANTRY_ANGLE'], plot_kws=dict(s=10, edgecolor=None, alpha=0.3), corner=True)      
        # except:
        #     g = sns.pairplot(joint_df, vars=['DELTA_X(mm)', 'DELTA_Y(mm)', 'LAYER_ENERGY(MeV)', 'MU', 'GANTRY_ANGLE'], plot_kws=dict(s=10, edgecolor=None, alpha=0.3))      
        # # g._legend.remove()
        # plt.legend(bbox_to_anchor=(1, 1))
        # plt.tight_layout()
        # plt.savefig(f'{output_dir}/{self.patient_id}_pairplot.png', dpi=1000)


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
        beam_list = self.patient_record_df['BEAM_ID'].drop_duplicates()
        indices = beam_list.index.to_list()

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
                key = int(input('\n Select beam key: '))
                if key > len(beam_list) or key <= 0:
                    print('Key out of bounds, select another..')
                    continue
                else:
                    break
            except:
                print('Invalid input, try again..')

        beam_id = str(beam_list[key - 1])

        fig, axs = plt.subplots(4, 1, sharex=True, figsize=(10, 9), gridspec_kw={'height_ratios': [3, 1, 1, 1]})
        # fig.subplots_adjust(hspace=0.0, wspace=0.0)
        ax0 = fig.add_subplot(111, frameon=False)
        ax0.set_xticks([]), ax0.set_yticks([])
        ax0.set_ylabel('Time [s]', labelpad=30, fontweight='bold')
        ax0.set_xlabel('Fraction-ID', labelpad=30, fontweight='bold')
        # ax0.set_title(f'Beam timings for patient-ID {self.patient_id}', fontweight='bold')
        axs.flatten()
        x_axis = self.fraction_list
        # for ax in axs:
        #     ax.set_xticks(range(len(x_axis)))
        #     ax.set_xticklabels(range(1, len(x_axis) + 1))

        print('\nGenerating beam timing plot..')
        global_drills, global_spot_switches, global_energy_switches, totals = [], [], [], []
        for x, fx_id in enumerate(x_axis):  # remember type(fx_id) = <str>
            beam_df = self.patient_record_df.loc[(self.patient_record_df['BEAM_ID'] == beam_id) & (self.patient_record_df['FRACTION_ID'] == fx_id)]
            beam_tuning_df = self.patient_tuning_df.loc[(self.patient_tuning_df['BEAM_ID'] == beam_id) & (self.patient_tuning_df['FRACTION_ID'] == fx_id)]
            if beam_df.empty or beam_tuning_df.empty:
                print('empty')
                global_drills.append(0.0), global_spot_switches.append(0.0), global_energy_switches.append(0.0)
                continue

            total_drill_time = (beam_df['DRILL_TIME(ms)'].sum() + beam_tuning_df['DRILL_TIME(ms)'].sum()) / 1000
            total_layer_time, total_energy_switch = 0.0, 0.0
            layer_dfs = [beam_df.loc[beam_df['LAYER_ID'] == lid] for lid in beam_df['LAYER_ID'].drop_duplicates()]
            layer_tuning_dfs = [beam_tuning_df.loc[beam_tuning_df['LAYER_ID'] == lid] for lid in beam_tuning_df['LAYER_ID'].drop_duplicates()]
            for layer_id, layer_df in enumerate(layer_dfs):
                start_this = pd.to_datetime(layer_tuning_dfs[layer_id].first_valid_index())
                end_this = pd.to_datetime(layer_df.last_valid_index())
                layer_time = end_this - start_this
                total_layer_time += layer_time.total_seconds()
                
                if layer_id > 0:
                    end_previous = pd.to_datetime(layer_dfs[layer_id - 1].last_valid_index())
                    energy_switch = start_this - end_previous
                    total_energy_switch += energy_switch.total_seconds()
            
            total_spot_switch = total_layer_time - total_drill_time
            global_drills.append(total_drill_time), global_spot_switches.append(total_spot_switch), global_energy_switches.append(total_energy_switch)
            totals.append(total_drill_time + total_spot_switch + total_energy_switch)
       
        print(f'\n\t\tMedian\tMin\tMax\tSigma')
        print(f'Drill\t\t{np.median(global_drills):.3f}\t{min(global_drills):.3f}\t{max(global_drills):.3f}\t{np.std(global_drills):.3f}')
        print(f'Spot switch\t{np.median(global_spot_switches):.3f}\t{min(global_spot_switches):.3f}\t{max(global_spot_switches):.3f}\t{np.std(global_spot_switches):.3f}')
        print(f'Energy switch\t{np.median(global_energy_switches):.3f}\t{min(global_energy_switches):.3f}\t{max(global_energy_switches):.3f}\t{np.std(global_energy_switches):.3f}')
        print('_____________________________________________________')
        print(f'Total beamtime\t{np.median(totals):.3f}\t{min(totals):.3f}\t{max(totals):.3f}\t{np.std(totals):.3f}\n')
        
        global_interlocks = []
        for fx, spot_switch in enumerate(global_spot_switches):
            interlock_time = 0.0
            if spot_switch > np.median(global_spot_switches) + 5:
                this_interlock = spot_switch - np.median(global_spot_switches)
                interlock_time += this_interlock
                global_spot_switches[fx] = global_spot_switches[fx] - this_interlock
                print(f'  /!\ Possible spot switch interlock detected in fraction #{fx + 1}')
            if global_energy_switches[fx] > np.median(global_energy_switches) + 5:
                this_interlock = global_energy_switches[fx] - np.median(global_energy_switches)
                interlock_time += this_interlock
                global_energy_switches[fx] = global_energy_switches[fx] - this_interlock
                print(f'  /!\ Possible energy switch interlock detected in fraction #{fx + 1}')

            global_interlocks.append(interlock_time)
        
        global_drills = np.array(global_drills)
        global_spot_switches = np.array(global_spot_switches)
        global_energy_switches = np.array(global_energy_switches)
        global_interlocks = np.array(global_interlocks)

        ec = 'none'
        axs[0].bar(range(len(global_drills)), global_drills, label='Drill', color='yellow', edgecolor=ec)
        axs[0].bar(range(len(x_axis)), global_spot_switches, bottom=global_drills, label='Spot switch', color='tab:orange', edgecolor=ec)
        axs[0].bar(range(len(x_axis)), global_energy_switches, bottom=global_spot_switches + global_drills, label='Energy switch', color='tab:red', edgecolor=ec)
        axs[0].bar(range(len(x_axis)), global_interlocks, bottom=global_spot_switches + global_drills + global_energy_switches, label='Interlock', color='tab:purple', edgecolor=ec)
        axs[3].plot(range(len(x_axis)), global_drills, marker='o', color='black', markerfacecolor='yellow', label='Drill')
        axs[3].axhline(np.median(global_drills), ls='--', color='black', lw=1.0, zorder=-1, label='$\widetilde{x} \pm \sigma =$' + f' ({np.median(global_drills):.2f} $\pm$ {np.std(global_drills):.2f}) s')
        axs[2].plot(range(len(x_axis)), global_spot_switches, marker='o', color='black', markerfacecolor='tab:orange', label='Spot switch')  
        axs[2].axhline(np.median(global_spot_switches), ls='--', color='black', lw=1.0, zorder=-1, label='$\widetilde{x} \pm \sigma =$' + f' ({np.median(global_spot_switches):.2f} $\pm$ {np.std(global_spot_switches):.2f}) s')
        axs[1].plot(range(len(x_axis)), global_energy_switches, marker='o', color='black', markerfacecolor='tab:red', label='Energy switch')  
        axs[1].axhline(np.median(global_energy_switches), ls='--', color='black', lw=1.0, zorder=-1, label='$\widetilde{x} \pm \sigma =$' + f' ({np.median(global_energy_switches):.2f} $\pm$ {np.std(global_energy_switches):.2f}) s')

        for ax in axs:
            # ax.set_ylabel(f'Beam {beam_id}')
            # ax.set_yscale('log')
            ax.legend(loc='upper left')

        # print(f'\n\t\tMedian\tMin\tMax\tSigma')
        # print(f'Drill\t\t{np.median(global_drills):.3f}\t{np.min(global_drills):.3f}\t{np.max(global_drills):.3f}\t{np.std(global_drills):.3f}')
        # print(f'Spot switch\t{np.median(global_spot_switches):.3f}\t{np.min(global_spot_switches):.3f}\t{np.max(global_spot_switches):.3f}\t{np.std(global_spot_switches):.3f}')
        # print(f'Energy switch\t{np.median(global_energy_switches):.3f}\t{np.min(global_energy_switches):.3f}\t{np.max(global_energy_switches):.3f}\t{np.std(global_energy_switches):.3f}')
        # print('_____________________________________________________')
        # totals = global_drills + global_energy_switches + global_interlocks + global_spot_switches
        # print(f'Total beamtime\t{np.median(totals):.3f}\t{np.min(totals):.3f}\t{np.max(totals):.3f}\t{np.std(totals):.3f}\n')
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/{self.patient_id}_beam_timings_KS.png', dpi=300)        
        plt.show()


    def plan_creator(self, fraction, mode):
        if fraction == 'last':
            fx_list = [self.fraction_list[-1]]
        elif fraction == 'all':
            fx_list = [fx for fx in self.fraction_list]
        else: 
            root = Tk()
            fx_list = [filedialog.askdirectory(initialdir=self.logfile_dir).split('/')[-1]]
            root.destroy()
        
        for fx_id in fx_list:
            target_record = self.patient_record_df.loc[self.patient_record_df['FRACTION_ID'] == fx_id]
            target_tuning = self.patient_tuning_df.loc[self.patient_tuning_df['FRACTION_ID'] == fx_id]
            beam_list = [str(i) for i in target_record['BEAM_ID'].drop_duplicates()]
            plan_dcm, beam_ds = self.dicom_finder(fx_id, beam_list[0])
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
                    layer_energy = layer_df['LAYER_ENERGY(MeV)'].drop_duplicates().iloc[0]
                    plan_xy = plan_beam.IonControlPointSequence[layer_id * 2].ScanSpotPositionMap

                    # generate spot position list equal to DICOM tag --> [x1, y1, x2, y2, ..]
                    for log_spot in range(len(layer_df['X_POSITION(mm)'].to_list())):
                        layer_xy.append(layer_df['X_POSITION(mm)'][log_spot])
                        layer_xy.append(layer_df['Y_POSITION(mm)'][log_spot])
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
                        elif layer_id != 0 and spot_diff < 0:
                            print(f'  /x\ Critical: {abs(spot_diff)} less spots recorded than planned in layer-ID {layer_id}, fraction-ID {fx_id}\n      > Review log-files')

                        plan_beam.IonControlPointSequence[layer_id * 2].NumberOfScanSpotPositions = plan_beam.IonControlPointSequence[layer_id * 2 + 1].NumberOfScanSpotPositions = n_log_spots
                        plan_beam.IonControlPointSequence[layer_id * 2].ScanSpotPositionMap = plan_beam.IonControlPointSequence[layer_id * 2 + 1].ScanSpotPositionMap = layer_xy
                        plan_beam.IonControlPointSequence[layer_id * 2].CumulativeMetersetWeight = cumulative_mu
                        cumulative_mu += sum(layer_mu)
                        plan_beam.IonControlPointSequence[layer_id * 2 + 1].CumulativeMetersetWeight = cumulative_mu
                        plan_beam.IonControlPointSequence[layer_id * 2].ScanSpotMetersetWeights = layer_mu
                        plan_beam.IonControlPointSequence[layer_id * 2 + 1].ScanSpotMetersetWeights = [0.0 for _ in range(len(layer_mu))]
                        plan_beam.IonControlPointSequence[layer_id * 2].NominalBeamEnergy = plan_beam.IonControlPointSequence[layer_id * 2 + 1].NominalBeamEnergy = layer_energy
                    
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
            sop_uid = str(ds.SOPInstanceUID).split('.')
            sop_uid[-1] = '99999'
            new_sop_uid = '.'.join(sop_uid)
            ds.SOPInstanceUID = new_sop_uid
            if not ds.RTPlanName.__contains__('log'):
                ds.RTPlanName += f'_LOG_{fx_id}_{mode}'
            ds.RTPlanLabel = ds.RTPlanName
            ds.ReferringPhysicianName = 'Wolter^Lukas'
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

                    date = pd.to_datetime(fx_df.index[0]).date()
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
                self.patient_sss_df = pd.read_csv(os.path.join(self.df_destination, file), index_col='UNIQUE_INDEX', dtype={'BEAM_ID':str})
                single_stats_exists = True
            if delta_exists and single_stats_exists:
                re_init = input(f'''Re-initialize [y/N]? ''')
                if re_init == 'y':
                    break
                else:
                    return None
        
        if not delta_exists:
            print(f'  /x\ No deltaframe detected for patient_id {self.patient_id}, please run prepare_deltaframe()')
        
        print(f'\nInitializing spot statistics dataframe for patient-ID {self.patient_id}..')
        spots_data = {beam_id:[] for beam_id in self.patient_delta_df.BEAM_ID.drop_duplicates()}
        operation_columns = ['DELTA_X(mm)', 'DELTA_Y(mm)', 'DELTA_MU', 'DELTA_E(MeV)']

        # get spot counts per beam and fraction
        for fx_no, beams in enumerate(self.beam_list):
            for beam_id in beams:
                beam_df = self.patient_delta_df.loc[(self.patient_delta_df['FRACTION_ID'] == self.fraction_list[fx_no]) & (self.patient_delta_df['BEAM_ID'] == beam_id)]
                if not beam_df.empty: 
                    spots_data[beam_id].append(len(beam_df))
                else:
                    spots_data[beam_id].append(np.nan)
        
        to_concat = []
        for beam_id in spots_data.keys():
            beam_df = self.patient_delta_df.loc[self.patient_delta_df['BEAM_ID'] == beam_id]
            if beam_df.isnull().values.any() or beam_df.empty: continue

            # set reference fraction (first fx holding the most common spot count)
            most_common = max(set(spots_data[beam_id]), key=spots_data[beam_id].count)
            reference_fx = self.fraction_list[spots_data[beam_id].index(most_common)]
            print(f'  Referencing all spots for beam-ID {beam_id} to fx-ID {reference_fx}..')
            reference_df = beam_df.loc[beam_df['FRACTION_ID'] == reference_fx].reset_index()
            if reference_df.empty: 
                print('empty')
                continue

            beam_sss_df = pd.DataFrame()
            for layer_id in reference_df.LAYER_ID.drop_duplicates():
                # print(f'    Starting layer-ID {layer_id}')
                layer_df = reference_df.loc[reference_df.LAYER_ID == layer_id]
                fx_dfs = []
                for fx_id in self.fraction_list:
                    fx_df = beam_df.loc[(beam_df.FRACTION_ID == fx_id) & (beam_df.LAYER_ID == layer_id)]
                    if fx_df.empty: continue
                    
                    fx_df = fx_df.set_index('SPOT_ID')[operation_columns]

                    # if same spot count -> no sorting needed
                    if len(fx_df) == len(layer_df):
                        fx_dfs.append(fx_df)

                    # find and drop extra spot(s) from statistics
                    else:
                        layer_record = self.patient_record_df.loc[(self.patient_record_df.BEAM_ID == beam_id) & (self.patient_record_df.LAYER_ID == layer_id)] 
                        record_reference, record_current = layer_record.loc[layer_record.FRACTION_ID == reference_fx], layer_record.loc[layer_record.FRACTION_ID == fx_id]
                        record_reference['X_REF'], record_reference['Y_REF'] = record_reference['X_POSITION(mm)'].round(0), record_reference['Y_POSITION(mm)'].round(0)
                        record_current['X_CURRENT'], record_current['Y_CURRENT'] = record_current['X_POSITION(mm)'].round(0), record_current['Y_POSITION(mm)'].round(0)
                        record_reference.set_index('SPOT_ID', inplace=True), record_current.set_index('SPOT_ID', inplace=True)
                        while len(record_current) != len(record_reference):
                            compare_df = pd.concat([record_reference, record_current[['X_CURRENT', 'Y_CURRENT']]], axis=1)
                            indices = compare_df.loc[(abs(compare_df['X_REF'] - compare_df['X_CURRENT']) > 1) | (abs(compare_df['Y_REF'] - compare_df['Y_CURRENT']) > 1) | (compare_df['X_CURRENT'].apply(np.isnan)) | (compare_df['X_REF'].apply(np.isnan))].index
                            if len(record_current) > len(record_reference):
                                record_current.drop(record_current.index[indices[0]], inplace=True)
                                fx_df.drop(fx_df.index[indices[0]], inplace=True)
                            else:
                                record_reference.drop(record_reference.index[indices[0]], inplace=True)
                                fx_df.index = fx_df.index * 2 + 2
                                fx_df.loc[indices[0] * 2 + 1] = np.nan

                            fx_df = fx_df.sort_index().reset_index()
                            record_reference.reset_index(inplace=True), record_current.reset_index(inplace=True)
                        
                        fx_df = fx_df.drop(columns='SPOT_ID')
                        fx_df.index.name = 'SPOT_ID'
                        fx_dfs.append(fx_df)
                    
                layer_stacked = pd.concat(fx_dfs)
                layer_means, layer_stds = layer_stacked.groupby('SPOT_ID').mean(), layer_stacked.groupby('SPOT_ID').std()
                layer_means.rename(columns={col:col + '_MEAN' for col in operation_columns}, inplace=True)
                layer_stds.rename(columns={col:col + '_STD' for col in operation_columns}, inplace=True)
                sss_data = pd.concat([layer_means, layer_stds], axis=1)
                beam_sss_df = pd.concat([beam_sss_df, sss_data])

            beam_sss_df.reset_index(inplace=True)
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
                



if __name__ == '__main__':
    root_dir = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\Logfiles\converted'
    # root_dir = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\01_SpotShape\Logfiles_Spotshape_QA\converted'
    # root_dir = r'N:\fs4-HPRT\HPRT-Docs\Lukas\Logfile_Extraction\Logfiles'
    # root_dir = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\4D-PBS-LogFileBasedRecalc\Patient_dose_reconstruction\MOBILTest04_665914\Logfiles'
    # root_dir = r'/home/luke/Logfile_Extraction/1676348/Logfiles'
    # erroneous = [1230180, 1625909, 1627648, 1660835, 1698000, 1700535]
    # for errid in erroneous:
    #     dir = os.path.join(root_dir, str(errid))
    #     log = MachineLog(dir)
    #     log.prepare_dataframe()
    #     log.prepare_dataframe()

    ponaqua_qualified = [id.strip('\n') for id in open(r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\qualified_IDs.txt', 'r').readlines()]
    for id in ponaqua_qualified:
        # if id != "1632783": continue
        log = MachineLog(os.path.join(root_dir, id))
        log.prepare_sss_dataframe()
        # x_mean = log.patient_sss_df['DELTA_Y(mm)_MEAN']
        # x_std = log.patient_sss_df['DELTA_Y(mm)_STD']
        # plt.hist(x_mean, bins=100)
        # plt.hist(x_std, bins=50)
        # plt.show()
        # log.delta_dependencies()
        
        # df = log.patient_record_df.loc[(log.patient_record_df['FRACTION_ID'] == '20210617') & (log.patient_record_df['BEAM_ID'] == '2') & (log.patient_record_df['LAYER_ID'] == 5)]
        # log.prepare_dataframe()
        # log.prepare_deltaframe()
        # log.delta_dependencies()
        # log.plot_beam_layers()
    
    # patients = {}
    # print('Searching for log-file directories with existent plans..')
    # for root, dir, files in os.walk(root_dir):
    #     if dir.__contains__('Logfiles') and dir.__contains__('DeliveredPlans'):
    #         patient_id = root.split('\\')[-1]
    #         print(f'''  Found {patient_id}''')
    #         patients[patient_id] = os.path.join(root, 'Logfiles')
    # for patient_id, log_dir in patients.items():
    #     print(f'\n...STARTING PATIENT {patient_id}...\n')
    #     log = MachineLog(log_dir)
    #     log.prepare_dataframe()
    #     # log.prepare_deltaframe()
    #     # log.delta_dependencies()

    # log = MachineLog(root_dir)
    # df = log.patient_record_df
    # for fx in log.fraction_list:
    #     control = df.loc[(df['FRACTION_ID'] == fx) & (df['BEAM_ID'] == '2')]
    #     print(fx, control['LAYER_ID'].max())
    # df1 = df.loc[(df['FRACTION_ID'] == log.fraction_list[0]) & (df['BEAM_ID'] == '2')]
    # df2 = df.loc[(df['FRACTION_ID'] == '20200723') & (df['BEAM_ID'] == '2')]
    # for lid in range(50):
    #     dfslice1, dfslice2 = df1.loc[(df['LAYER_ID'] == lid)], df2.loc[(df['LAYER_ID'] == lid)]
    #     plt.scatter(dfslice1['X_POSITION(mm)'], dfslice1['Y_POSITION(mm)'], alpha=0.5)
    #     plt.scatter(dfslice2['X_POSITION(mm)'], dfslice2['Y_POSITION(mm)'], alpha=0.5)
    #     plt.show()

    # log.prepare_dataframe()
    # log.plot_beam_layers()
    # log.prepare_sss_dataframe()
    # log.prepare_qa_dataframe()
    # log.prepare_deltaframe()
    # log.delta_dependencies()
    # log.fractional_evolution(all=True)
    # log.plan_creator(fraction='all', mode='all')
    # log.delta_correlation_matrix(gtr_only=False)
    pass

else:
    print('>> Module', __name__, 'loaded')
    