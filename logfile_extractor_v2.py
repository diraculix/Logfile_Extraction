'''Encoding: UTF-8'''
__author__ = 'Lukas C. Wolter, OncoRay ZIK, Dresden, Germany'
__project__ = 'Logfile-based dose calculation & beam statistics'
__version__ = 2.0

from cProfile import label
import os
import sys
import pydicom
import pandas as pd
import numpy as np
from scipy import stats
from tkinter import Tk, filedialog
from matplotlib import pyplot as plt


class MachineLog:
    def __init__(self) -> None:
        valid_dir = False
        while not valid_dir:
            root = Tk()
            self.logfile_dir = filedialog.askdirectory(title='Select logfile root directory')  # open logfile root directory, which contains one dir per fraction
            root.destroy()
            if self.logfile_dir == '':
                sys.exit('No directory selected, exiting..')
            for index, element in enumerate(os.listdir(self.logfile_dir)):
                if not os.path.isdir(os.path.join(self.logfile_dir, element)):
                    print(f'''Chosen path '{self.logfile_dir}' may only contain directories (one per fraction). Please retry..''')
                    break
                elif index == len(os.listdir(self.logfile_dir)) - 1:
                    valid_dir = True
        
        self.df_destination = r'N:\fs4-HPRT\HPRT-Docs\Lukas\Logfile_Extraction\dataframes'  # TO BE CHANGED
        self.patient_record_df, self.patient_tuning_df = pd.DataFrame(), pd.DataFrame()  # initialize empty patient dataframes
        self.fraction_list = os.listdir(self.logfile_dir)
        self.num_fractions = len(self.fraction_list)
        self.beam_list = []
        for f in self.fraction_list:
            beams_in_frac = os.listdir(os.path.join(self.logfile_dir, f))
            self.beam_list.append(beams_in_frac)
        
        for fraction_no, fraction_id in enumerate(self.fraction_list):
            for beam_no, beam_id in enumerate(self.beam_list[fraction_no]):
                current_beam_path = os.path.join(self.logfile_dir, fraction_id, beam_id)
                os.chdir(current_beam_path)
                while len(os.listdir('.')) <= 3:  # log-file dirs may be nested irregularly
                    os.chdir(os.listdir('.')[0])
                
                for file in os.listdir('.'):
                    if file.__contains__('beam.'):
                        beam_file = file
                with open(beam_file, 'r') as beam_file:
                    for line in beam_file.readlines():
                            if line.__contains__('PatientId'):
                                self.patient_id = int(line.split('>')[1].split('<')[0])
        
        record_df_exists = False
        tuning_df_exists = False
        for dirpath, dirnames, filenames in os.walk(os.path.join(self.df_destination, '..')):  # sweep over root directory, read existing dataframes
            for fname in filenames:
                if fname.__contains__(f'{self.patient_id}_records') and fname.endswith('.csv'):
                    self.record_df_name = fname
                    print(f'''Found patient record dataframe '{self.record_df_name}', reading in..''')
                    self.patient_record_df = pd.read_csv(os.path.join(dirpath, fname), index_col='TIME', dtype={'BEAM_ID': str})
                    record_df_exists = True
                elif fname.__contains__(f'{self.patient_id}_tunings') and fname.endswith('.csv'):
                    self.tuning_df_name = fname
                    print(f'''Found patient tuning dataframe '{self.tuning_df_name}', reading in..''')
                    self.patient_tuning_df = pd.read_csv(os.path.join(dirpath, fname), index_col='TIME', dtype={'BEAM_ID': str})
                    tuning_df_exists = True
        
        if not record_df_exists or not tuning_df_exists:
            print(f'''\nUnable to locate patient record/tuning dataframes for patient-ID {self.patient_id}.\nPlease run MachineLog.prepare_dataframe()''')                
        
        os.chdir(self.logfile_dir)
        
    def summarize_beams(self):
        print(f'Total {self.num_fractions} fractions in {self.logfile_dir}\n')
        for index, fraction in enumerate(self.fraction_list):
            print(f'[{index + 1}] - {fraction}')
            for beam in self.beam_list[index]:
                print('\t', beam)
    
    def prepare_dataframe(self, sparse=True):  # disabling sparse will keep every 250us entry, enables single spot averaging, but results in larger dataframe
        if not self.patient_record_df.empty and not self.patient_tuning_df.empty:  # function call obsolete if record/tuning dataframes exist, re-init possible
            print('Already read existing patient record/tuning dataframes:')
            print(f'  {self.record_df_name}\n  {self.tuning_df_name}\n')
            re_init = input('Re-initialize patient dataframe [y/n]? ')
            if re_init != 'y':
                return None

        print(f'Initializing dataframe for patient-ID {self.patient_id}..')
        self.patient_record_df, self.patient_tuning_df = pd.DataFrame(), pd.DataFrame()  # overwrite stored df's
        for fraction_no, fraction_id in enumerate(self.fraction_list):
            for beam_no, beam_id in enumerate(self.beam_list[fraction_no]):
                current_beam_path = os.path.join(self.logfile_dir, fraction_id, beam_id)
                os.chdir(current_beam_path)
                while len(os.listdir('.')) <= 3:
                    os.chdir(os.listdir('.')[0])

                map_records, tunings, record_specifs, tuning_specifs = [], [], [], []  # get file lists
                for file in os.listdir('.'):
                    if file.__contains__('beam.'):
                        beam_file = file
                    elif file.__contains__('beam_config.'):
                        beam_config = file
                    elif file.__contains__('map_record') and file.__contains__('part'):
                        map_records.append(file)
                    elif file.__contains__('map_record') and file.__contains__('tuning'):
                        tunings.append(file)
                    elif file.__contains__('map_specif') and file.__contains__('part'):
                        record_specifs.append(file)
                    elif file.__contains__('map_specif') and file.__contains__('tuning'):
                        tuning_specifs.append(file)
                    
                if len(map_records) == 0 or len(tunings) == 0:
                    raise LookupError('No logfiles found!') 

                num_layers = max([int(fname.split('_')[2].split('_')[0]) for fname in map_records])

                with open(beam_file, 'r') as beam_file:  # draw beam specs from *_beam.csv
                    for line in beam_file.readlines():
                        if line.__contains__('GantryAngle'):
                            gantry_angle = float(line.split('>')[1].split('<')[0])
                        elif line.__contains__('pressure'):
                            pressure = float(line.split(',')[-1])
                        elif line.__contains__('temperature'):
                            temperature = float(line.split(',')[-1])
                        elif line.__contains__('doseCorrectionFactor'):
                            chamber_correction = float(line.split(',')[-1])

                with open(beam_config, 'r') as beam_config:  # draw machine parameters from *beam_config.csv
                    for line in beam_config.readlines():
                        if line.__contains__('SAD parameter'):
                            sad_x, sad_y = float(line.split(';')[-1].split(',')[1]), float(line.split(';')[-1].split(',')[0])  # coordinate system switch x <--> y
                        elif line.__contains__('distanceFromIcToIsocenter'):
                            ictoiso_x, ictoiso_y = float(line.split(';')[-1].split(',')[1]), float(line.split(';')[-1].split(',')[0])  # coordinate system switch x <--> y
                        elif line.__contains__('chargePerMUIC2'):
                            charge_per_mu = float(line.split(';')[-1])
                
                    ref_pressure, ref_temperature = 1013., 293.
                    correction_factor = chamber_correction * (ref_pressure * temperature) / (ref_temperature * pressure)  # classic chamber-specific * (p_0*T) / (T_0*p)
                    beam_file.close()
                
                layer_exceptions = []
                for layer_id in range(num_layers):
                    to_do_layers, to_do_tunings = [], []
                    no_exceptions = True
                    for record_file in map_records:  # actual (processed) log-file analysis
                        if int(record_file.split('_')[2].split('_')[0]) == layer_id:
                            record_file_df = pd.read_csv(record_file, delimiter=',', skiprows=10, skipfooter=11, engine='python')
                            record_file_df['TIME'] = pd.to_datetime(record_file_df['TIME'])     # datetime index --> chronological order
                            record_file_df.index = record_file_df['TIME']                                 
                            charge_col = pd.Series(record_file_df['DOSE_PRIM(C)'])              # ion dose [C], will be translated to MU later
                            record_file_df = record_file_df.loc[:, :'Y_POSITION(mm)']           # slice dataframe, drop redundant columns
                            record_file_df['DOSE_PRIM(C)'] = charge_col
                            record_file_df.drop(columns=['TIME'], inplace=True)
                            record_file_df = record_file_df.loc[(record_file_df['X_WIDTH(mm)'] != -10000.0) & (record_file_df['Y_WIDTH(mm)'] != -10000.0) & (record_file_df['X_POSITION(mm)'] != -10000.0) & (record_file_df['Y_POSITION(mm)'] != -10000.0)]
                                                                                                # drop all rows where no spot is targeted (machine stuff going on)
                            if record_file_df.empty:  # in case no usable information left in log-file, skip these occasions (very unlikely)
                                no_exceptions = False
                                continue
                            
                            current_spot_submap = record_file_df['SUBMAP_NUMBER'].min()
                            current_spot_id = 0
                            record_file_df['SPOT_ID'] = 0
                            while current_spot_submap <= record_file_df['SUBMAP_NUMBER'].max(): # sweep over all spots (submaps), SUBMAP_NUMBER is locked whenever a spot is active
                                record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['SPOT_ID']] = current_spot_id  
                                                                                                                                # assign new column 'SPOT_ID'
                                if sparse:
                                    accumulated_charge = record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['DOSE_PRIM(C)']].abs().sum().mean()
                                    record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['ACC_CHARGE(C)']] = accumulated_charge
                                                                                                                                # accumulate released charge and assign to new column

                                current_spot_submap = record_file_df.loc[record_file_df['SUBMAP_NUMBER'] > current_spot_submap]['SUBMAP_NUMBER'].min()
                                current_spot_id += 1

                            if sparse:
                                record_file_df.drop(columns=['DOSE_PRIM(C)'], inplace=True)
                                record_file_df.drop_duplicates(subset=['SUBMAP_NUMBER'], keep='last', inplace=True)  # keep only last entries for each spot
                            
                            for specif_file in record_specifs:  # draw machine parameters from *map_specif*.csv
                                if int(specif_file.split('_')[2].split('_')[0]) == layer_id:
                                    with open(specif_file, 'r') as specif_file:
                                        line = specif_file.readlines()[3]
                                        ic_offset_x, ic_offset_y = float(line.split(',')[3]), float(line.split(',')[2])  # coordinate system switch x <--> y
                                        del line
                            
                            new_x_series = (pd.Series(record_file_df['Y_POSITION(mm)']) - ic_offset_x) * sad_x / (sad_x - ictoiso_x)  # coordinate system transform iba <-> raystation (x <-> y)
                            new_y_series = (pd.Series(record_file_df['X_POSITION(mm)']) - ic_offset_y) * sad_y / (sad_y - ictoiso_y)
                            record_file_df['X_POSITION(mm)'], record_file_df['Y_POSITION(mm)'] = new_x_series, new_y_series
                            del new_x_series, new_y_series
                        
                            if sparse:
                                record_file_df['MU'] = record_file_df['ACC_CHARGE(C)'] / correction_factor / charge_per_mu # calculate accumulated MUs and assign new column 'MU'
                                record_file_df.drop(columns=['ACC_CHARGE(C)'], inplace=True)
                            else:
                                record_file_df['MU'] = abs(record_file_df['DOSE_PRIM(C)'] / correction_factor / charge_per_mu)  # calculate single MUs and assign new column 'MU'
                                record_file_df.drop(columns=['DOSE_PRIM(C)'], inplace=True)

                            record_file_df.reindex()  # make sure modified layer df is consistent with indexing
                            to_do_layers.append(record_file_df)

                    for tuning_file in tunings:  # do the same for all tuning files
                        if int(tuning_file.split('_')[2].split('_')[0]) == layer_id:
                            try:
                                tuning_file_df = pd.read_csv(tuning_file, delimiter=',', skiprows=10, skipfooter=11, engine='python')
                                tuning_file_df['TIME'] = pd.to_datetime(tuning_file_df['TIME'])
                                tuning_file_df.index = tuning_file_df['TIME']
                                charge_col = pd.Series(tuning_file_df['DOSE_PRIM(C)'])
                                tuning_file_df = tuning_file_df.loc[:, :'Y_POSITION(mm)']
                                tuning_file_df['DOSE_PRIM(C)'] = charge_col
                                tuning_file_df.drop(columns=['TIME'], inplace=True)
                                tuning_file_df = tuning_file_df.loc[(tuning_file_df['X_WIDTH(mm)'] != -10000.0) & (tuning_file_df['Y_WIDTH(mm)'] != -10000.0) & (tuning_file_df['X_POSITION(mm)'] != -10000.0) & (tuning_file_df['Y_POSITION(mm)'] != -10000.0)]
                            except KeyError:
                                print(f'''\n  /!\ Key error occured while handling tuning file '{tuning_file}' dataframe (layer {str(layer_id).zfill(2)})''')
                                continue

                            if tuning_file_df.empty:  # in case of low-weighted tuning spots, no usable information will be left in log-file, skip these occasions
                                no_exceptions = False
                                continue

                            current_spot_submap = tuning_file_df['SUBMAP_NUMBER'].min()
                            current_spot_id = 0
                            tuning_file_df['SPOT_ID'] = 0
                            tuning_file_df.reindex()
                            while current_spot_submap <= tuning_file_df['SUBMAP_NUMBER'].max():
                                tuning_file_df.loc[tuning_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['SPOT_ID']] = current_spot_id

                                if sparse:
                                    accumulated_charge = tuning_file_df.loc[tuning_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['DOSE_PRIM(C)']].abs().sum().mean()
                                    tuning_file_df.loc[tuning_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['ACC_CHARGE(C)']] = accumulated_charge

                                current_spot_submap = tuning_file_df.loc[tuning_file_df['SUBMAP_NUMBER'] > current_spot_submap]['SUBMAP_NUMBER'].min()
                                current_spot_id += 1
                            
                            if sparse:
                                tuning_file_df.drop(columns=['DOSE_PRIM(C)'], inplace=True)
                                tuning_file_df.drop_duplicates(subset=['SUBMAP_NUMBER'], keep='last', inplace=True)
                            
                            for specif_file in tuning_specifs:
                                if int(specif_file.split('_')[2].split('_')[0]) == layer_id:
                                    with open(specif_file, 'r') as specif_file:
                                        line = specif_file.readlines()[3]
                                        ic_offset_x, ic_offset_y = float(line.split(',')[2]), float(line.split(',')[3])
                                        del line                

                            new_x_series = (pd.Series(tuning_file_df['Y_POSITION(mm)']) - ic_offset_x) * sad_x / (sad_x - ictoiso_x)  # coordinate system transform iba <-> raystation (x <-> y)
                            new_y_series = (pd.Series(tuning_file_df['X_POSITION(mm)']) - ic_offset_y) * sad_y / (sad_y - ictoiso_y)
                            tuning_file_df['X_POSITION(mm)'], tuning_file_df['Y_POSITION(mm)'] = new_x_series, new_y_series
                            del new_x_series, new_y_series            
                            
                            if sparse:
                                tuning_file_df['MU'] = tuning_file_df['ACC_CHARGE(C)'] / correction_factor / charge_per_mu
                                tuning_file_df.drop(columns=['ACC_CHARGE(C)'], inplace=True)
                            else:
                                tuning_file_df['MU'] = abs(tuning_file_df['DOSE_PRIM(C)'] / correction_factor / charge_per_mu)
                                tuning_file_df.drop(columns=['DOSE_PRIM(C)'], inplace=True)

                            tuning_file_df.reindex()
                            to_do_tunings.append(tuning_file_df)
                    
                    for i in range(len(to_do_tunings)):  # in case of multiple layer parts: enable continuous spot indexing
                        if i > 0:  
                            to_do_tunings[i]['SPOT_ID'] += (to_do_tunings[i - 1]['SPOT_ID'].max() + 1)
                    for j in range(len(to_do_layers)):
                        if j > 0:
                            to_do_layers[j]['SPOT_ID'] += (to_do_layers[j - 1]['SPOT_ID'].max() + 1)
                    
                    layer_df = pd.concat(to_do_layers)  # concatenate layers, assign additional columns
                    layer_df['LAYER_ID'] = layer_id
                    layer_df['TOTAL_LAYERS'] = num_layers
                    layer_df['BEAM_ID'] = beam_id
                    layer_df['GANTRY_ANGLE'] = gantry_angle
                    layer_df['FRACTION_ID'] = fraction_id
                    # layer_df['PATIENT_ID'] = self.patient_id
                    layer_df.drop(columns=['SUBMAP_NUMBER'], inplace=True)
                    tuning_df = pd.concat(to_do_tunings)
                    tuning_df['LAYER_ID'] = layer_id
                    tuning_df['TOTAL_LAYERS'] = num_layers
                    tuning_df['BEAM_ID'] = beam_id
                    tuning_df['GANTRY_ANGLE'] = gantry_angle
                    tuning_df['FRACTION_ID'] = fraction_id
                    # tuning_df['PATIENT_ID'] = self.patient_id
                    tuning_df.drop(columns=['SUBMAP_NUMBER'], inplace=True)
                    del to_do_layers, to_do_tunings
                    
                    if self.patient_record_df.empty:
                        self.patient_record_df = layer_df
                        self.patient_tuning_df = tuning_df
                    else:
                        self.patient_record_df = pd.concat([self.patient_record_df, layer_df], sort=True)    
                        self.patient_tuning_df = pd.concat([self.patient_tuning_df, tuning_df], sort=True)
                    
                    del layer_df, tuning_df
                    
                    if no_exceptions:
                        char = '#'
                    else:
                        char = '_'

                    if layer_id == (num_layers - 1):  # progress visualization
                        if no_exceptions:
                            print('  ', '[' + (layer_id + 1) * char + (num_layers - layer_id - 1) * '-' + ']', end=f' Beam {beam_id} complete\n')
                        else:    
                            print('  ', '[' + (layer_id + 1) * char + (num_layers - layer_id - 1) * '-' + ']', end=f' Beam {beam_id} complete (empty dataframe exception in layer(s) {layer_exceptions})\n')
                    else:
                        print('  ', '[' + (layer_id + 1) * char + (num_layers - layer_id - 1) * '-' + ']', end=f' Layer {str(layer_id + 1).zfill(2)}/{str(num_layers).zfill(2)}\r')
                    
                    no_exceptions = True

                os.chdir(self.logfile_dir) 

            print(f'  ..Fraction {str(fraction_no + 1).zfill(2)}/{str(self.num_fractions).zfill(2)} complete..') 

        os.chdir(self.df_destination)
        print(f'''  ..Writing dataframe to '{self.df_destination}' as .CSV.. ''')
        while True:   
            try:
                self.patient_record_df.to_csv(f'patient_{self.patient_id}_records_data.csv')
                self.patient_tuning_df.to_csv(f'patient_{self.patient_id}_tunings_data.csv')
                break
            except PermissionError:
                input('  Permission denied, close target file and press ENTER.. ')

        print('Complete')
    
    def plot_beam_layers(self):     # For all layers in one beam, averaged over all fractions:
                                    # derive scatter plot of (x,y)-positions from logfile, 
                                    # compare vs. planned positions extracted from RP*.dcm
        if self.patient_record_df.empty or self.patient_tuning_df.empty:
            print(f'''\nUnable to locate patient record/tuning dataframes for patient-ID {self.patient_id}. Calling MachineLog.prepare_dataframe()..''')                
            self.prepare_dataframe()

        beam_list = self.patient_record_df['BEAM_ID'].drop_duplicates()
        indices = beam_list.index.to_list()

        print('\nSelect beam key for layer-wise spot plotting:\n')
        print('Key\tBeam-ID\t\tGantry Angle\tLayers\tFrom Fraction\n')  # choose beam
        for choice, (index, beam) in enumerate(zip(indices, beam_list)):
            num_layers = int(self.patient_record_df['TOTAL_LAYERS'][index].mean())
            try:
                fraction = int(self.patient_record_df['FRACTION_ID'][index].mean())
            except:
                fraction = int(self.patient_record_df['FRACTION_ID'][index])
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

        print('\nTrying to auto-locate plan .dcm file..')  # read RT plan dicom via filedialog
        for path, dirnames, filenames in os.walk(os.path.join(self.logfile_dir, '..')):
            for fname in filenames:
                if fname.__contains__('RP') and fname.endswith('.dcm'):
                    ds = pydicom.read_file(os.path.join(path, fname))
                    for i, beam in enumerate(ds.IonBeamSequence):
                        if beam.BeamName == beam_id or beam.BeamDescription == beam_id and float(beam.IonControlPointSequence[0].GantryAngle) == self.patient_record_df.loc[self.patient_record_df['BEAM_ID'] == beam_id]['GANTRY_ANGLE'].mean():
                            if 'dcm_beam_id' and 'dcm_target' in locals():
                                prompt = input(f'''  <!> Conflict: file '{fname}' also contains matching beam-ID {ds.IonBeamSequence[i].BeamName} ({ds.IonBeamSequence[i].BeamDescription}), keep previous[k] or replace[r]? ''')
                                if prompt == 'r':
                                    continue
                                else:
                                    print('  Scanning for conflicts.. ', end='\r')
                                    break
                            
                            dcm_target, dcm_beam_id = ds, i  # locate beam according to previous selection
                            print(f'''  Detected matching plan '{fname}' for beam-ID {beam_id} with patient-ID {ds.PatientID}''')
                            print('  Scanning for conflicts.. ', end='\r')
        
        while True:
            try:
                dcm_beam = dcm_target.IonBeamSequence[dcm_beam_id]
                print('  Scanning for conflicts.. OK')
                break
            except UnboundLocalError:
                print('Failed to auto-locate plan .dcm file. Please choose file manually..')
                root = Tk()
                try:
                    dcm_target = pydicom.read_file(filedialog.askopenfilename())
                except FileNotFoundError:
                    print('ERROR: File not found or none selected, try again..')
                root.destroy()
        
        del dcm_target, dcm_beam_id
        print('\nGenerating layer plot..')
        fig, axs = plt.subplots(6, 8, sharex=True, sharey=True, figsize=(24, 24 * 6/8), dpi=150)  # initiate matrix-like layer plot
        ax0 = fig.add_subplot(111, frameon=False)
        fig.subplots_adjust(hspace=0.0, wspace=0.0)
        axs = axs.ravel()

        pos_x_delta, pos_y_delta, pos_abs_delta = [], [], []
        
        for layer_id in scope_record_df['LAYER_ID'].drop_duplicates():
            x_positions, y_positions = [], []
            x_tunings, y_tunings, = [], []
            dcm_layer = dcm_beam.IonControlPointSequence[layer_id * 2]
            plan_spotmap = dcm_layer.ScanSpotPositionMap
            plan_x_positions, plan_y_positions = [], []
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
            xy_positions_sorted = [(0.0, 0.0) for i in range(max([len(spot_points_log), len(spot_points_plan)]))]
            for xy_log in spot_points_log:
                x_differences, y_differences = [], []
                distances = []
                for xy_plan in spot_points_plan:  # calculate position difference and euclidian distance for every spot pair, find minimum to reindex (return sorted spots)
                    diff = np.array(xy_plan) - np.array(xy_log)
                    sqdist = abs((diff).dot(diff))
                    distances.append(sqdist)
                    x_differences.append(diff[0]), y_differences.append(diff[1])

                min_dist = min(distances)
                pos_abs_delta.append(min_dist)
                index = distances.index(min_dist)
                pos_x_delta.append(x_differences[index]), pos_y_delta.append(y_differences[index])
                xy_positions_sorted[index] = xy_log
            
            for element in xy_positions_sorted:
                if element == (0.0, 0.0):
                    xy_positions_sorted.remove(element)

            axs[layer_id].plot(plan_x_positions, plan_y_positions, marker='x', linestyle='-', markersize=2.0, markeredgewidth=0.2, linewidth=0.2, label='Planned')
            axs[layer_id].plot(*zip(*xy_positions_sorted), marker='o', markerfacecolor='None', linestyle='-', color='black', markersize=2.0, markeredgewidth=0.2, linewidth=0.2, label='Log-file sorted')
            axs[layer_id].plot(x_positions, y_positions, marker='o', markerfacecolor='None', linestyle='--', color='black', markersize=2.0, markeredgewidth=0.2, linewidth=0.2, label='Log-file')
            axs[layer_id].plot(x_tunings, y_tunings, marker='o', markerfacecolor='None', linestyle='None', markersize=2.0, markeredgewidth=0.2, color='limegreen')
            axs[layer_id].annotate(f'Layer #{str(layer_id + 1).zfill(2)} | $\Delta$ = {abs(len(plan_x_positions) - len(x_positions))}', xy=(1.0, 1.0), xycoords='axes points', fontsize=8)
            if len(x_tunings) > 1:
                axs[layer_id].annotate(f'Layer #{str(layer_id + 1).zfill(2)} | $\Delta$ = {abs(len(plan_x_positions) - len(x_positions))} | $n_t$ > 1', xy=(1.0, 1.0), xycoords='axes points', fontsize=8) 
            axs[layer_id].legend(loc='upper right', fontsize=8)
        plt.suptitle(f'Spot Positions for Patient-ID {self.patient_id}, Beam-ID: {beam_id}', fontweight='bold', y=0.9)
        ax0.set_xlabel('X [mm]', fontweight='bold', labelpad=10)
        ax0.set_ylabel('Y [mm]', fontweight='bold', labelpad=10)
        plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
        output_dir = r'N:\fs4-HPRT\HPRT-Docs\Lukas\Logfile_Extraction\output'  # TO BE CHANGED
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        while True:
            try:
                plt.savefig(f'{output_dir}/{self.patient_id}_{beam_id}_spots.pdf')
                break
            except PermissionError:
                input('  Permission denied, close target file and press ENTER.. ')
        plt.close()
        plt.clf()
        
        print('Generating spot position histograms..')
        fig, axs = plt.subplots(1, 2, figsize=(10, 5), dpi=150)
        fig.subplots_adjust(hspace=4.0)
        x_bins, y_bins = int(max(pos_x_delta) - min(pos_x_delta)) * 10, int(max(pos_y_delta) - min(pos_y_delta) * 10)
        axs[0].hist(pos_x_delta, bins=x_bins, color='tab:blue', alpha=0.5, edgecolor='black', linewidth=1.0, label='$\Delta x$')
        axs[0].hist(pos_y_delta, bins=y_bins, color='tab:orange', alpha=0.5, edgecolor='black', linewidth=1.0, label='$\Delta y$')
        axs[0].set_ylabel('Counts')
        axs[0].set_xlabel('Plan - Log-file [mm]')
        axs[0].legend()
        axs[1].hist(pos_abs_delta, bins=int(max(pos_abs_delta) - min(pos_abs_delta)) * 10, density=True, color='white', edgecolor='black', linewidth=1.0, label='Abs. distance')
        axs[1].set_xlabel('$\Delta(x, y)$ [mm]')
        axs[1].set_ylabel('Counts')
        axs[1].legend()
        plt.suptitle('Spot Position Histogram Plan vs. Log-file', fontweight='bold', y=0.95)
        while True:
            try:
                plt.savefig(f'{output_dir}/{self.patient_id}_{beam_id}_pos_stats.pdf')
                break
            except PermissionError:
                input('  Permission denied, close target file and press ENTER.. ')
        plt.clf()

    def plot_spot_statistics(self):     # For one beam over all fractions:
                                        # derive histograms and MAE-evolution by comparing log vs. plan spot positions and MU's
        if self.patient_record_df.empty or self.patient_tuning_df.empty:
            print(f'''\nUnable to locate patient record/tuning dataframes for patient-ID {self.patient_id}. Calling MachineLog.prepare_dataframe()..''')                
            self.prepare_dataframe()

        beam_list = [str(i) for i in self.patient_record_df['BEAM_ID'].drop_duplicates()]

        print('\nTrying to auto-locate patient plan dicoms..')  # read RT plan dicom via filedialog
        patient_dicoms = []
        for path, dirnames, filenames in os.walk(os.path.join(self.logfile_dir, '..')):
            for fname in filenames:
                if fname.__contains__('RP') and fname.endswith('.dcm'):
                    ds = pydicom.read_file(os.path.join(path, fname))
                    for i, beam in enumerate(ds.IonBeamSequence):
                        if beam.BeamName in beam_list or beam.BeamDescription in beam_list and float(beam.IonControlPointSequence[0].GantryAngle) in self.patient_record_df['GANTRY_ANGLE'].drop_duplicates():
                            patient_dicoms.append(os.path.join(path, fname))

        print('Collecting data..')
        total_beam_maes, layer_axis = [], []  # mean absolute errors
        for fraction_id in self.fraction_list:
            for beam_id in beam_list:
                scope_record_df = self.patient_record_df.loc[(self.patient_record_df['BEAM_ID'] == beam_id) & (self.patient_record_df['FRACTION_ID'] == int(fraction_id))]  # slice currently needed portion from patient df
                scope_tuning_df = self.patient_tuning_df.loc[(self.patient_tuning_df['BEAM_ID'] == beam_id) & (self.patient_tuning_df['FRACTION_ID'] == int(fraction_id))]

                for dcm in patient_dicoms:
                    ds = pydicom.read_file(dcm)
                    for i, beam in enumerate(ds.IonBeamSequence):
                        if beam.BeamName == beam_id or beam.BeamDescription == beam_id and float(beam.IonControlPointSequence[0].GantryAngle) == scope_record_df.loc[scope_record_df['BEAM_ID'] == beam_id]['GANTRY_ANGLE'].mean():
                            dcm_beam = ds.IonBeamSequence[i]

                for layer_id in scope_record_df['LAYER_ID'].drop_duplicates():
                    x_positions, y_positions = [], []
                    x_tunings, y_tunings, = [], []
                    dcm_layer = dcm_beam.IonControlPointSequence[layer_id * 2]
                    plan_spotmap = dcm_layer.ScanSpotPositionMap
                    plan_x_positions, plan_y_positions = [], []
                    layer_abs_x_diffs, layer_abs_y_diffs = [], []
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
                        layer_abs_x_diffs.append(abs(x_differences[index])), layer_abs_y_diffs.append(abs(y_differences[index]))
                    
                    total_beam_maes.append((np.mean(layer_abs_x_diffs) + np.mean(layer_abs_y_diffs)) / 2)
                    layer_axis.append(str(layer_id))
            break
        
        print('Generating layer MAE evolution plot..')
        plt.plot(total_beam_maes, label='layer mean absolute error')
        plt.legend()
        plt.show()


if __name__ == '__main__':
    log = MachineLog()
    # log.prepare_dataframe()
    # log.summarize_beams()
    log.plot_beam_layers()    
    log.plot_spot_statistics()
