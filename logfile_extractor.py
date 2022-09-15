'''Encoding: UTF-8'''
__author__ = 'Lukas C. Wolter, OncoRay ZIK, Dresden, Germany'
__project__ = 'Logfile-based dose calculation & beam statistics'
__version__ = 5.0

import os
import sys
import pydicom
import pandas as pd
import numpy as np
from tkinter import Tk, filedialog
from matplotlib import pyplot as plt


# output_dir = r'N:\fs4-HPRT\HPRT-Docs\Lukas\Logfile_Extraction\output'  # TO BE CHANGED
output_dir = r'C:\Users\lukas\Documents\OncoRay HPRT\Logfile_Extraction_mobile\output'


class MachineLog:
    def __init__(self):
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
        
        # self.df_destination = r'N:\fs4-HPRT\HPRT-Docs\Lukas\Logfile_Extraction\dataframes'  # TO BE CHANGED
        self.df_destination = r'C:\Users\lukas\Documents\OncoRay HPRT\Logfile_Extraction_mobile\dataframes'
        self.fraction_list = os.listdir(self.logfile_dir)
        self.num_fractions = len(self.fraction_list)
        self.beam_list = []
        for f in self.fraction_list:
            beams_in_frac = []
            for dir in os.listdir(os.path.join(self.logfile_dir, f)):
                if os.path.isdir(os.path.join(self.logfile_dir, f, dir)):
                    beams_in_frac.append(dir)
            self.beam_list.append(beams_in_frac)

        self.missing_beams = []
        for fraction_no, fraction_id in enumerate(self.fraction_list):
            for beam_no, beam_id in enumerate(self.beam_list[fraction_no]):
                current_beam_path = os.path.join(self.logfile_dir, fraction_id, beam_id)
                os.chdir(current_beam_path)
                valid_beam = True
                while len(os.listdir('.')) <= 3:  # log-file dirs may be nested irregularly
                    try:
                        os.chdir(os.listdir('.')[0])
                    except:
                        print(f'  /x\ No log-files detected for beam {beam_id} in fraction {fraction_id}, continuing..')
                        self.missing_beams.append((fraction_id, beam_id))
                        valid_beam = False
                        break
                
                if not valid_beam:
                    continue

                for file in os.listdir('.'):
                    if file.__contains__('beam.'):
                        beam_file = file
                with open(beam_file, 'r') as beam_file:
                    for line in beam_file.readlines():
                            if line.__contains__('PatientId'):
                                self.patient_id = int(line.split('>')[1].split('<')[0])
        
        record_df_exists, tuning_df_exists = False, False
        for dirpath, dirnames, filenames in os.walk(os.path.join(self.df_destination, '..')):  # sweep over root directory, read existing dataframes
            for fname in filenames:
                if fname.__contains__(f'{self.patient_id}_records') and fname.endswith('.csv'):
                    self.record_df_name = fname
                    print(f'''Found patient record dataframe '{self.record_df_name}', reading in..''')
                    self.patient_record_df = pd.read_csv(os.path.join(dirpath, fname), index_col='TIME', dtype={'BEAM_ID':str, 'FRACTION_ID':str})
                    record_df_exists = True
                elif fname.__contains__(f'{self.patient_id}_tunings') and fname.endswith('.csv'):
                    self.tuning_df_name = fname
                    print(f'''Found patient tuning dataframe '{self.tuning_df_name}', reading in..''')
                    self.patient_tuning_df = pd.read_csv(os.path.join(dirpath, fname), index_col='TIME', dtype={'BEAM_ID':str, 'FRACTION_ID':str})
                    tuning_df_exists = True
        
        if not record_df_exists or not tuning_df_exists:
            print(f'''\nUnable to locate patient record/tuning dataframes for patient-ID {self.patient_id}. Calling prepare_dataframe()..''')
            self.patient_record_df, self.patient_tuning_df = pd.DataFrame(), pd.DataFrame()
            self.prepare_dataframe()
        
        os.chdir(self.logfile_dir)
    

    def prepare_dataframe(self):  # disabling sparse will keep every 250us entry, results in larger dataframes (deprecated)
        if not self.patient_record_df.empty and not self.patient_tuning_df.empty:  # function call obsolete if record/tuning dataframes exist, re-init possible
            print('Already read existing patient record/tuning dataframes:')
            print(f'  {self.record_df_name}\n  {self.tuning_df_name}\n')
            re_init = input('Re-initialize patient dataframe [y/n]? ')
            if re_init != 'y':
                return None

        print(f'Initializing dataframe for patient-ID {self.patient_id}..')
        self.patient_record_df, self.patient_tuning_df = pd.DataFrame(), pd.DataFrame()  # overwrite stored df's
        for fraction_no, fraction_id in enumerate(self.fraction_list):
            for beam_id in self.beam_list[fraction_no]:
                if (fraction_id, beam_id) in self.missing_beams:
                    print(f'  /!\ Skipping missing beam {beam_id} in fraction {fraction_id}')
                    continue

                current_beam_path = os.path.join(self.logfile_dir, fraction_id, beam_id)
                os.chdir(current_beam_path)
                while len(os.listdir('.')) <= 3:
                    try:
                        os.chdir(os.listdir('.')[0])  # navigate through nested logfile dir structure (possibly risky)
                    except OSError:
                        print(f'  /!\ No directory to enter in {os.getcwd()}, staying here..')
                        break

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
                    raise LookupError('  /x\ No logfiles found!') 

                num_layers = max([int(fname.split('_')[2].split('_')[0]) + 1 for fname in map_records])

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

                    beam_file.close()

                with open(beam_config, 'r') as beam_config:  # draw machine parameters from *beam_config.csv
                    for line in beam_config.readlines():
                        if line.__contains__('SAD parameter'):
                            sad_x, sad_y = float(line.split(';')[-1].split(',')[1]), float(line.split(';')[-1].split(',')[0])  # coordinate system switch x <--> y
                        elif line.__contains__('distanceFromIcToIsocenter'):
                            ictoiso_x, ictoiso_y = float(line.split(';')[-1].split(',')[1]), float(line.split(';')[-1].split(',')[0])  # coordinate system switch x <--> y
                        elif line.__contains__('chargePerMUIC2'):
                            charge_per_mu = float(line.split(';')[-1])
                        elif line.__contains__('Nozzle WET polynomial'):
                            nozzle_wet_poly = [float(x) for x in line.split(';')[-1].split(',')]
                
                    ref_pressure, ref_temperature = 1013., 293.  # [hPa, K] standard reference, can also be read from same file
                    correction_factor = (1 / chamber_correction) * (ref_pressure * temperature) / (ref_temperature * pressure)  # why inverse chamber correction?
                    beam_config.close()
                
                # Source: IBA Particle Therapy 08/22 (Jozef Bokor), universal nozzle WET polynomial coefficients
                iba_gtr2_poly = [0.001684756748152, -0.00490089228886989, 0.561372013469097, 3.46404838890297]
                
                layer_exceptions = []
                for layer_id in range(num_layers):
                    to_do_layers, to_do_tunings = [], []
                    no_exceptions = True
                    for record_file in map_records:  # actual (processed) log-file analysis
                        if int(record_file.split('_')[2].split('_')[0]) == layer_id:
                            try:
                                record_file_df = pd.read_csv(record_file, delimiter=',', skiprows=10, skipfooter=11, engine='python')
                            except:
                                print('Read CSV error:', fraction_id, beam_id, record_file)
                                continue
                            try:
                                record_file_df['TIME'] = pd.to_datetime(record_file_df['TIME'])     # datetime index --> chronological order
                                record_file_df.index = record_file_df['TIME']                              
                                charge_col = pd.Series(record_file_df['DOSE_PRIM(C)'])              # ion dose [C], to be converted in MU
                                record_file_df = record_file_df.loc[:, :'Y_POSITION(mm)']           # slice dataframe, drop redundant columns
                                record_file_df['DOSE_PRIM(C)'] = charge_col
                                record_file_df.drop(columns=['TIME'], inplace=True)
                                record_file_df.drop(record_file_df[record_file_df['SUBMAP_NUMBER'] < 0].index, inplace=True)
                                record_file_df = record_file_df[record_file_df.groupby('SUBMAP_NUMBER')['SUBMAP_NUMBER'].transform('count') > 1]  # drop all rows without plan-relevant data
                                record_file_df.drop(columns=['X_WIDTH(mm)', 'Y_WIDTH(mm)'], inplace=True)  # not needed for plan, uniform spot sizes in beam model
                            except:  # unlikely event of unusable information in log-file (possible if split into parts)
                                no_exceptions = False
                                layer_exceptions.append(layer_id)
                                continue
                            
                            current_spot_submap = record_file_df['SUBMAP_NUMBER'].min()
                            current_spot_id = 0
                            record_file_df['SPOT_ID'] = 0
                            record_file_df.reindex()
                            while current_spot_submap <= record_file_df['SUBMAP_NUMBER'].max(): # sweep over all spots (submaps), SUBMAP_NUMBER is locked whenever a spot is active
                                record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['SPOT_ID']] = current_spot_id  # assign new column
                                spot_drill_time = len(record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == current_spot_submap]) * 0.25  # in [ms]
                                record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['DRILL_TIME(ms)']] = spot_drill_time  # assign new column
                                accumulated_charge = record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['DOSE_PRIM(C)']].abs().sum().mean()  # accumulate charge released per spot
                                record_file_df.loc[record_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['ACC_CHARGE(C)']] = accumulated_charge  # assign new column

                                current_spot_submap = record_file_df.loc[record_file_df['SUBMAP_NUMBER'] > current_spot_submap]['SUBMAP_NUMBER'].min()
                                current_spot_id += 1

                            record_file_df.drop(columns=['DOSE_PRIM(C)'], inplace=True)
                            record_file_df.drop_duplicates(subset=['SUBMAP_NUMBER'], keep='last', inplace=True)  # keep only last entries for each spot (most accurate)
                            record_file_df = record_file_df.loc[(record_file_df['X_POSITION(mm)'] != -10000.0) & (record_file_df['Y_POSITION(mm)'] != -10000.0)]  # drop redundant rows
                            
                            for specif_file in record_specifs:  # draw machine parameters from *map_specif*.csv
                                if int(specif_file.split('_')[2].split('_')[0]) == layer_id:
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
                            
                            # coordinate system transform iba <-> raystation (x <-> y)
                            new_x_series = (pd.Series(record_file_df['Y_POSITION(mm)']) - ic_offset_x) * sad_x / (sad_x - ictoiso_x) 
                            new_y_series = (pd.Series(record_file_df['X_POSITION(mm)']) - ic_offset_y) * sad_y / (sad_y - ictoiso_y)
                            record_file_df['X_POSITION(mm)'], record_file_df['Y_POSITION(mm)'] = new_x_series, new_y_series
                            del new_x_series, new_y_series
                        
                            # charge to MU conversion using correction factor
                            record_file_df['MU'] = record_file_df['ACC_CHARGE(C)'] * correction_factor / charge_per_mu
                            record_file_df.drop(columns=['ACC_CHARGE(C)'], inplace=True)
                            record_file_df.reindex()  # make sure modified layer df is consistent with indexing
                            to_do_layers.append(record_file_df)
                
                    for tuning_file in tunings:  # do the same for all tuning files
                        if int(tuning_file.split('_')[2].split('_')[0]) == layer_id:
                            try:
                                tuning_file_df = pd.read_csv(tuning_file, delimiter=',', skiprows=10, skipfooter=11, engine='python')
                            except:
                                print('Read CSV error:', fraction_id, beam_id, record_file)
                                continue
                            try:
                                tuning_file_df['TIME'] = pd.to_datetime(tuning_file_df['TIME'])
                                tuning_file_df.index = tuning_file_df['TIME']
                                charge_col = pd.Series(tuning_file_df['DOSE_PRIM(C)'])
                                tuning_file_df = tuning_file_df.loc[:, :'Y_POSITION(mm)']
                                tuning_file_df['DOSE_PRIM(C)'] = charge_col
                                tuning_file_df.drop(columns=['TIME'], inplace=True)
                                tuning_file_df.drop(tuning_file_df[tuning_file_df['SUBMAP_NUMBER'] < 0].index, inplace=True)
                                tuning_file_df = tuning_file_df[tuning_file_df.groupby('SUBMAP_NUMBER')['SUBMAP_NUMBER'].transform('count') > 1]
                                tuning_file_df.drop(columns=['X_WIDTH(mm)', 'Y_WIDTH(mm)'], inplace=True)
                            except KeyError:
                                print(f'''\n  /!\ Key error occured while handling '{tuning_file}' (layer {str(layer_id).zfill(2)}), continuing..''')
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
                                tuning_file_df.loc[tuning_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['SPOT_ID']] = current_spot_id
                                spot_drill_time = len(tuning_file_df.loc[tuning_file_df['SUBMAP_NUMBER'] == current_spot_submap]) * 0.25
                                tuning_file_df.loc[tuning_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['DRILL_TIME(ms)']] = spot_drill_time
                                accumulated_charge = tuning_file_df.loc[tuning_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['DOSE_PRIM(C)']].abs().sum().mean()
                                tuning_file_df.loc[tuning_file_df['SUBMAP_NUMBER'] == current_spot_submap, ['ACC_CHARGE(C)']] = accumulated_charge

                                current_spot_submap = tuning_file_df.loc[tuning_file_df['SUBMAP_NUMBER'] > current_spot_submap]['SUBMAP_NUMBER'].min()
                                current_spot_id += 1
                            
                            tuning_file_df.drop(columns=['DOSE_PRIM(C)'], inplace=True)
                            tuning_file_df.drop_duplicates(subset=['SUBMAP_NUMBER'], keep='last', inplace=True)
                            tuning_file_df = tuning_file_df.loc[(tuning_file_df['X_POSITION(mm)'] != -10000.0) & (tuning_file_df['Y_POSITION(mm)'] != -10000.0)]
                            
                            for specif_file in tuning_specifs:
                                if int(specif_file.split('_')[2].split('_')[0]) == layer_id:
                                    with open(specif_file, 'r') as specif_file:
                                        lines = specif_file.readlines()
                                        ic_offsets = lines[3]
                                        ic_offset_x, ic_offset_y = float(ic_offsets.split(',')[3]), float(ic_offsets.split(',')[2])
                                        range_at_degrader = float(lines[1].split(',')[1])

                            nozzle_wet = np.polyval(nozzle_wet_poly, range_at_degrader)  # [mm]
                            range_at_iso = range_at_degrader - nozzle_wet
                            layer_energy = np.exp(np.polyval(iba_gtr2_poly, np.log(range_at_iso)))  # (MeV)
                            tuning_file_df['LAYER_ENERGY(MeV)'] = layer_energy

                            new_x_series = (pd.Series(tuning_file_df['Y_POSITION(mm)']) - ic_offset_x) * sad_x / (sad_x - ictoiso_x)  # coordinate system transform iba <-> raystation (x <-> y)
                            new_y_series = (pd.Series(tuning_file_df['X_POSITION(mm)']) - ic_offset_y) * sad_y / (sad_y - ictoiso_y)
                            tuning_file_df['X_POSITION(mm)'], tuning_file_df['Y_POSITION(mm)'] = new_x_series, new_y_series
                            del new_x_series, new_y_series            
                            
                            tuning_file_df['MU'] = tuning_file_df['ACC_CHARGE(C)'] * correction_factor / charge_per_mu
                            tuning_file_df.drop(columns=['ACC_CHARGE(C)'], inplace=True)

                            tuning_file_df.reindex()
                            to_do_tunings.append(tuning_file_df)
                    
                    for i in range(len(to_do_tunings)):  # in case of multiple layer parts: enable continuous spot indexing
                        if i > 0:  
                            to_do_tunings[i]['SPOT_ID'] += (to_do_tunings[i - 1]['SPOT_ID'].max() + 1)
                    for j in range(len(to_do_layers)):
                        if j > 0:
                            to_do_layers[j]['SPOT_ID'] += (to_do_layers[j - 1]['SPOT_ID'].max() + 1)
                    
                    if len(to_do_layers) > 0:  # can be zero, if only one spot in layer and omitted by high-weighted tuning
                        layer_df = pd.concat(to_do_layers)  # concatenate layers, assign additional columns
                        layer_df['LAYER_ID'] = layer_id
                        layer_df['TOTAL_LAYERS'] = num_layers
                        layer_df['BEAM_ID'] = beam_id
                        layer_df['GANTRY_ANGLE'] = gantry_angle
                        layer_df['FRACTION_ID'] = fraction_id
                        # layer_df['PATIENT_ID'] = self.patient_id
                        layer_df.drop(columns=['SUBMAP_NUMBER'], inplace=True)
                    else:
                        print(f'  /!\ No record found for layer-ID {layer_id} in beam {beam_id}, continuing..')
                        continue
                    
                    if len(to_do_tunings) > 0:
                        tuning_df = pd.concat(to_do_tunings)
                        tuning_df['LAYER_ID'] = layer_id
                        tuning_df['TOTAL_LAYERS'] = num_layers
                        tuning_df['BEAM_ID'] = beam_id
                        tuning_df['GANTRY_ANGLE'] = gantry_angle
                        tuning_df['FRACTION_ID'] = fraction_id
                        # tuning_df['PATIENT_ID'] = self.patient_id
                        tuning_df.drop(columns=['SUBMAP_NUMBER'], inplace=True)
                    else:
                        print(f'  /!\ No tunings found for layer-ID {layer_id} and beam {beam_id}, continuing..')
                        continue
                    
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

        # write out as .csv
        os.chdir(self.df_destination)
        print(f'''  ..Writing dataframe to '{self.df_destination}' as .CSV.. ''')
        while True:   
            try:
                self.patient_record_df.to_csv(f'patient_{self.patient_id}_records_data.csv')
                self.patient_tuning_df.to_csv(f'patient_{self.patient_id}_tunings_data.csv')
                break
            except PermissionError:
                input('  Permission denied, close target file and press ENTER.. ')
        self.patient_record_df = self.patient_record_df.astype(dtype={'BEAM_ID':str, 'FRACTION_ID':str})
        self.patient_tuning_df = self.patient_tuning_df.astype(dtype={'BEAM_ID':str, 'FRACTION_ID':str})
        self.record_df_name = f'patient_{self.patient_id}_records_data.csv'
        self.tuning_df_name = f'patient_{self.patient_id}_tuning_data.csv'
        print('Complete')
    

    def dicom_finder(self, fraction_id, beam_id):
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
                            
        if not found:
            print(f'  /!\ Auto-location failed for beam-ID {beam_id} in fraction-ID {fraction_id}, select plan DICOM..')
            pydicom.read(filedialog.askopenfilename(initialdir=os.path.join(self.logfile_dir, '..')))
        
        return beam_dcm


    def plot_beam_layers(self):     # For all layers in one beam, averaged over all fractions:
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

        dcm_beam = self.dicom_finder(fraction_id=scope_record_df['FRACTION_ID'].iloc[0], beam_id=beam_id)

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
            
            if beam_id == '03' and layer_id == 10:
                fig2, ax2 = plt.subplots(figsize=(6, 6))
                ax2.plot(plan_x_positions, plan_y_positions, marker='x', linestyle='-', label='Planned')
                ax2.plot(x_positions, y_positions, marker='o', markerfacecolor='None', linestyle='--', color='grey', label='Log-file original')
                ax2.plot(*zip(*xy_positions_sorted), marker='o', markerfacecolor='None', linestyle='-', color='black', label='Log-file sorted')
                ax2.plot(x_tunings, y_tunings, marker='o', markerfacecolor='None', linestyle='None', color='limegreen', label='Tuning spot(s)')
                # ax2.annotate(f'Beam {beam_id} | Layer #{str(layer_id + 1).zfill(2)} | $\Delta$ = {abs(len(plan_x_positions) - len(x_positions))}', xy=(1.0, 1.0), xycoords='axes points')
                ax2.set_xlabel('X [mm]')
                ax2.set_ylabel('Y [mm]')
                ax2.legend()
                fig2.tight_layout()
                fig2.savefig(f'{output_dir}/{self.patient_id}_{beam_id}_spotmap_KS.png', dpi=600)

            axs[layer_id].plot(plan_x_positions, plan_y_positions, marker='x', linestyle='-', markersize=2.0, markeredgewidth=0.2, linewidth=0.2, label='Planned')
            axs[layer_id].plot(*zip(*xy_positions_sorted), marker='o', markerfacecolor='None', linestyle='-', color='black', markersize=2.0, markeredgewidth=0.2, linewidth=0.2, label='Log-file sorted')
            axs[layer_id].plot(x_positions, y_positions, marker='o', markerfacecolor='None', linestyle='--', color='black', markersize=2.0, markeredgewidth=0.2, linewidth=0.2, label='Log-file original')
            axs[layer_id].plot(x_tunings, y_tunings, marker='o', markerfacecolor='None', linestyle='None', markersize=2.0, markeredgewidth=0.2, color='limegreen', label='Tuning spot(s)')
            axs[layer_id].annotate(f'Layer #{str(layer_id + 1).zfill(2)} | $\Delta$ = {abs(len(plan_x_positions) - len(x_positions))}', xy=(1.0, 1.0), xycoords='axes points', fontsize=8)
            if len(x_tunings) > 1:
                axs[layer_id].annotate(f'Layer #{str(layer_id + 1).zfill(2)} | $\Delta$ = {abs(len(plan_x_positions) - len(x_positions))} | $n_t$ > 1', xy=(1.0, 1.0), xycoords='axes points', fontsize=8) 
            axs[layer_id].legend(loc='upper right', fontsize=8)
        plt.suptitle(f'Spot Positions for Patient-ID {self.patient_id}, Beam-ID: {beam_id}', fontweight='bold', y=0.9)
        ax0.set_xlabel('X [mm]', fontweight='bold', labelpad=10)
        ax0.set_ylabel('Y [mm]', fontweight='bold', labelpad=10)
        plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        while True:
            try:
                fig.savefig(f'{output_dir}/{self.patient_id}_{beam_id}_spots.pdf')
                break
            except PermissionError:
                input('  Permission denied, close target file and press ENTER.. ')
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
        for df_file in os.listdir(self.df_destination):
            if df_file.__contains__('delta') and df_file.__contains__(str(self.patient_id)) and df_file.endswith('.csv'):
                re_init = input(f'''Found existing patient deltaframe '{df_file}', re-initialize [y/n]? ''')
                if re_init == 'y':
                    break
                else:
                    return None

        # initialize deltaframe
        self.patient_delta_df = pd.DataFrame(columns=self.patient_record_df.columns, index=self.patient_record_df.index)
        self.patient_delta_df.rename(columns={  'X_POSITION(mm)':'DELTA_X(mm)',
                                                'Y_POSITION(mm)':'DELTA_Y(mm)',
                                                'MU':'DELTA_MU',
                                                'LAYER_ENERGY(MeV)':'DELTA_E(MeV)'
                                             }, inplace=True)
        self.patient_delta_df['UNIQUE_INDEX'] = np.nan
        self.patient_delta_df.index = self.patient_delta_df['UNIQUE_INDEX']
        self.patient_delta_df = self.patient_delta_df.iloc[0:0]

        # find corresponding plan dicoms
        beam_list = [str(i) for i in self.patient_record_df['BEAM_ID'].drop_duplicates()]
        beam_dcm_dict = {id:[] for id in beam_list}

        print('\nTrying to auto-locate patient plan dicoms..')  # find corresponding plan dicoms
        for path, dirnames, filenames in os.walk(os.path.join(self.logfile_dir, '..')):
            for fname in filenames:
                if fname.__contains__('RP') and fname.endswith('.dcm') and not fname.__contains__('log'):
                    ds = pydicom.read_file(os.path.join(path, fname))
                    for i, dcm_beam in enumerate(ds.IonBeamSequence):
                        if float(dcm_beam.IonControlPointSequence[0].GantryAngle) in self.patient_record_df['GANTRY_ANGLE'].drop_duplicates().to_list():
                            if dcm_beam.BeamName in beam_list:
                                beam_dcm_dict[dcm_beam.BeamName].append(os.path.join(path, fname))
                            if dcm_beam.BeamDescription in beam_list:
                                beam_dcm_dict[dcm_beam.BeamDescription].append(os.path.join(path, fname))
        
        valid_beams = 0
        invalid_keys = []
        for key in beam_dcm_dict:
            if beam_dcm_dict[key] == []:
                print(f'''  /!\ Beam-ID '{key}' not found in patient dicoms''')
                invalid_keys.append(key)
            else:
                valid_beams += 1
        
        for key in invalid_keys:
            del beam_dcm_dict[key]
        
        print(f'{valid_beams}/{len(beam_list)} recorded beams found in patient dicoms')
    	
        unique_index = 0
        for f, fx_id in enumerate(self.patient_record_df['FRACTION_ID'].drop_duplicates()):
            for beam_id in beam_list:
                for i in range(len(beam_dcm_dict[beam_id])):
                    try:
                        ds = pydicom.read_file(beam_dcm_dict[beam_id][i])
                        for dcm_beam in ds.IonBeamSequence:
                            if (dcm_beam.BeamName == beam_id or dcm_beam.BeamDescription == beam_id) and float(dcm_beam.IonControlPointSequence[0].GantryAngle) == self.patient_record_df.loc[self.patient_record_df['BEAM_ID'] == beam_id]['GANTRY_ANGLE'].drop_duplicates().mean():
                                plan_beam = dcm_beam
                                break                        
                        
                        beam_df = self.patient_record_df.loc[(self.patient_record_df['BEAM_ID'] == beam_id) & (self.patient_record_df['FRACTION_ID'] == fx_id)]
                        num_layers = beam_df['TOTAL_LAYERS'].iloc[0]
                        for layer_id in beam_df['LAYER_ID'].drop_duplicates():
                            try:
                                plan_layer = plan_beam.IonControlPointSequence[layer_id * 2]
                            except IndexError:
                                print(f'''  /!\ Layer-ID mismatch, skipping layer #{layer_id + 1} in beam {beam_id}, fraction {fx_id}''')
                                continue

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

                            log_layer = beam_df.loc[(beam_df['LAYER_ID'] == layer_id) & (beam_df['FRACTION_ID'] == fx_id)]
                            log_xy = [tup for tup in zip(log_layer['X_POSITION(mm)'], log_layer['Y_POSITION(mm)'])]
                            log_mu, log_e = log_layer['MU'], log_layer['LAYER_ENERGY(MeV)']
                            log_xy_sorted = [(np.nan, np.nan) for _ in range(max(len(log_xy), len(plan_xy)))]
                            log_mu_sorted = [np.nan for _ in log_xy_sorted]
                            delta_x, delta_y, delta_mu = [], [], []
                            delta_e = plan_e - log_e.mean()

                            # match (x,y)-positions to plan, transform MU list equally
                            for i, log_spot in enumerate(log_xy):
                                shifts = [np.array(plan_spot) - np.array(log_spot) for plan_spot in plan_xy]
                                dists = [np.abs((shift).dot(shift)) for shift in shifts]
                                index = dists.index(min(dists))
                                dx, dy = shifts[index]
                                delta_x.append(dx), delta_y.append(dy)
                                delta_mu.append(plan_mu[index] - log_mu[i])
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
                                plt.plot(*zip(*plan_xy), marker='x', ls='-', color='tab:blue', label='plan')
                                plt.plot(*zip(*log_xy), marker='o', ls='--', color='tab:grey', label='log')
                                plt.plot(*zip(*log_xy_sorted), marker='o', ls='-', color='black', label='sorted')
                                plt.legend()
                                plt.draw()
                                raise ValueError(f'  /!\ Log beam {beam_id} does not match plan beam {plan_beam.BeamName}, retrying..')

                            # calculate deltas
                            # delta_mu = [mu_p - mu_l for mu_p, mu_l in zip(plan_mu, log_mu_sorted)]
                            delta_e = plan_e - log_e.mean()

                            # generate new dataframe 
                            fx_delta_df = pd.DataFrame(columns=self.patient_delta_df.columns)
                            try:
                                fx_delta_df['UNIQUE_INDEX'] = [log_mu_sorted.index(i) + unique_index for i in log_mu]
                            except ValueError:
                                print('  /!\ Sorting error: Lengths are (planned, log, sorted):', len(plan_mu), len(log_mu), len(log_mu_sorted))
                                continue

                            fx_delta_df.index = fx_delta_df['UNIQUE_INDEX']
                            fx_delta_df.sort_index(inplace=True)
                            fx_delta_df.drop(columns=['UNIQUE_INDEX'], inplace=True)

                            fx_delta_df['DELTA_X(mm)'] = delta_x
                            fx_delta_df['DELTA_Y(mm)'] = delta_y
                            fx_delta_df['DELTA_MU'] = delta_mu
                            fx_delta_df['DELTA_E(MeV)'] = delta_e
                            fx_delta_df['LAYER_ID'] = layer_id
                            fx_delta_df['SPOT_ID'] = [id for id in range(len(delta_x))]
                            fx_delta_df['FRACTION_ID'] = fx_id
                            fx_delta_df['BEAM_ID'] = beam_id
                            fx_delta_df['GANTRY_ANGLE'] = float(dcm_beam.IonControlPointSequence[0].GantryAngle)
                            fx_delta_df['TOTAL_LAYERS'] = int(beam_df['TOTAL_LAYERS'].mean())

                            if self.patient_delta_df.empty:
                                self.patient_delta_df = fx_delta_df
                            else:
                                self.patient_delta_df = pd.concat([self.patient_delta_df, fx_delta_df], sort=True)
                            
                            unique_index += len(log_xy_sorted)

                            if layer_id == (num_layers - 1):  # progress visualization
                                print('  ', '[' + (layer_id + 1) * '#' + (num_layers - layer_id - 1) * '-' + ']', end=f' Beam {beam_id} complete\n')
                            else:
                                print('  ', '[' + (layer_id + 1) * '#' + (num_layers - layer_id - 1) * '-' + ']', end=f' Layer {str(layer_id + 1).zfill(2)}/{str(num_layers).zfill(2)}\r')
                    
                    except:
                        continue

            print(f'  ..Fraction {f + 1}/{self.num_fractions} complete..')

        # write deltaframe
        os.chdir(self.df_destination)
        print(f'''  ..Writing deltaframe to '{self.df_destination}' as .CSV.. ''')
        while True:   
            try:
                self.patient_delta_df.to_csv(f'patient_{self.patient_id}_plan-log_delta.csv')
                break
            except PermissionError:
                input('  Permission denied, close target file and press ENTER.. ')
        plt.show()
        print('Complete')


    def delta_dependencies(self):
        for file in os.listdir(self.df_destination):
            if file.__contains__(str(self.patient_id)) and file.__contains__('delta') and file.endswith('.csv'):
                print(f'''Found patient deltaframe '{file}', reading in..''')
                self.patient_delta_df = pd.read_csv(os.path.join(self.df_destination, file), index_col='UNIQUE_INDEX', dtype={'BEAM_ID':str, 'FRACTION_ID':str})
                break            

        try:
            if self.patient_record_df.shape != self.patient_delta_df.shape:
                print(f'  /!\ Dataframe shape does not match deltaframe shape [{self.patient_record_df.shape} vs. {self.patient_delta_df.shape}], proceed with caution..')
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
            fig, axs = plt.subplots(2, 2, figsize=(9, 8))
            ax0 = fig.add_subplot(111, frameon=False)
            ax0.set_xticks([])
            ax0.set_yticks([])
            axs = axs.flatten()
            axs[0].hist(self.patient_delta_df['DELTA_X(mm)'], bins=80, alpha=0.7, edgecolor='black', label=f'''$\mu_x =$ {self.patient_delta_df['DELTA_X(mm)'].mean():.3f} mm\n$\sigma_x =$ {self.patient_delta_df['DELTA_X(mm)'].std():.3f} mm''')
            axs[0].axvline(self.patient_delta_df['DELTA_X(mm)'].mean(), ls='-', color='black', lw=0.5)
            axs[0].axvline(0.0, ls='--', color='black', lw=0.5)
            axs[0].set_xlabel('$\Delta x$ [mm]')
            axs[0].set_xlim(-2, 2)
            axs[1].hist(self.patient_delta_df['DELTA_Y(mm)'], bins=110, alpha=0.7, edgecolor='black', label=f'''$\mu_y =$ {self.patient_delta_df['DELTA_Y(mm)'].mean():.3f} mm\n$\sigma_y =$ {self.patient_delta_df['DELTA_Y(mm)'].std():.3f} mm''')
            axs[1].axvline(self.patient_delta_df['DELTA_Y(mm)'].mean(), ls='-', color='black', lw=0.5)
            axs[1].axvline(0.0, ls='--', color='black', lw=0.5)
            axs[1].set_xlabel('$\Delta y$ [mm]')
            axs[1].set_xlim(-2, 2)
            axs[2].hist(self.patient_delta_df['DELTA_MU'], bins=1000, color='tab:green', alpha=0.7, edgecolor='black', label=f'''$\mu_D =$ {self.patient_delta_df['DELTA_MU'].mean():.3f} MU\n$\sigma_D =$ {self.patient_delta_df['DELTA_MU'].std():.3f} MU''')
            axs[2].axvline(self.patient_delta_df['DELTA_MU'].mean(), ls='-', color='black', lw=0.5)
            axs[2].axvline(0.0, ls='--', color='black', lw=0.5)
            axs[2].set_xlabel('$\Delta D$ [MU]')
            axs[2].set_xlim(-0.005, 0.005)
            # axs[2].set_ylim(0, 30000)
            axs[3].hist(self.patient_delta_df['DELTA_E(MeV)'].drop_duplicates(), bins=15, color='tab:red', alpha=0.7, edgecolor='black', label=f'''$\mu_E =$ {self.patient_delta_df['DELTA_E(MeV)'].mean():.3f} MeV\n$\sigma_E =$ {self.patient_delta_df['DELTA_E(MeV)'].std():.3f} MeV''')
            axs[3].axvline(self.patient_delta_df['DELTA_E(MeV)'].mean(), ls='-', color='black', lw=0.5)
            axs[3].axvline(0.0, ls='--', color='black', lw=0.5)
            axs[3].set_xlabel('$\Delta E$ [MeV]')
            axs[3].set_xlim(-0.04, 0.04)
            for ax in axs:
                ax.grid(axis='y', zorder=-1)
                ax.set_axisbelow(True)
                # ax.set_yticklabels([])
                ax.legend()
            # ax0.set_title(f'Delta Histograms for Patient-ID {self.patient_id}', fontweight='bold')
            plt.tight_layout()
            plt.savefig(f'{output_dir}/{self.patient_id}_histograms.png')
            plt.show()            
        
        if choice == 2:  # gantry angle vs. delta(x,y)
            other_dfs = [pd.read_csv(os.path.join(self.df_destination, file), index_col='UNIQUE_INDEX', dtype={'BEAM_ID':str, 'FRACTION_ID':str}) for file in os.listdir(self.df_destination) if file.__contains__('delta') and file.endswith('.csv')]
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
        plt.savefig(f'{output_dir}/{self.patient_id}_beam_timings_KS.png', dpi=600)        
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
            beam_dcm_dict = {id:[] for id in beam_list}

            # auto-locate corresponding plan
            print('\nTrying to auto-locate patient plan dicoms..')
            dicoms = []
            for path, dirnames, filenames in os.walk(os.path.join(self.logfile_dir, '..')):
                for fname in filenames:
                    if fname.__contains__('RP') and fname.endswith('.dcm') and not fname.__contains__('log'):
                        ds = pydicom.read_file(os.path.join(path, fname))
                        for i, dcm_beam in enumerate(ds.IonBeamSequence):
                            if float(dcm_beam.IonControlPointSequence[0].GantryAngle) in self.patient_record_df['GANTRY_ANGLE'].drop_duplicates().to_list():
                                if dcm_beam.BeamName in beam_list or dcm_beam.BeamDescription in beam_list:
                                    dicoms.append(os.path.join(path, fname))
            
            for i, dci in enumerate(dicoms):
                for j, dcj in enumerate(dicoms):
                    if i != j and dci == dcj:
                        dicoms.remove(dcj)

            if len(dicoms) == 1:
                plan_dcm = dicoms[0]
                print(f'''  Found matching dicom file {dicoms[0]}''')
            elif len(dicoms) > 1:
                print(f'''  /!\ Multiple dicom files matching beam combination in fraction-ID {fx_id}\n      Will proceed with {dicoms[0]}''')
                plan_dcm = dicoms[0]
            elif len(dicoms) == 0:
                print(f'  /x\ Fx-ID {fx_id}: No matching dicom files detected, please select manually..')
                root = Tk()
                plan_dcm = filedialog.askopenfilename()
            
            dcm_path = os.path.dirname(plan_dcm)
            
            print('Manipulating dicom..')
            ds = pydicom.read_file(plan_dcm)
            for i, (beam_id, plan_beam) in enumerate(zip(beam_list, ds.IonBeamSequence)):
                beam_df = target_record.loc[target_record['BEAM_ID'] == beam_id]
                total_layers = beam_df['TOTAL_LAYERS'][0]
                cumulative_mu = 0
                for layer_id in target_record.loc[target_record['BEAM_ID'] == beam_id]['LAYER_ID'].drop_duplicates():
                    layer_df = pd.concat([target_record.loc[(target_record['BEAM_ID'] == beam_id) & (target_record['LAYER_ID'] == layer_id)], target_tuning.loc[(target_tuning['BEAM_ID'] == beam_id) & (target_tuning['LAYER_ID'] == layer_id)]])
                    layer_xy, layer_mu = [], []
                    layer_energy = layer_df['LAYER_ENERGY(MeV)'].drop_duplicates().mean()
                    for spot in range(len(layer_df['X_POSITION(mm)'].to_list())):
                        layer_xy.append(layer_df['X_POSITION(mm)'][spot])
                        layer_xy.append(layer_df['Y_POSITION(mm)'][spot])
                        layer_mu.append(layer_df['MU'][spot])
                    
                    n_spots = len(layer_mu)

                    plan_beam.IonControlPointSequence[layer_id * 2].NumberOfScanSpotPositions = plan_beam.IonControlPointSequence[layer_id * 2 + 1].NumberOfScanSpotPositions = n_spots
                    plan_beam.NumberOfControlPoints = total_layers * 2

                    if mode == 'all':
                        if layer_id == 0:
                            print(f'Will overwrite planned positions and metersets for beam-iD {beam_id}..')
                        plan_beam.IonControlPointSequence[layer_id * 2].ScanSpotPositionMap = plan_beam.IonControlPointSequence[layer_id * 2 + 1].ScanSpotPositionMap = layer_xy
                        plan_beam.IonControlPointSequence[layer_id * 2].CumulativeMetersetWeight = cumulative_mu
                        cumulative_mu += sum(layer_mu)
                        plan_beam.IonControlPointSequence[layer_id * 2 + 1].CumulativeMetersetWeight = cumulative_mu
                        plan_beam.IonControlPointSequence[layer_id * 2].ScanSpotMetersetWeights = layer_mu
                        plan_beam.IonControlPointSequence[layer_id * 2 + 1].ScanSpotMetersetWeights = [0.0 for _ in range(len(layer_mu))]
                        # plan_beam.IonControlPointSequence[layer_id * 2].NominalBeamEnergy = plan_beam.IonControlPointSequence[layer_id * 2 + 1].NominalBeamEnergy = layer_energy
                    
                    else:
                        if layer_id == 0:
                            print(f'Sorting spots for beam-ID {beam_id}..')
                        plan_spotmap = plan_beam.IonControlPointSequence[layer_id * 2].ScanSpotPositionMap
                        plan_x, plan_y = [], []
                        new_plan_x, new_plan_y, new_plan_mu = [], [], []
                        for i, spot in enumerate(plan_spotmap):
                            if i % 2 == 0:
                                plan_x.append(spot)
                            else:
                                plan_y.append(spot)

                        plan_xy = [tup for tup in zip(plan_x, plan_y)]
                        plan_mu = plan_beam.IonControlPointSequence[layer_id * 2].ScanSpotMetersetWeights
                        log_xy = [tup for tup in zip(layer_df['X_POSITION(mm)'].to_list(), layer_df['Y_POSITION(mm)'].to_list())]
                        for i, plan_spot in enumerate(plan_xy):
                            shifts = [np.array(plan_spot) - np.array(log_spot) for log_spot in log_xy]
                            dists = [np.abs((shift).dot(shift)) for shift in shifts]
                            index = dists.index(min(dists))
                            new_plan_x.append(plan_xy[index][0]), new_plan_y.append(plan_xy[index][1]), new_plan_mu.append(plan_mu[index])
                        
                        # plt.plot(*zip(*log_xy), label='log', linestyle='-', marker='o')
                        # plt.plot(new_plan_x, new_plan_y, label='plan', linestyle='-', marker='o')
                        # plt.scatter(plan_mu, new_plan_mu, label='MU')
                        # plt.xlabel('Plan MU')
                        # plt.ylabel('Log MU')
                        # plt.legend()
                        # plt.show()

                ds.FractionGroupSequence[0].ReferencedBeamSequence[i].BeamMeterset = cumulative_mu
                ds.FractionGroupSequence[0].NumberOfFractionsPlanned = 1
                plan_beam.FinalCumulativeMetersetWeight = cumulative_mu
                
                if len(plan_beam.IonControlPointSequence) > 2 * total_layers:
                    diff = len(plan_beam.IonControlPointSequence) - 2 * total_layers
                    print(f'  /!\ Deleting {diff} entries from beam {beam_id} control point sequence..')
                    plan_beam.IonControlPointSequence = plan_beam.IonControlPointSequence[: - diff]
                
                ds.IonBeamSequence[i] = plan_beam
            
            sop_uid = str(ds.SOPInstanceUID).split('.')
            sop_uid[-1] = '99999'
            new_sop_uid = '.'.join(sop_uid)
            ds.SOPInstanceUID = new_sop_uid
            if not ds.RTPlanName.__contains__('log'):
                ds.RTPlanName += f'_log_{mode}'
            ds.RTPlanLabel = ds.RTPlanName
            ds.ReferringPhysicianName = 'Wolter^Lukas'
            ds.ApprovalStatus = 'UNAPPROVED'

            print('Writing dicom..')
            destination = os.path.join(dcm_path, f'RP{ds.SOPInstanceUID}_fx_{fx_id}_log_{mode}.dcm')
            pydicom.write_file(destination, ds)
            print(f'Wrote log-based plan to RP{ds.SOPInstanceUID}_fx_{fx_id}_log_{mode}.dcm')
                

if __name__ == '__main__':
    # log = MachineLog()
    # log.prepare_dataframe()
    # log.plot_beam_layers()    
    # log.plot_spot_statistics()
    # log.prepare_deltaframe()
    # log.delta_dependencies()
    # log.plan_creator(fraction='last', mode='all')
    # log.beam_timings()
    pass

else:
    print('>> Module', __name__, 'loaded')

# %%
from logfile_extractor import MachineLog
log = MachineLog()