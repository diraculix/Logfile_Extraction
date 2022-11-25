import os, time
import pandas as pd
from influxdb_client import InfluxDBClient
from influxdb_client.client.write_api import ASYNCHRONOUS

# db access config
token = "vHLCAHK5mURBcgOJane0P6GsOyBkaDg1H9cHTGwgEwWlRwS7UuFitjIewPhf_o9AZVwinuzHyDD6c9bL6HN2Ug=="
org = "OncoRay"
bucket = "GTR2"
url = "http://g40invivodossrv:8086"

# dataframe dump location
df_dump = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\InfluxDB\_push'


def stage_patients() -> None:
    global df_dump
    candidates = []
    for df in os.listdir(df_dump):
        pat_id = df.split('_')[1]
        if df.__contains__('_pushed') or df.__contains__('_staged'):
            continue
        
        if not pat_id in candidates and df.endswith('.csv'):
            candidates.append(pat_id)    
    
    return candidates


def push_to_influx(patients) -> None:
    global token, org, bucket, url, df_dump
    for pat_id in patients:
        count = 0  # complete dataset if count==3 (record+tuning+delta)
        for df in os.listdir(df_dump):
            if df.__contains__(pat_id):
                os.rename(os.path.join(df_dump, df), os.path.join(df_dump, f'''{df.replace('.', '_staged.')}'''))
                pass
        
        for df in os.listdir(df_dump):
            if df.__contains__(pat_id) and df.__contains__('_staged'):
                if df.__contains__('record'):
                    record = os.path.join(df_dump, df)
                elif df.__contains__('tuning'):
                    tuning = os.path.join(df_dump, df)
                elif df.__contains__('delta'):
                    delta = os.path.join(df_dump, df)

                count += 1
        
        if count == 3:
            print(f'[{time.ctime()}] > Influx: Push request is valid for patient-ID {pat_id}')
        else:
            print(f'[{time.ctime()}] > Influx: Invalid push request for patient-ID {pat_id}, dataframe(s) missing')
            return None

        # dataframe column renaming
        log_col_dict = {
            'X_POSITION(mm)':'LOG_X_POS_TPS(mm)', 'Y_POSITION(mm)':'LOG_Y_POS_TPS(mm)',
            'X_POS_IC23(mm)':'LOG_X_POS_IC23(mm)', 'Y_POS_IC23(mm)':'LOG_Y_POS_IC23(mm)',
            'X_POSITION_CORR(mm)':'LOG_X_POS_CORR(mm)', 'Y_POSITION_CORR(mm)':'LOG_Y_POS_CORR(mm)',
            'X_WIDTH(mm)':'LOG_X_WID_TPS(mm)', 'Y_WIDTH(mm)':'LOG_Y_WID_TPS(mm)',
            'X_WID_IC23(mm)':'LOG_X_WID_IC23(mm)', 'Y_WID_IC23(mm)':'LOG_Y_WID_IC23(mm)',
            'LAYER_ENERGY(MeV)':'LOG_ENERGY(MeV)', 'MU':'LOG_DOSE(MU)',
            'DRILL_TIME(ms)':'LOG_DRILL(ms)', 'CHARGE(C)':'LOG_CHARGE(C)'            
        }
        tune_col_dict = {
            'X_POSITION(mm)':'TUNE_X_POS_TPS(mm)', 'Y_POSITION(mm)':'TUNE_Y_POS_TPS(mm)',
            'X_POS_IC23(mm)':'TUNE_X_POS_IC23(mm)', 'Y_POS_IC23(mm)':'TUNE_Y_POS_IC23(mm)',
            'X_POSITION_CORR(mm)':'TUNE_X_POS_CORR(mm)', 'Y_POSITION_CORR(mm)':'TUNE_Y_POS_CORR(mm)',
            'X_WIDTH(mm)':'TUNE_X_WID_TPS(mm)', 'Y_WIDTH(mm)':'TUNE_Y_WID_TPS(mm)',
            'X_WID_IC23(mm)':'TUNE_X_WID_IC23(mm)', 'Y_WID_IC23(mm)':'TUNE_Y_WID_IC23(mm)',
            'LAYER_ENERGY(MeV)':'TUNE_ENERGY(MeV)', 'MU':'TUNE_DOSE(MU)',
            'DRILL_TIME(ms)':'TUNE_DRILL(ms)', 'CHARGE(C)':'TUNE_CHARGE(C)'            
        }
        delta_col_dict = {
            'DELTA_X(mm)':'DELTA_X_POS_TPS(mm)', 'DELTA_Y(mm)':'DELTA_Y_POS_TPS(mm)',
            'DELTA_E(MeV)':'DELTA_ENERGY(MeV)', 'DELTA_MU':'DELTA_DOSE(MU)'            
        }

        # read dataframe csv, rename cols
        print(f'[{time.ctime()}] > Influx:   Staging dataframes for pat.-ID {pat_id}')
        record_df = pd.read_csv(record, index_col='TIME', dtype={'FRACTION_ID':str, 'BEAM_ID':str}, parse_dates=True, dayfirst=True)
        tuning_df = pd.read_csv(tuning, index_col='TIME', dtype={'FRACTION_ID':str, 'BEAM_ID':str}, parse_dates=True, dayfirst=True)
        record_df.index = pd.to_datetime(record_df.index, format='%Y-%m-%dT%H:%M:%SZ')
        tuning_df.index = pd.to_datetime(tuning_df.index, format='%Y-%m-%dT%H:%M:%SZ')
        delta_df = pd.read_csv(delta, index_col='UNIQUE_INDEX', dtype={'FRACTION_ID':str, 'BEAM_ID':str})
        delta_df.index = record_df.index
        record_df.rename(log_col_dict, axis=1, inplace=True)
        tuning_df.rename(tune_col_dict, axis=1, inplace=True)
        delta_df.rename(delta_col_dict, axis=1, inplace=True)

        # concat record and delta dataframes
        joint_df = pd.concat([record_df, delta_df], axis=1)
        joint_df = joint_df.loc[:,~joint_df.columns.duplicated()]

        with open(f'{df_dump}/tuning_df_cols.txt', 'w+') as file:
            for col in tuning_df.columns.to_list():
                file.write(f'\'{col}\','+'\n')
            
            file.close()

        # change column sequence to match influx tag hierarchy
        joint_df = joint_df[[
            'PATIENT_ID',
            'FRACTION_ID',
            'BEAM_ID',
            'GANTRY_ANGLE',
            'PRESSURE(hPa)',
            'TEMPERATURE(K)',
            'TOTAL_LAYERS',
            'LAYER_ID',
            'SPOT_ID',
            'LOG_X_POS_IC23(mm)',
            'LOG_X_POS_TPS(mm)',
            'LOG_X_POS_CORR(mm)',
            'LOG_Y_POS_IC23(mm)',
            'LOG_Y_POS_TPS(mm)',
            'LOG_Y_POS_CORR(mm)',
            'LOG_X_WID_IC23(mm)',
            'LOG_X_WID_TPS(mm)',
            'LOG_Y_WID_TPS(mm)',
            'LOG_Y_WID_IC23(mm)',
            'LOG_ENERGY(MeV)',
            'LOG_CHARGE(C)',
            'LOG_DOSE(MU)',
            'LOG_DRILL(ms)',
            'DELTA_X_POS_TPS(mm)',
            'DELTA_Y_POS_TPS(mm)',
            'DELTA_ENERGY(MeV)',
            'DELTA_DOSE(MU)'
        ]]
        tuning_df = tuning_df[[
            'PATIENT_ID',
            'FRACTION_ID',
            'BEAM_ID',
            'GANTRY_ANGLE',
            # 'PRESSURE(hPa)',
            # 'TEMPERATURE(K)',
            'LAYER_ID',
            'SPOT_ID',
            'TUNE_X_POS_IC23(mm)',
            'TUNE_X_POS_TPS(mm)',
            'TUNE_Y_POS_IC23(mm)',
            'TUNE_Y_POS_TPS(mm)',
            'TUNE_X_WID_IC23(mm)',
            'TUNE_X_WID_TPS(mm)',
            'TUNE_Y_WID_TPS(mm)',
            'TUNE_Y_WID_IC23(mm)',
            'TUNE_ENERGY(MeV)',
            'TUNE_CHARGE(C)',
            'TUNE_DOSE(MU)',
            'TUNE_DRILL(ms)'
        ]]

        # push dataframe to db
        with InfluxDBClient(url=url, token=token, org=org) as client:
            print(f'[{time.ctime()}] > Influx:   Connected to {bucket}@{org}')
            write_client = client.write_api(write_options=ASYNCHRONOUS)
            print(f'[{time.ctime()}] > Influx:   Writing pat.-ID {pat_id} to database')
            write_client.write(
                bucket, 
                org, 
                record=joint_df, 
                data_frame_measurement_name=f'{joint_df.PATIENT_ID.iloc[0]}_testing',
                data_frame_tag_columns=['PATIENT_ID', 'FRACTION_ID', 'BEAM_ID', 'GANTRY_ANGLE', 'PRESSURE(hPa)', 'TEMPERATURE(K)', 'TOTAL_LAYERS', 'LAYER_ID', 'SPOT_ID']
            )
            write_client.write(
                bucket,
                org,
                record=tuning_df,
                data_frame_measurement_name=f'{tuning_df.PATIENT_ID.iloc[0]}_testing',
                data_frame_tag_columns=['PATIENT_ID', 'FRACTION_ID', 'BEAM_ID', 'GANTRY_ANGLE', 'PRESSURE(hPa)', 'TEMPERATURE(K)', 'LAYER_ID', 'SPOT_ID']
            )

        print(f'[{time.ctime()}] > Influx:   Push successful for pat.-ID {pat_id}')
        for df in [record, tuning, delta]:
                os.rename(os.path.join(df_dump, df), os.path.join(df_dump, f'''{df.replace('_staged.', '_pushed.')}'''))


if __name__=='__main__':
    patients = stage_patients()
    push_to_influx(patients)
