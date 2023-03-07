'''
Tested with Python 3.6.3 & pymedphys 0.36.1
'''

import pymedphys, pydicom, os, time
from pymedphys import gamma
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import multiprocessing


def gamma_3d(dpt, dta, path_to_ref, path_to_eval, cutoff=10, interp=10, max=1.1):
    # Load the two DICOM dose files to be compared
    try:
        ref = pymedphys.dicom.zyx_and_dose_from_dataset(pydicom.read_file(path_to_ref))
        eval = pymedphys.dicom.zyx_and_dose_from_dataset(pydicom.read_file(path_to_eval))
    except:
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
        quiet=True
    )

    # Filter the gamma index matrix to remove noise
    filtered_gamma = gamma_analysis[~np.isnan(gamma_analysis)]

    # Generate return values
    passed = np.sum(filtered_gamma <= 1) / len(filtered_gamma)
    max_gamma = np.max(filtered_gamma)
    mean_gamma = np.mean(filtered_gamma)
    std_gamma = np.std(filtered_gamma)

    return passed, max_gamma, mean_gamma, std_gamma


def voxel_wise(dpt, path_to_ref, path_to_eval, cutoff=10):
    # Load the two DICOM dose files to be compared
    ref = pymedphys.dicom.zyx_and_dose_from_dataset(pydicom.read_file(path_to_ref))
    eval = pymedphys.dicom.zyx_and_dose_from_dataset(pydicom.read_file(path_to_eval))

    # Get the dose grids and coordinate arrays from the DICOM files
    dose_ref_grid, coord_ref = ref[:2]
    dose_eval_grid, coord_eval = eval[:2]

    ref_max = np.nanmax(coord_ref)
    mask = (coord_ref < (cutoff / 100) * ref_max) & (~np.isnan(coord_ref))
    # coord_ref[mask] = np.nan
    percent_difference = np.abs(coord_ref - coord_eval) / ref_max * 100
    passed = np.count_nonzero((percent_difference < dpt) & (~np.isnan(percent_difference))) / np.count_nonzero((~np.isnan(percent_difference)))
    print(passed)


def run_chunk(chunk):
        for arg in chunk:
            full_eval(arg)


def parallelize(args):
    n_cores = multiprocessing.cpu_count() - 2
    chunks = [args[i::n_cores] for i in range(n_cores)]
    
    print('Queueing jobs on', n_cores, 'threads...')
    with multiprocessing.Pool(n_cores) as p:
        p.map(run_chunk, chunks)


def full_eval(id):
        root_dir = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\Logfiles\converted'
        plan_dir = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\Logfiles\DeliveredPlans'

        print(f'>>> START ID {id}')
        doses = os.path.join(plan_dir, id, 'Doses')
        log_dir = os.path.join(root_dir, '..', 'extracted')
        for file in os.listdir(log_dir):
            if file.__contains__(str(id)) and file.__contains__('record') and file.endswith('csv'):
                record_df = pd.read_csv(os.path.join(log_dir, file), index_col='TIME', dtype={'FRACTION_ID':str, 'BEAM_ID':str})
            elif file.__contains__(str(id)) and file.__contains__('delta') and file.endswith('csv'):
                delta_df = pd.read_csv(os.path.join(log_dir, file), index_col='UNIQUE_INDEX', dtype={'FRACTION_ID':str, 'BEAM_ID':str})
        
        if len(record_df) == len(delta_df):
            record_df.reset_index(inplace=True), delta_df.reset_index(inplace=True)
            df = pd.concat([record_df, delta_df], axis=1)
            df = df.loc[:,~df.columns.duplicated()].copy()
        else:
            print('   /x\ DF length mismatch!')
            return None
        
        # accept only more than 15 fx
        for beam in df.BEAM_ID.drop_duplicates():
            beam_df = df.loc[df.BEAM_ID == beam]
            if len(beam_df.FRACTION_ID.drop_duplicates()) < 15:
                df.drop(beam_df.index, inplace=True)

        # prepare output
        beams = df.BEAM_ID.drop_duplicates().to_list()
        fractions = df.FRACTION_ID.drop_duplicates().to_list()
        # target_df = pd.DataFrame(index=fractions, columns=[f'Beam_{bid}_LOGPASS' for bid in beams] + ['TOTAL_LOGPASS'] + [f'Beam_{bid}_GAMMAPASS' for bid in beams] + ['TOTAL_GAMMAPASS'])
        gamma_crits = [(3, 3), (2, 2), (1, 1)]  # dpt, dta
        log_crits = [2, 1.5, 1]  # mm distance to plan spot
        target_df = pd.DataFrame(columns=['FRACTION_ID'] + [f'LOGPASS_{crit}_mm' for crit in log_crits] + [f'GAMMAPASS_{crit[0]}_{crit[1]}' for crit in gamma_crits])
        df['DISTANCE'] = np.sqrt(df['DELTA_X(mm)'].copy() ** 2 + df['DELTA_Y(mm)'].copy() ** 2)

        # calculate physical spot pass rate
        for i, fx in enumerate(fractions):
            if fx != fractions[-1]:
                print(f'  > Calculating physical passrates ({i + 1}/{len(fractions)})..', end='\r')  
            else:
                print(f'  > Calculating physical passrates ({i + 1}/{len(fractions)})..', end='\n')   
            for crit in log_crits:
                fx_df = df.loc[df.FRACTION_ID == fx]
                total_log_pass = fx_df.loc[fx_df.DISTANCE <= crit].MU.sum() / fx_df.MU.sum()
                target_df.loc[i, 'FRACTION_ID'] = fx
                target_df.loc[i, f'LOGPASS_{crit}_mm'] = total_log_pass
        
        # calculate gamma dose pass rate
        for i, fx in enumerate(fractions):
            if fx != fractions[-1]:
                print(f'  > Calculating gamma-3D passrates ({i + 1}/{len(fractions)})..', end='\r') 
            else:
                print(f'  > Calculating gamma-3D passrates ({i + 1}/{len(fractions)})..', end='\n') 
                
            for crit in gamma_crits:
                for bid in beams:
                    has_ref_plan, has_ref_beam, has_eval_plan, has_eval_beam = False, False, False, False
                    # get reference doses
                    for file in os.listdir(doses):
                        if file.startswith('RD') and file.endswith('.dcm'):  # accept only RS export dose dicoms (these are the only dose files in /Doses)
                            ds = pydicom.read_file(os.path.join(doses, file))
                            if ds.DoseSummationType == 'PLAN':
                                ref_plan_dose = ds
                                has_ref_plan = True
                            elif ds.DoseSummationType == 'BEAM':
                                if ds.InstanceNumber == int(bid):
                                    ref_beam_dose = ds
                                    has_ref_beam = True

                        if has_ref_plan and has_ref_beam: 
                            break

                    # get evaluation doses
                    for file in os.listdir(os.path.join(doses, 'eval')):
                        if file.startswith('RD') and file.endswith('.dcm'):
                            ds = pydicom.read_file(os.path.join(doses, 'eval', file))
                            if ds.SeriesDescription.__contains__(str(fx)):
                                if ds.DoseSummationType == 'PLAN':
                                    eval_plan_dose = ds
                                    has_eval_plan = True
                                elif ds.DoseSummationType == 'BEAM':
                                    if ds.InstanceNumber == int(bid):
                                        eval_beam_dose = ds
                                        has_eval_beam = True
                        
                        if has_eval_plan and has_eval_beam:
                            break

                    # dpt, dta = crit
                    # beam_gamma_pass, _, _, _ = gamma_3d(dpt, dta, ref_beam_dose, eval_beam_dose)
                    # target_df.loc[fx, f'Beam_{bid}_GAMMAPASS'] = beam_gamma_pass
                
                dpt, dta = crit
                target_df.loc[i, 'FRACTION_ID'] = fx
                if has_eval_plan and has_ref_plan:
                    plan_gamma_pass, _, _, _ = gamma_3d(dpt, dta, ref_plan_dose, eval_plan_dose)
                    target_df.loc[i, f'GAMMAPASS_{dpt}_{dta}'] = plan_gamma_pass
                else:
                    if crit == gamma_crits[0]: print(f'    /!\ Missing evaluation dose for fraction {fx}')
                    target_df.loc[i, f'GAMMAPASS_{dpt}_{dta}'] = np.nan

                written = False
                while not written:
                    try: 
                        target_df.to_csv(os.path.join(doses, f'{id}_3Dgamma.csv'))               
                        break
                    except PermissionError:
                        input(f'  /!\ Permission denied for target CSV, close file..')
                    except FileNotFoundError:
                        print(f'  /!\ Target location {doses} not found')
                        print('   ... Waiting for drive mapping ...')
                        time.sleep(5)
            

if __name__ == '__main__':
    # root_dir = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\Logfiles\DeliveredPlans\1663630\Doses'
    # eval_dir = os.path.join(root_dir, 'eval')
    # criteria = [(2, 2), (1, 1), (1, 0.5), (0.5, 1), (0.5, 0.5)]
    ponaqua_qualified = [id.strip('\n') for id in open(r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\qualified_IDs.txt', 'r').readlines()]
    for id in ponaqua_qualified:
        if int(id) in [671075,1230180,1635796,1683480,1676596,280735,1367926]:
            ponaqua_qualified.remove(id) 
    to_do = ["1632783", "1180747"]        
    parallelize(to_do)

    # ref = os.path.join(root_dir, 'RD_R1_HNO.dcm')
    # evals = [os.path.join(eval_dir, dcm) for dcm in os.listdir(eval_dir) if dcm.__contains__('RD') and dcm.endswith('.dcm')]
    # data = ['Fx\tPass\tMax\tMean\tStd\n']
    # for fx, eval in enumerate(evals):
    #     print(f'>> STARTING 3D gamma evaluation for  {fx+1}/{len(evals)} with {dpt}%/{dta}mm..')
    #     voxel_wise(1, ref, eval)
    #     passed, max, mean, std = gamma_3d(dpt, 1000, ref, eval)
    #     print(passed)
    #     # data.append(f'{fx+1}\t{np.round(passed, 5)}\t{np.round(max, 5)}\t{np.round(mean, 5)}\t{np.round(std, 5)}\n')
    
    # with open(os.path.join(root_dir, f'gamma3D_{dpt}p_{dta}mm.txt'), 'w+') as file:
    #     file.writelines(data)
    #     file.close()
