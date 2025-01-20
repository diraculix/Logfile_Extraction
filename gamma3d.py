'''
Tested with Python 3.11.1 & pymedphys 0.36.1

'''

import pymedphys, pydicom, os, time
from pymedphys import gamma
from scipy.spatial.distance import cdist
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing


def gamma_3d(dpt, dta, path_to_ref, path_to_eval, cutoff=10, interp=10, max=1.1):
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
        quiet=False
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
    '''
    Voxel-wise dose pass rate [%]
    Not used in manuscript
    '''
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
            calculate_passrates(arg)


def parallelize(args):
    n_cores = multiprocessing.cpu_count() - 2
    chunks = [args[i::n_cores] for i in range(n_cores)]
    
    print('Queueing jobs on', n_cores, 'threads...')
    with multiprocessing.Pool(n_cores) as p:
        p.map(run_chunk, chunks)


def lambda_to_first_fx(df, fraction_id):
    '''
    Prepare a dataframe for Lambda calculation vs first Fx instead of vs. plan. Requested during MedPhys Review 01
    '''
    # Step 1: Filter the first fraction
    first_fraction_df = df[df['FRACTION_ID'] == df['FRACTION_ID'].iloc[0]]
    first_fraction = first_fraction_df[['X_POSITION(mm)', 'Y_POSITION(mm)']].values

    # Step 2: Get the (X, Y) positions for the given fraction
    fraction_df = df[df['FRACTION_ID'] == fraction_id]
    fraction = fraction_df[['X_POSITION(mm)', 'Y_POSITION(mm)']].values
    
    # Step 3: Calculate pairwise Euclidean distances between points of the two fractions
    distances = cdist(fraction, first_fraction, metric='euclidean')
    
    # Step 4: If the number of points differs, drop the farthest points from the larger set
    if len(fraction) > len(first_fraction):
        # More points in selected fraction, drop the farthest ones
        farthest_indices = np.argsort(np.min(distances, axis=1))[-(len(fraction) - len(first_fraction)):]
        fraction_df = fraction_df.drop(fraction_df.index[farthest_indices])
        fraction = fraction_df[['X_POSITION(mm)', 'Y_POSITION(mm)']].values
    elif len(first_fraction) > len(fraction):
        # More points in the first fraction, drop the farthest ones
        farthest_indices = np.argsort(np.min(distances, axis=0))[-(len(first_fraction) - len(fraction)):]
        first_fraction_df = first_fraction_df.drop(first_fraction_df.index[farthest_indices])
        first_fraction = first_fraction_df[['X_POSITION(mm)', 'Y_POSITION(mm)']].values
    
    # Step 5: Recalculate pairwise distances after adjusting sizes
    distances = cdist(fraction, first_fraction, metric='euclidean')
    
    # Step 6: Find the closest points in the adjusted sets
    sorted_indices = np.argmin(distances, axis=1)
    min_distances = distances[np.arange(len(fraction)), sorted_indices]
    
    # Step 7: Add the distances to the original dataframe (carefully update the right rows)
    df.loc[fraction_df.index, 'DIST_TO_FIRST_FX'] = min_distances

    return df


def calculate_passrates(id):
    root_dir = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\Logfiles\converted'
    plan_dir = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\02_cCTPatients\Logfiles\DeliveredPlans'

    # print(f'>>> START ID {id}')
    doses = os.path.join(plan_dir, id + '\Doses')
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
    for fx in fractions:
        df = lambda_to_first_fx(df, fx)
    
    # target_df = pd.DataFrame(index=fractions, columns=[f'Beam_{bid}_LOGPASS' for bid in beams] + ['TOTAL_LOGPASS'] + [f'Beam_{bid}_GAMMAPASS' for bid in beams] + ['TOTAL_GAMMAPASS'])
    gamma_crits = [(1, 1)] # [(3, 3), (2, 2), (1, 1)]  # dpt, dta
    log_crits = [1]  # mm distance to plan spot
    target_df = pd.DataFrame(columns=['FRACTION_ID'] + [f'LOGPASS_{crit}_mm' for crit in log_crits] + [f'GAMMAPASS_{crit[0]}_{crit[1]}' for crit in gamma_crits])
    df['DISTANCE'] = np.sqrt(df['DELTA_X(mm)'].copy() ** 2 + df['DELTA_Y(mm)'].copy() ** 2)

    # calculate physical spot pass rate
    for i, fx in enumerate(fractions):
        if fx != fractions[-1]:
            print(f'  > Calculating Lambda passrates ({i + 1}/{len(fractions)})..', end='\r')  
        else:
            print(f'  > Calculating Lambda passrates ({i + 1}/{len(fractions)})..', end='\n')   
        for crit in log_crits:
            fx_df = df.loc[df.FRACTION_ID == fx]
            total_log_pass = fx_df.loc[fx_df.DISTANCE <= crit].MU.sum() / fx_df.MU.sum() * 100
            log_pass_to_first_fx = fx_df.loc[fx_df.DIST_TO_FIRST_FX <= crit].MU.sum() / fx_df.MU.sum() * 100
            target_df.loc[i, 'FRACTION_ID'] = fx
            target_df.loc[i, f'LOGPASS_{crit}_mm'] = total_log_pass
            target_df.loc[i, f'LOGPASS_TO_FIRST_FX{crit}_mm'] = log_pass_to_first_fx
    
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
                for file in os.listdir(os.path.join(doses, 'Fx01')):
                    if file.startswith('RD') and file.endswith('.dcm'):  # accept only RS export dose dicoms (these are the only dose files in /Doses)
                        ds = pydicom.read_file(os.path.join(doses, 'Fx01', file))
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
                target_df.loc[i, f'GAMMAPASS_{dpt}_{dta}'] = plan_gamma_pass * 100
            else:
                if crit == gamma_crits[0]: print(f'    /!\ Missing evaluation dose for fraction {fx}')
                target_df.loc[i, f'GAMMAPASS_{dpt}_{dta}'] = np.nan

    written = False
    if target_df.LOGPASS_TO_FIRST_FX1_mm.min() < 99.:
        print(f'/!\ {id}: Minimum Lambda = {target_df.LOGPASS_TO_FIRST_FX1_mm.min() }')
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
    all_fx = 0
    for id in ponaqua_qualified:
        pat_fx = calculate_passrates(id)
        continue
        print(id, pat_fx)
        all_fx += pat_fx
    
    print('ALL', all_fx)
    #     if int(id) in [671075,1230180,1635796,1683480,1676596,280735,1367926]:
    #         ponaqua_qualified.remove(id) 
            
    # to_do = ["1669130"]        
    # parallelize(to_do)

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
