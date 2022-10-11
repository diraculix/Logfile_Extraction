import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse

try:
    df_destination = r'N:\fs4-HPRT\HPRT-Docs\Lukas\Logfile_Extraction\dataframes'  # TO BE CHANGED
    os.chdir(df_destination)
except:
    df_destination = r'C:\Users\lukas\Documents\OncoRay HPRT\Logfile_Extraction_mobile\dataframes'
    os.chdir(df_destination)


def plot_qa_beams(qa_df):
    beam_ids = qa_df['BEAM_ID'].drop_duplicates().to_list()
    n_beams = len(beam_ids)
    for i in range(int(np.sqrt(n_beams)) - 1, n_beams):
        if n_beams % i == 0:
            rows = i
            cols = int(n_beams / i)
            break

    # fig, axs = plt.subplots(4, 3, sharex=True, sharey=True, figsize=(20/4, 20/3), dpi=100)  # initiate matrix-like layer plot
    fig, axs = plt.subplots(3, 4, sharex=True, sharey=True, figsize=(30/3, 30/4), dpi=100)  # initiate matrix-like layer plot
    ax0 = fig.add_subplot(111, frameon=False)
    ax0.set_xticks([]), ax0.set_yticks([])
    fig.subplots_adjust(hspace=0.0, wspace=0.0)
    axs = axs.ravel()
    for beam_no, beam_id in enumerate(beam_ids):
        beam_df = qa_df.loc[qa_df['BEAM_ID'] == beam_id]
        axs[beam_no].set_axisbelow(True)
        axs[beam_no].grid(which='major', axis='both')
        axs[beam_no].scatter(beam_df['X_POSITION(mm)'], beam_df['Y_POSITION(mm)'], marker='o', c=beam_df['X_WIDTH(mm)'], s=beam_df['X_WIDTH(mm)'], vmin=7, vmax=19, cmap='inferno')
        lim = 200
        axs[beam_no].set_xlim(-lim, lim)
        axs[beam_no].set_ylim(-lim, lim)
        axs[beam_no].annotate(f'''{beam_df['FRACTION_ID'].iloc[0]} - {beam_df['LAYER_ENERGY(MeV)'].mean():.3f} MeV''', xy=(1.0, 1.0), xycoords='axes points')
        
    
    plt.show()
    plt.clf()

    # qa_dataframe.plot('LAYER_ENERGY(MeV)', 'X_WIDTH(mm)', kind='scatter')
    # qa_dataframe.plot('LAYER_ENERGY(MeV)', 'Y_WIDTH(mm)', kind='scatter')
    
    energies = qa_dataframe['LAYER_ENERGY(MeV)'].drop_duplicates().to_list()

    for df in [qa_dataframe.loc[qa_dataframe['LAYER_ENERGY(MeV)'] == energy] for energy in energies]:
        bp1 = plt.boxplot(df['X_WIDTH(mm)'], positions=[df['LAYER_ENERGY(MeV)'].iloc[0] - 6], widths=[10], patch_artist=True)
        bp2 = plt.boxplot(df['Y_WIDTH(mm)'], positions=[df['LAYER_ENERGY(MeV)'].iloc[0] + 6], widths=[10], patch_artist=True)
        plt.setp(bp1['boxes'], facecolor='tab:blue', alpha=0.5)
        plt.setp(bp1['medians'], color='tab:blue')
        plt.setp(bp2['boxes'], facecolor='tab:orange', alpha=0.5)
        plt.setp(bp2['medians'], color='tab:orange')
    
    plt.xticks(energies, [f'{energy:.3f}' for energy in energies])
    plt.legend([bp1["boxes"][0], bp2["boxes"][0]], ['X_WIDTH(mm)', 'Y_WIDTH(mm)'], loc='upper right')
    plt.xlabel('ENERGY [MeV]')
    plt.ylabel('SPOT WIDTH [mm]')
    # plt.grid(axis='y', zorder=-1)
    plt.show()

    return None
        

if __name__ == '__main__':
    qa_dataframe = pd.read_csv('QA_2017-2022_records_data.csv', dtype={'BEAM_ID':str, 'FRACTION_ID':str})
    print(qa_dataframe['X_WIDTH(mm)'].min(), qa_dataframe['X_WIDTH(mm)'].max())
    print(qa_dataframe['Y_WIDTH(mm)'].min(), qa_dataframe['Y_WIDTH(mm)'].max())
    plot_qa_beams(qa_dataframe)