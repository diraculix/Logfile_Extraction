import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

try:
    df_destination = r'N:\fs4-HPRT\HPRT-Docs\Lukas\Logfile_Extraction\dataframes'  # TO BE CHANGED
    os.chdir(df_destination)
except:
    df_destination = r'/home/luke/Logfile_Extraction/dataframes'
    os.chdir(df_destination)


def plot_qa_beams(qa_df):
    beam_ids = qa_df['BEAM_ID'].drop_duplicates().to_list()
    sns.set()
    g = sns.scatterplot(data=qa_df, x='X_POSITION(mm)', y='Y_POSITION(mm)', hue='LAYER_ENERGY(MeV)', size='X_WIDTH(mm)')
    plt.show()

    g = sns.histplot(data=qa_df[['DELTA_X(mm)','DELTA_Y(mm)']], kde=True)
    plt.xlabel('$\Delta$ (log - plan) [mm]')
    leg = g.axes.get_legend()
    new_title = 'Median diff to plan'
    leg.set_title(new_title)
    new_labels = [f'''$\\tilde x =$ {qa_df['DELTA_X(mm)'].median().round(3)}''', f'''$\\tilde y =$ {qa_df['DELTA_Y(mm)'].median().round(3)}''']
    for t, l in zip(leg.texts, new_labels):
        t.set_text(l)   
    plt.show()

    only_one_angle = qa_dataframe.loc[qa_dataframe['GANTRY_ANGLE'] == 0.0]
    only_one_energy = qa_dataframe.loc[qa_dataframe['LAYER_ENERGY(MeV)'] == 226.7].astype({'GANTRY_ANGLE':float})
    g = sns.pairplot(data=only_one_angle[['DELTA_X(mm)', 'DELTA_Y(mm)', 'X_WIDTH(mm)', 'Y_WIDTH(mm)', 'LAYER_ENERGY(MeV)', 'GANTRY_ANGLE']], vars=['DELTA_X(mm)', 'DELTA_Y(mm)', 'X_WIDTH(mm)', 'Y_WIDTH(mm)'], hue='LAYER_ENERGY(MeV)', corner=True)
    plt.suptitle('Gantry Angle 0.0°', fontweight='bold')
    plt.show()
    g = sns.pairplot(data=only_one_energy[['DELTA_X(mm)', 'DELTA_Y(mm)', 'X_WIDTH(mm)', 'Y_WIDTH(mm)', 'LAYER_ENERGY(MeV)', 'GANTRY_ANGLE']], vars=['DELTA_X(mm)', 'DELTA_Y(mm)', 'X_WIDTH(mm)', 'Y_WIDTH(mm)'], hue='GANTRY_ANGLE', corner=True)
    plt.suptitle('Beam Energy 226.7 MeV', fontweight='bold')
    plt.show()

    return None

    # fig, axs = plt.subplots(4, 3, sharex=True, sharey=True, figsize=(20/4, 20/3), dpi=100)  # initiate matrix-like layer plot
    fig, axs = plt.subplots(5, 5, sharex=True, sharey=True, figsize=(8, 8), dpi=100)  # initiate matrix-like layer plot
    ax0 = fig.add_subplot(111, frameon=False)
    ax0.set_xticks([]), ax0.set_yticks([])
    fig.subplots_adjust(hspace=0.0, wspace=0.0)
    axs = axs.ravel()
    for beam_no, beam_id in enumerate(beam_ids):
        beam_df = qa_df.loc[(qa_df['BEAM_ID'] == beam_id) & (np.round(qa_df['LAYER_ENERGY'], 1) == qa_energies[0])]
        axs[beam_no].set_axisbelow(True)
        axs[beam_no].grid(which='major', axis='both')
        axs[beam_no].scatter(beam_df['X_POSITION(mm)'], beam_df['Y_POSITION(mm)'], marker='o', c=beam_df['X_WIDTH(mm)'], s=beam_df['X_WIDTH(mm)'], vmin=7, vmax=19, cmap='inferno')
        lim = 200
        axs[beam_no].set_xlim(-lim, lim)
        axs[beam_no].set_ylim(-lim, lim)
        axs[beam_no].annotate(f'''{beam_df['FRACTION_ID'].iloc[0]} - {beam_df['LAYER_ENERGY(MeV)'].mean():.3f} MeV''', xy=(1.0, 1.0), xycoords='axes points')
        
    
    plt.show()
    plt.clf()

    qa_dataframe.plot('LAYER_ENERGY(MeV)', 'X_WIDTH(mm)', kind='scatter')
    qa_dataframe.plot('LAYER_ENERGY(MeV)', 'Y_WIDTH(mm)', kind='scatter')
    
    
    qa_angles = np.linspace(0., 360., 8, endpoint=False)        

    for energy in qa_energies:
        for df in [qa_dataframe.loc[(qa_dataframe['GANTRY_ANGLE'] == gtr) & (abs(qa_dataframe['LAYER_ENERGY(MeV)'] - energy) < 0.1 )] for gtr in qa_angles]:
            print(df['GANTRY_ANGLE'].iloc[0])
            try:
                bp1 = plt.boxplot(df['X_WIDTH(mm)'], positions=[df['GANTRY_ANGLE'].iloc[0] - 6], widths=[10], patch_artist=True)
                bp2 = plt.boxplot(df['Y_WIDTH(mm)'], positions=[df['GANTRY_ANGLE'].iloc[0] + 6], widths=[10], patch_artist=True)
            except:
                continue
            plt.setp(bp1['boxes'], facecolor='tab:blue', alpha=0.5)
            plt.setp(bp1['medians'], color='tab:blue')
            plt.setp(bp2['boxes'], facecolor='tab:orange', alpha=0.5)
            plt.setp(bp2['medians'], color='tab:orange')

        plt.xticks(qa_angles, [f'{gtr:.1f}' for gtr in qa_angles])
        plt.legend([bp1["boxes"][0], bp2["boxes"][0]], ['X_WIDTH(mm)', 'Y_WIDTH(mm)'], loc='upper right')
        plt.xlabel('GANTRY_ANGLE[°]')
        plt.ylabel('SPOT WIDTH [mm]')
        # plt.grid(axis='y', zorder=-1)
        plt.title(str(energy) + ' MeV')
        plt.show()

    return None
        

def post_process(qa_df):
    qa_energies = [100., 140., 165., 185., 205., 226.7]
    qa_xy = [0., 30., 60., 120.]
    qa_xy += [- xy for xy in qa_xy]
    
    # qa_df['XY'] = [np.array(_).round(-1) for _ in zip(qa_df['X_POSITION(mm)'], qa_df['Y_POSITION(mm)'])]
    qa_df = qa_df[np.round(qa_df['X_POSITION(mm)'], -1).isin(qa_xy) & np.round(qa_df['Y_POSITION(mm)'], -1).isin(qa_xy)]
    qa_df['DELTA_X(mm)'] = qa_df['X_POSITION(mm)'] - np.round(qa_df['X_POSITION(mm)'], -1)
    qa_df['DELTA_Y(mm)'] = qa_df['Y_POSITION(mm)'] - np.round(qa_df['Y_POSITION(mm)'], -1)
    qa_df['LAYER_ENERGY(MeV)'] = np.round(qa_df['LAYER_ENERGY(MeV)'], 1)

    return qa_df
    
    

if __name__ == '__main__':
    qa_dataframe = pd.read_csv('QA_2017-2022_records_data.csv', dtype={'BEAM_ID':str, 'FRACTION_ID':int})
    # print(qa_dataframe['X_WIDTH(mm)'].min(), qa_dataframe['X_WIDTH(mm)'].max())
    # print(qa_dataframe['Y_WIDTH(mm)'].min(), qa_dataframe['Y_WIDTH(mm)'].max())
    qa_dataframe = post_process(qa_dataframe)
    plot_qa_beams(qa_dataframe)
    # get_delta(qa_dataframe)
