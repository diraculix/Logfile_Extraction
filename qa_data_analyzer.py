print('Read dataframe .. ', end='\r')

import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize
import seaborn as sns


try:
    df_destination = r'N:\fs4-HPRT\HPRT-Docs\Lukas\Logfile_Extraction\dataframes'
    output_dir = r'N:\fs4-HPRT\HPRT-Docs\Lukas\Logfile_Extraction\output'
    os.chdir(df_destination)
except:
    df_destination = r'/home/luke/Logfile_Extraction/dataframes'
    output_dir = r'/home/luke/Logfile_Extraction/output'
    os.chdir(df_destination)


def plot_qa_beams(qa_df):
    beam_ids = qa_df['BEAM_ID'].drop_duplicates().to_list()
    qa_df = qa_df.astype({'FRACTION_ID':str})

    sns.set_theme()
    fig = plt.figure(figsize=(11, 8))
    g = sns.scatterplot(data=qa_df, x='X_POSITION(mm)', y='Y_POSITION(mm)', hue='LAYER_ENERGY(MeV)', size='X_WIDTH(mm)', edgecolor='black')
    sns.move_legend(g, "upper left", bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.savefig(f'''{output_dir}/QA_scatter_{qa_df['LAYER_ENERGY(MeV)'].mean()}_{qa_df['GANTRY_ANGLE'].mean()}.png''', dpi=300)
    plt.show()
    plt.clf()

    g = sns.histplot(data=qa_df[['DELTA_X(mm)','DELTA_Y(mm)']], kde=True)
    plt.xlabel('$\Delta$ (log - plan) [mm]')
    leg = g.axes.get_legend()
    new_title = 'Mean difference to plan'
    leg.set_title(new_title)
    new_labels = [f'''$\Delta x =$ {np.round(qa_df['DELTA_X(mm)'].median(), 3)} $\pm$ {np.round(qa_df['DELTA_X(mm)'].std(), 3)}''', f'''$\Delta y =$ {np.round(qa_df['DELTA_Y(mm)'].median(), 3)} $\pm$ {np.round(qa_df['DELTA_Y(mm)'].std(), 3)}''']
    for t, l in zip(leg.texts, new_labels):
        t.set_text(l)   
    plt.tight_layout()
    plt.savefig(f'''{output_dir}/QA_hist_{qa_df['LAYER_ENERGY(MeV)'].mean()}_{qa_df['GANTRY_ANGLE'].mean()}.png''', dpi=300)
    plt.show()
    plt.clf()

    sns.color_palette('hls', n_colors=8)
    g = sns.pairplot(data=qa_df, vars=['DELTA_X(mm)', 'DELTA_Y(mm)', 'GANTRY_ANGLE', 'LAYER_ENERGY(MeV)'], hue='GANTRY_ANGLE', corner=True)
    plt.tight_layout()
    plt.show()
    plt.clf()

    # only_one_angle = qa_dataframe.loc[qa_dataframe['GANTRY_ANGLE'] == 0.0]
    # only_one_energy = qa_dataframe.loc[qa_dataframe['LAYER_ENERGY(MeV)'] == 100.].astype({'GANTRY_ANGLE':float})
    # g = sns.pairplot(data=only_one_angle[['DELTA_X(mm)', 'DELTA_Y(mm)', 'X_WIDTH(mm)', 'Y_WIDTH(mm)', 'LAYER_ENERGY(MeV)', 'GANTRY_ANGLE']], vars=['DELTA_X(mm)', 'DELTA_Y(mm)', 'X_WIDTH(mm)', 'Y_WIDTH(mm)'], hue='LAYER_ENERGY(MeV)', corner=True)
    # plt.suptitle('Gantry Angle 0.0°', fontweight='bold')
    # plt.savefig(f'{output_dir}/QA_pairplot_one_angle.png', dpi=600)
    # plt.clf()
    
    # g = sns.pairplot(data=only_one_energy[['DELTA_X(mm)', 'DELTA_Y(mm)', 'X_WIDTH(mm)', 'Y_WIDTH(mm)', 'LAYER_ENERGY(MeV)', 'GANTRY_ANGLE']], vars=['DELTA_X(mm)', 'DELTA_Y(mm)', 'X_WIDTH(mm)', 'Y_WIDTH(mm)'], hue='GANTRY_ANGLE', corner=True)
    # plt.suptitle('Beam Energy 100.0 MeV', fontweight='bold')
    # plt.savefig(f'{output_dir}/QA_pairplot_one_energy.png', dpi=600)
    # plt.clf()

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


def fit_sin(t, y):
    tt, yy = np.array(t), np.array(y)
    
    def sinfunc(t, A, w, p, c):  return A * np.sin(w*t + p) + c

    guess_freq = 1 / 360
    guess_amp = np.std(yy) * 2.**0.5
    guess_offset = np.mean(yy)
    guess = np.array([guess_amp, 2.*np.pi*guess_freq, 0., guess_offset])

    popt, pcov = optimize.curve_fit(sinfunc, tt, yy, p0=guess)
    A, w, p, c = popt
    f = w/(2.*np.pi)
    fitfunc = lambda t: A * np.sin(w*t + p) + c

    print('Amplitude:', A, '\nPeriod:', 1/f, '\nPhase:', p, '\nOffset:', c)
    return fitfunc


def correct_gtr(qa_df):
    angles = sorted(qa_df.GANTRY_ANGLE.drop_duplicates().to_list())
    x_diff_scatter = zip(qa_df.GANTRY_ANGLE, qa_df['DELTA_X(mm)'])
    y_diff_scatter = zip(qa_df.GANTRY_ANGLE, qa_df['DELTA_Y(mm)'])
    x_mean_diffs = [qa_df.loc[qa_df.GANTRY_ANGLE == angle]['DELTA_X(mm)'].mean() for angle in angles]
    y_mean_diffs = [qa_df.loc[qa_df.GANTRY_ANGLE == angle]['DELTA_Y(mm)'].mean() for angle in angles]

    output_df = pd.DataFrame(columns=['GANTRY_ANGLE', 'DELTA_X[mm]', 'DELTA_Y[mm]'])
    output_df['GANTRY_ANGLE'] = angles
    output_df['DELTA_X[mm]'] = x_mean_diffs
    output_df['DELTA_Y[mm]'] = y_mean_diffs
    output_df.to_csv(f'{output_dir}/QA_angular_dependence.csv')

    angles.append(360.)
    zero_x, zero_y = x_mean_diffs[0], y_mean_diffs[0]
    x_mean_diffs.append(zero_x), y_mean_diffs.append(zero_y)

    fig, axs = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)
    x_axis = np.linspace(0, 360, 1000)
    sine_fit_x = fit_sin(angles, x_mean_diffs)
    sine_fit_y = fit_sin(angles, y_mean_diffs)

    axs[0].scatter(*zip(*x_diff_scatter), alpha=0.2, label='$\Delta x$ to plan', zorder=-1)
    axs[0].scatter(angles, x_mean_diffs, edgecolors='black', c='white', label='mean shift', zorder=1)
    axs[0].plot(x_axis, sine_fit_x(x_axis), c='black', label='Sine fit', zorder=0)

    axs[1].scatter(*zip(*y_diff_scatter), alpha=0.2, label='$\Delta y$ to plan', zorder=-1)
    axs[1].scatter(angles, y_mean_diffs, edgecolors='black', c='white', label='mean shift', zorder=1)
    axs[1].plot(x_axis, sine_fit_y(x_axis), c='black', label='Sine fit', zorder=0)
    
    axs[0].set_ylabel('Positional error [mm]')
    axs[0].set_xlabel('Gantry angle [°]')
    axs[0].legend()
    axs[1].set_xlabel('Gantry angle [°]')
    axs[1].legend()

    plt.tight_layout()
    plt.savefig(f'{output_dir}/QA_gantry_sine_fit.png', dpi=300)
    plt.show()
      

if __name__ == '__main__':
    qa_dataframe = pd.read_csv('QA_2017-2022_records_data.csv', dtype={'BEAM_ID':str, 'FRACTION_ID':int})
    # qa_dataframe = qa_dataframe.loc[qa_dataframe['FRACTION_ID'] > 20220000]
    print(qa_dataframe.FRACTION_ID.drop_duplicates())
    print('Read dataframe .. DONE')
    plot_qa_beams(qa_dataframe)
    correct_gtr(qa_dataframe)
