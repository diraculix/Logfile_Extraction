print('Read dataframe .. ', end='\r')

import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize
import seaborn as sns
plt.rcParams['axes.grid'] = True
plt.rcParams['axes.grid.axis'] = 'y'


try:
    df_destination = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\AutoPatSpecQA\01_SpotShape\Logfiles_Spotshape_QA'
    output_dir = r'N:\fs4-HPRT\HPRT-Docs\Lukas\Logfile_Extraction\output'
    os.chdir(df_destination)
except:
    df_destination = r'/home/luke/Logfile_Extraction/dataframes'
    output_dir = r'/home/luke/Logfile_Extraction/output'
    os.chdir(df_destination)

def rotate(qa_df):
    new_x = qa_df['Y_POSITION(mm)']
    new_y = - qa_df['X_POSITION(mm)']
    new_dx = qa_df['DELTA_Y(mm)']
    new_dy = - qa_df['DELTA_X(mm)']
    qa_df['X_POSITION(mm)'] = new_x
    qa_df['Y_POSITION(mm)'] = new_y
    qa_df['DELTA_X(mm)'] = new_dx
    qa_df['DELTA_Y(mm)'] = new_dy

    return qa_df

def plot_qa_beams(qa_df):
    beam_ids = qa_df['BEAM_ID'].drop_duplicates().to_list()
    qa_df = qa_df.astype({'FRACTION_ID':str})

    # sns.set_theme()
    # fig = plt.figure(figsize=(11, 8))
    # g = sns.scatterplot(data=qa_df, x='X_POSITION(mm)', y='Y_POSITION(mm)', hue='LAYER_ENERGY(MeV)', size='X_WIDTH(mm)', edgecolor='black')
    # sns.move_legend(g, "upper left", bbox_to_anchor=(1, 1))
    # plt.tight_layout()
    # plt.savefig(f'''{output_dir}/QA_scatter_{qa_df['LAYER_ENERGY(MeV)'].mean()}_{qa_df['GANTRY_ANGLE'].mean()}.png''', dpi=300)
    # plt.show()
    # plt.clf()

    # g = sns.histplot(data=qa_df[['DELTA_X(mm)','DELTA_Y(mm)']], kde=True)
    # plt.xlabel('$\Delta$ (log - plan) [mm]')
    # leg = g.axes.get_legend()
    # new_title = 'Mean difference to plan'
    # leg.set_title(new_title)
    # new_labels = [f'''$\Delta x =$ {np.round(qa_df['DELTA_X(mm)'].median(), 3)} $\pm$ {np.round(qa_df['DELTA_X(mm)'].std(), 3)}''', f'''$\Delta y =$ {np.round(qa_df['DELTA_Y(mm)'].median(), 3)} $\pm$ {np.round(qa_df['DELTA_Y(mm)'].std(), 3)}''']
    # for t, l in zip(leg.texts, new_labels):
    #     t.set_text(l)   
    # plt.tight_layout()
    # plt.savefig(f'''{output_dir}/QA_hist_{qa_df['LAYER_ENERGY(MeV)'].mean()}_{qa_df['GANTRY_ANGLE'].mean()}.png''', dpi=300)
    # plt.show()
    # plt.clf()

    energies, angles = qa_df['LAYER_ENERGY(MeV)'].drop_duplicates(), qa_df['GANTRY_ANGLE'].drop_duplicates()
    fig, ax = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    inv_width = 2
    for energy in energies:
        count = len(qa_df.loc[qa_df['LAYER_ENERGY(MeV)'] == energy])
        ax[0].bar(energy, count, color='tab:blue', width=(max(energies) - min(energies)) / len(energies) / inv_width)
        ax[0].set_ylabel('Count')
        ax[0].set_xlabel('Energy [MeV]')
    for gantry in angles:
        count = len(qa_df.loc[qa_df['GANTRY_ANGLE'] == gantry])
        ax[1].bar(gantry, count, color='tab:orange', width=(max(angles) - min(angles)) / len(angles) / inv_width)
        ax[1].set_xlabel('Gantry angle [°]')
    
    plt.suptitle('Spot shape QA energy/angle frequency', fontweight='bold')
    plt.show()
    # sns.color_palette('hls', n_colors=8)
    # g = sns.pairplot(data=qa_df, vars=['DELTA_X(mm)', 'DELTA_Y(mm)', 'X_WIDTH(mm)', 'Y_WIDTH(mm)'], hue='GANTRY_ANGLE', corner=True)
    # plt.tight_layout()
    # plt.show()
    # plt.clf()

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
    # fig, axs = plt.subplots(5, 5, sharex=True, sharey=True, figsize=(8, 8), dpi=100)  # initiate matrix-like layer plot
    # ax0 = fig.add_subplot(111, frameon=False)
    # ax0.set_xticks([]), ax0.set_yticks([])
    # fig.subplots_adjust(hspace=0.0, wspace=0.0)
    # axs = axs.ravel()
    # for beam_no, beam_id in enumerate(qa_df['FRACTION_ID'].drop_duplicates()):
    #     beam_df = qa_df.loc[(qa_df['BEAM_ID'] == beam_id)]
    #     axs[beam_no].set_axisbelow(True)
    #     axs[beam_no].grid(which='major', axis='both')
    #     axs[beam_no].scatter(beam_df['X_POSITION(mm)'], beam_df['Y_POSITION(mm)'], marker='o', c=beam_df['X_WIDTH(mm)'], s=beam_df['X_WIDTH(mm)'], vmin=7, vmax=19, cmap='inferno')
    #     lim = 200
    #     axs[beam_no].set_xlim(-lim, lim)
    #     axs[beam_no].set_ylim(-lim, lim)
    #     axs[beam_no].annotate(f'''{beam_df['FRACTION_ID'].iloc[0]} - {beam_df['LAYER_ENERGY(MeV)'].mean():.3f} MeV''', xy=(1.0, 1.0), xycoords='axes points')
        
    
    # plt.show()
    # plt.clf()

    # qa_dataframe.plot('LAYER_ENERGY(MeV)', 'X_WIDTH(mm)', kind='scatter')
    # qa_dataframe.plot('LAYER_ENERGY(MeV)', 'Y_WIDTH(mm)', kind='scatter')
    
    
    # qa_angles = np.linspace(0., 360., 8, endpoint=False)        

    # for energy in qa_energies:
    #     for df in [qa_dataframe.loc[(qa_dataframe['GANTRY_ANGLE'] == gtr) & (abs(qa_dataframe['LAYER_ENERGY(MeV)'] - energy) < 0.1 )] for gtr in qa_angles]:
    #         print(df['GANTRY_ANGLE'].iloc[0])
    #         try:
    #             bp1 = plt.boxplot(df['X_WIDTH(mm)'], positions=[df['GANTRY_ANGLE'].iloc[0] - 6], widths=[10], patch_artist=True)
    #             bp2 = plt.boxplot(df['Y_WIDTH(mm)'], positions=[df['GANTRY_ANGLE'].iloc[0] + 6], widths=[10], patch_artist=True)
    #         except:
    #             continue
    #         plt.setp(bp1['boxes'], facecolor='tab:blue', alpha=0.5)
    #         plt.setp(bp1['medians'], color='tab:blue')
    #         plt.setp(bp2['boxes'], facecolor='tab:orange', alpha=0.5)
    #         plt.setp(bp2['medians'], color='tab:orange')

    #     plt.xticks(qa_angles, [f'{gtr:.1f}' for gtr in qa_angles])
    #     plt.legend([bp1["boxes"][0], bp2["boxes"][0]], ['X_WIDTH(mm)', 'Y_WIDTH(mm)'], loc='upper right')
    #     plt.xlabel('GANTRY_ANGLE[°]')
    #     plt.ylabel('SPOT WIDTH [mm]')
    #     # plt.grid(axis='y', zorder=-1)
    #     plt.title(str(energy) + ' MeV')
    #     plt.show()

    # return None


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


def angle_plot(qa_df):
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    import seaborn.categorical
    seaborn.categorical._Old_Violin = seaborn.categorical._ViolinPlotter

    class _My_ViolinPlotter(seaborn.categorical._Old_Violin):

        def __init__(self, *args, **kwargs):
            super(_My_ViolinPlotter, self).__init__(*args, **kwargs)
            self.gray='black'

    seaborn.categorical._ViolinPlotter = _My_ViolinPlotter

    # qa_df = qa_df.loc[qa_df.FRACTION_ID.astype(int) > 20210000]  # no big difference
    angles = sorted(qa_df.GANTRY_ANGLE.drop_duplicates().to_list())
    x_diff_scatter = zip(qa_df.GANTRY_ANGLE, qa_df['DELTA_X(mm)'])
    y_diff_scatter = zip(qa_df.GANTRY_ANGLE, qa_df['DELTA_Y(mm)'])
    x_median_diffs = [qa_df.loc[qa_df.GANTRY_ANGLE == angle]['DELTA_X(mm)'].median() for angle in angles]
    y_median_diffs = [qa_df.loc[qa_df.GANTRY_ANGLE == angle]['DELTA_Y(mm)'].median() for angle in angles]

    output_df = pd.DataFrame(columns=['GANTRY_ANGLE', 'DELTA_X[mm]', 'DELTA_Y[mm]'])
    output_df['GANTRY_ANGLE'] = angles
    output_df['DELTA_X[mm]'] = x_median_diffs
    output_df['DELTA_Y[mm]'] = y_median_diffs
    output_df.to_csv(f'{output_dir}/QA_angular_dependence.csv')

    angles.append(360.)
    zero_x, zero_y = x_median_diffs[0], y_median_diffs[0]
    x_median_diffs.append(zero_x), y_median_diffs.append(zero_y)

    fig, axs = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)
    x_axis = np.linspace(0, 315, 1000)
    sine_fit_x = fit_sin(angles, x_median_diffs)
    sine_fit_y = fit_sin(angles, y_median_diffs)

    sns.set(style='whitegrid', context='talk', font_scale=0.8)
    # axs[0].scatter(*zip(*x_diff_scatter), alpha=0.2, label='$\Delta x$ to plan', zorder=10)
    # axs[0].scatter(angles, x_mean_diffs, edgecolors='black', c='white', label='mean shift', zorder=30)
    axs[0].plot(x_axis, sine_fit_x(x_axis), c='black', label='Sine fit')
    sns.violinplot(data=qa_df, x=qa_df.GANTRY_ANGLE, y=qa_df['DELTA_X(mm)'], ax=axs[0], color=sns.color_palette()[0], lc='black', order=np.arange(360), width=20, scale='width', cut=0)

    # axs[1].scatter(*zip(*y_diff_scatter), alpha=0.2, label='$\Delta y$ to plan', zorder=10)
    # axs[1].scatter(angles, y_mean_diffs, edgecolors='black', c='white', label='mean shift', zorder=30)
    axs[1].plot(x_axis, sine_fit_y(x_axis), c='black', label='Sine fit')
    sns.violinplot(data=qa_df, x=qa_df.GANTRY_ANGLE, y=qa_df['DELTA_Y(mm)'], ax=axs[1], color=sns.color_palette()[1], lc='black', order=np.arange(360), width=20, scale='width', cut=0)
    
    legend1 = [Patch(facecolor=sns.color_palette()[0], edgecolor='None', label='$\Delta x$ distribution'),
              Line2D([0], [0], marker='o', linewidth=8, color='black', markerfacecolor='white', label='IQR & Median'),
              Line2D([0], [0], color='black', label='Sine fit')]
    legend2 = [Patch(facecolor=sns.color_palette()[1], edgecolor='None', label='$\Delta y$ distribution'),
              Line2D([0], [0], marker='o', linewidth=8, color='black', markerfacecolor='white', label='IQR & Median'),
              Line2D([0], [0], color='black', label='Sine fit')]

    for ax in axs:
        ax.set_xlabel('Gantry angle [°]')
        ax.set_xlim(-20, 335)
        ax.set_ylabel(None)
        # tick_every = 45
        # [tick.set_visible(False) for (i, tick) in enumerate(ax.get_xticklabels()) if i % tick_every != 0]
        ax.set_xticks([0, 45, 90, 135, 180, 225, 270, 315])
        ax.set_axisbelow(True)
        for collect in ax.collections:
                collect.set_edgecolor('None')
    
    axs[0].set_ylabel('Spot position error [mm]')
    axs[0].legend(handles=legend1, loc='lower right')
    axs[0].set_ylim(-3., 3.)
    axs[1].legend(handles=legend2, loc='lower right')

    plt.tight_layout()
    plt.savefig(f'{output_dir}/QA_gantry_sine_fit.png', dpi=300)
    plt.show()


def significance(qa_df):
    from scipy.stats import normaltest, ttest_ind, kstest
    zero_df = qa_df.loc[qa_df.GANTRY_ANGLE == 0.]
    angles = sorted(qa_df.GANTRY_ANGLE.drop_duplicates().to_list())
    for g in angles:
        df = qa_df.loc[qa_df.GANTRY_ANGLE == g]
        print(g)
        print('\tX:\t', ttest_ind(zero_df['DELTA_X(mm)'], df['DELTA_X(mm)']).pvalue)
        print('\tY:\t', ttest_ind(zero_df['DELTA_Y(mm)'], df['DELTA_Y(mm)']).pvalue)

    df90 = qa_df.loc[qa_df.GANTRY_ANGLE == 90.]
    df270 = qa_df.loc[qa_df.GANTRY_ANGLE == 270.]
    sns.histplot(df90['DELTA_X(mm)'])
    sns.histplot(df270['DELTA_X(mm)'])
    plt.show()


def histograms(qa_df):
    # statistical info
    print('Error\tMAX\tMEAN\tSTD')
    print(f"X:\t{qa_df['DELTA_X(mm)'].max().round(3)}\t{qa_df['DELTA_X(mm)'].mean().round(3)}\t{qa_df['DELTA_X(mm)'].std().round(3)}")
    print(f"Y:\t{qa_df['DELTA_Y(mm)'].max().round(3)}\t{qa_df['DELTA_Y(mm)'].mean().round(3)}\t{qa_df['DELTA_Y(mm)'].std().round(3)}")
    print(f'Delta > 2mm violations:')
    print(f"X:\t{len(qa_df.loc[abs(qa_df['DELTA_X(mm)']) > 2.])}")
    print(f"Y:\t{len(qa_df.loc[abs(qa_df['DELTA_Y(mm)']) > 2.])}")

    sns.set(style='ticks', context='talk')
    sns.histplot(data=qa_df['DELTA_X(mm)'], binwidth=0.05, label='$\Delta x$')
    sns.histplot(data=qa_df['DELTA_Y(mm)'], binwidth=0.05, label='$\Delta y$')
    plt.legend()
    plt.xlabel('Spot position error [mm]')
    plt.xlim(-3, 3)
    plt.grid()
    plt.tight_layout()
    plt.savefig(f'{output_dir}/QA_histograms.png', dpi=300)
    plt.yscale('log')
    plt.show()

def central_spot(qa_df):
    ctrl_df = qa_df.loc[(np.round(qa_df['X_POSITION(mm)'], 0) == 0.) & (np.round(qa_df['Y_POSITION(mm)'], 0) == 0.) & (qa_df.GANTRY_ANGLE.isin([0, 90, 180, 270]))]
    sns.set(style='whitegrid')
    for g in [0., 90., 180., 270.]:
        df = ctrl_df.loc[(ctrl_df.GANTRY_ANGLE == g) & (ctrl_df['LAYER_ENERGY(MeV)'] == 226.7)]
        plt.scatter(df['X_POSITION(mm)'].median(), df['Y_POSITION(mm)'].median(), label=g)
    
    plt.show()

def missing_beams(qa_df):
    qa_energies = [100., 120., 140., 165., 185., 205., 226.7]
    qa_angles = np.linspace(0., 360., 8, endpoint=False)
    fx_list = qa_df.FRACTION_ID.drop_duplicates()

    for fx in fx_list:
        fx_df = qa_df.loc[qa_df['FRACTION_ID'] == fx]
        missing_energies = []
        print(fx)
        for energy in qa_energies:
            energy_df = fx_df.loc[fx_df['LAYER_ENERGY(MeV)'] == energy]
            if energy_df.empty:
                missing_energies.append(energy)
                continue
           
            print('\t', energy) 
            angles, missing_angles = energy_df['GANTRY_ANGLE'].drop_duplicates().to_list(), []
            for gantry in qa_angles:
                if not gantry in angles:
                    missing_angles.append(gantry)

            print('\t\t', 'Missing angles:', missing_angles)
        
        print('\t', 'Missing energies:', missing_energies)  
                  

if __name__ == '__main__':
    qa_dataframe = pd.read_csv('QA_2017-2023_records_data.csv', dtype={'BEAM_ID':str, 'FRACTION_ID':int})
    # qa_dataframe = qa_dataframe.loc[qa_dataframe['FRACTION_ID'] > 20220000]
    # print(qa_dataframe.FRACTION_ID.drop_duplicates())
    print('Read dataframe .. DONE')
    # plot_qa_beams(qa_dataframe)
    qa_dataframe = rotate(qa_dataframe)
    angle_plot(qa_dataframe)
    # histograms(qa_dataframe)
    # significance(qa_dataframe)
    # central_spot(qa_dataframe)
    # missing_beams(qa_dataframe)
