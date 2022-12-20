import os
import numpy as np
import pandas as pd
from scipy import interpolate, optimize
from tkinter import filedialog
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams['axes.grid'] = True
plt.rcParams["axes.grid.axis"] = 'both'


beam_config_old = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\4D-PBS-LogFileBasedRecalc\Patient_dose_reconstruction\MOBILTest01_1588055\Logfiles\20191028\01-G190\Neuer Ordner\Neuer Ordner\beam_config.20191028_122536_576.csv'
beam_config_new = r'N:\fs4-HPRT\HPRT-Data\ONGOING_PROJECTS\4D-PBS-LogFileBasedRecalc\Patient_dose_reconstruction\MOBIL001_671075\Logfiles\20210901\01\NewFolder\24543\beam_config.20210901_180757_314.csv'


try:
    df_destination = r'N:\fs4-HPRT\HPRT-Docs\Lukas\Logfile_Extraction\dataframes'
    output_dir = r'N:\fs4-HPRT\HPRT-Docs\Lukas\Logfile_Extraction\output'
    os.chdir(df_destination)
except:
    df_destination = r'/home/luke/Logfile_Extraction/dataframes'
    output_dir = r'/home/luke/Logfile_Extraction/output'
    os.chdir(df_destination)


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

    print('\nAmplitude:', A, '\nPeriod:', 1/f, '\nPhase:', p, '\nOffset:', c)
    return fitfunc


with open(beam_config_old, 'r') as beam_config_old:  # draw machine parameters from *beam_config.csv
    iso_disp_x, iso_disp_y = [], []
    for no, line in enumerate(beam_config_old.readlines()):
        if line.__contains__('Isodisplacement_snout_extentions'):
            columns = [float(s) for s in line.split(';')[-1].split(',')]
        elif line.__contains__('Isodisplacement_gantry_angles'):
            rows = [float(g) for g in line.split(';')[-1].split(',')]
        elif line.__contains__('Isodisplacement_X_'):
            iso_disp_x.append([float(dx) for dx in line.split(';')[-1].split(',')])
        elif line.__contains__('Isodisplacement_Y_'):
            iso_disp_y.append([float(dy) for dy in line.split(';')[-1].split(',')])

    iso_disp_x_cubic = interpolate.interp2d(columns, rows, iso_disp_x, kind='cubic')
    iso_disp_y_cubic = interpolate.interp2d(columns, rows, iso_disp_y, kind='cubic')
    snout_axis = np.linspace(84, 573, 1000)
    gantry_axis = np.linspace(0, 360, 1000)

    qa_df = pd.read_csv(f'{output_dir}/QA_angular_dependence.csv')
    qa_angles = qa_df.GANTRY_ANGLE.to_list()
    qa_angles.append(360.)
    qa_median_x = qa_df['DELTA_X[mm]'].to_list()
    qa_median_x.append(qa_median_x[0])
    qa_median_y = qa_df['DELTA_Y[mm]'].to_list()
    qa_median_y.append(qa_median_y[0])
    sine_x = fit_sin(qa_angles, qa_median_x)
    sine_y = fit_sin(qa_angles, qa_median_y)

    fig, axs = plt.subplots(5, 2, figsize=(8, 8), sharex=True, sharey=True)
    ax0 = fig.add_subplot(111, frameon=False)
    ax0.set_xticks([]), ax0.set_yticks([])
    fig.subplots_adjust(hspace=0.15, wspace=0.05, top=0.9)
    snouts = [84., 200., 300., 400., 537.]
    gantries = [0,30,60,90,120,150,180,210,240,270,300,330,360]

    for i, s in enumerate(snouts):
        axs[i, 0].annotate(f'Snout {s} mm', xy=(5, 0.75),fontstyle='italic')

        if s != 300.:
            axs[i, 0].plot(gantries, iso_disp_x_cubic(s, gantries), 'v', linestyle='None', c='tab:grey', markerfacecolor='white')
            axs[i, 0].plot(gantry_axis, iso_disp_x_cubic(s, gantry_axis), c='tab:grey', linestyle=':')
            axs[i, 1].plot(gantries, iso_disp_y_cubic(s, gantries), 'v', linestyle='None', c='tab:grey', markerfacecolor='white', zorder=1)
            axs[i, 1].plot(gantry_axis, iso_disp_y_cubic(s, gantry_axis), c='tab:grey', linestyle=':', zorder=0)
        else:
            axs[i, 0].scatter(qa_angles, qa_median_x, edgecolors='black', c='white', zorder=3)
            axs[i, 0].plot(gantry_axis, sine_x(gantry_axis), color='black', zorder=2)
            axs[i, 0].plot(gantries, [0 for _ in gantries], 'v', linestyle='None', c='tab:grey', markerfacecolor='white')
            axs[i, 0].plot(gantry_axis, [0 for _ in gantry_axis], c='tab:grey', linestyle=':')
            axs[i, 1].plot(gantries, [0 for _ in gantries], 'v', linestyle='None', c='tab:grey', markerfacecolor='white' ,zorder=1)
            axs[i, 1].plot(gantry_axis, [0 for _ in gantry_axis], c='tab:grey', linestyle=':', zorder=0)
    
    with open(beam_config_new, 'r') as beam_config_new:
        iso_disp_x, iso_disp_y = [], []
        for no, line in enumerate(beam_config_new.readlines()):
            if line.__contains__('Isodisplacement_snout_extentions'):
                columns = [float(s) for s in line.split(';')[-1].split(',')]
            elif line.__contains__('Isodisplacement_gantry_angles'):
                rows = [float(g) for g in line.split(';')[-1].split(',')]
            elif line.__contains__('Isodisplacement_X_'):
                iso_disp_x.append([float(dx) for dx in line.split(';')[-1].split(',')])
            elif line.__contains__('Isodisplacement_Y_'):
                iso_disp_y.append([float(dy) for dy in line.split(';')[-1].split(',')])

        iso_disp_x_cubic = interpolate.interp2d(columns, rows, iso_disp_x, kind='cubic')
        iso_disp_y_cubic = interpolate.interp2d(columns, rows, iso_disp_y, kind='cubic')

        for i, s in enumerate(snouts):

            if s != 300.:
                axs[i, 0].plot(gantries, iso_disp_x_cubic(s, gantries), '^', linestyle='None', c='tab:grey', markerfacecolor='white')
                axs[i, 0].plot(gantry_axis, iso_disp_x_cubic(s, gantry_axis), c='tab:grey', linestyle='--')
                axs[i, 1].plot(gantries, iso_disp_y_cubic(s, gantries), '^', linestyle='None', c='tab:grey', markerfacecolor='white', zorder=1)
                axs[i, 1].plot(gantry_axis, iso_disp_y_cubic(s, gantry_axis), c='tab:grey', linestyle='--', zorder=0)
            else:
                axs[i, 1].scatter(qa_angles, qa_median_y, edgecolors='black', c='white', zorder=3)
                axs[i, 1].plot(gantry_axis, sine_y(gantry_axis), color='black', zorder=2)
                axs[i, 0].plot(gantries, [0 for _ in gantries], '^', linestyle='None', c='tab:grey', markerfacecolor='white')
                axs[i, 0].plot(gantry_axis, [0 for _ in gantry_axis], c='tab:grey', linestyle='--')
                axs[i, 1].plot(gantries, [0 for _ in gantries], '^', linestyle='None', c='tab:grey', markerfacecolor='white', zorder=1)
                axs[i, 1].plot(gantry_axis, [0 for _ in gantry_axis], c='tab:grey', linestyle='--', zorder=0)

        beam_config_new.close()
    
    axs[0, 0].set_title('$x$-coordinate')
    axs[0, 1].set_title('$y$-coordinate')
    ax0.set_xlabel('Gantry angle [°]', labelpad=25)
    ax0.set_ylabel('Positional shift [mm]',labelpad=40)

    legend = [
        Line2D([0], [0], c='tab:grey', marker='v', markerfacecolor='white', linestyle=':', label='Deprecated look-up'),
        Line2D([0], [0], c='tab:grey', marker='^', markerfacecolor='white', linestyle='--', label='Current look-up'),
        Line2D([0], [0], c='black', marker='o', markerfacecolor='white', label='Observed shift (sine fit)'),
        # Line2D([0], [0], c='black', label='Sine fit')
    ]

    ax0.legend(handles=legend, loc='upper center', bbox_to_anchor=(0.5, 1.07))
    plt.setp(axs, ylim=(-1, 1), xlim=(0, 360))
    # plt.suptitle('Look-up table for isocenter displacement', fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/lookup_vs_sine.png', dpi=300)

    beam_config_old.close()
