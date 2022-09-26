def beam_timings(self):
    fig, axs = plt.subplots(len(self.beam_list), 1, sharex=True, figsize=(10, 10))
    fig.subplots_adjust(hspace=0.0, wspace=0.0)
    ax0 = fig.add_subplot(111, frameon=False)
    ax0.set_xticks([]), ax0.set_yticks([])
    ax0.set_ylabel('Time [ms]', labelpad=20)
    ax0.set_xlabel('Fraction-ID')
    ax0.set_title(f'Beam timings for patient-ID {self.patient_id}', fontweight='bold')
    axs.flatten()
    x_axis = self.fraction_list
    
    for fx_id in x_axis:  # remember type(fx_id) = <str>
        for i, beam_id in enumerate(self.beam_list):
            beam_df = self.patient_record_df.loc[self.patient_record_df['BEAM_ID'] == beam_id]
            total_drill_time = beam_df['DRILL_TIME(ms)'].sum()
            total_layer_time, total_energy_switch = 0, 0
            layer_dfs = [beam_df.loc[beam_df['LAYER_ID'] = lid] for lid in beam_df['LAYER_ID'].drop_duplicates()]
            for layer_id, layer_df in enumerate(layer_dfs):
                start_this = layer_df.first_valid_index()
                end_this = layer_df.last_valid_index()
                layer_time = end_this - start_this
                total_layer_time += layer_time.milliseconds
                
                if layer_id > 0:
                    end_previous = layer_dfs[layer_id - 1].last_valid_index()
                    energy_switch = start_this - end_previous
                    total_energy_switch += energy_switch.milliseconds
            
            total_spot_switch = total_layer_time - total_drill_time
            
            axs[i].bar(fx_id, total_drill_time, label='Drill')
            axs[i].bar(fx_id, total_spot_switch, bottom=total_drill_time, label='Spot switch')
            axs[i].bar(fx_id, total_energy_switch, bottom=total_spot_switch, label='Energy switch')
    
    for ax, beam_id in zip(axs, self.beam_list):
        ax.set_ylabel(f'Beam {beam_id}')
        ax.legend()
    
    plt.show()
            
                
                
                    