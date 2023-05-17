import tkinter as tk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import profiles as pr
import numpy as np


def get_index(x, y):
    return (np.abs(x - y)).argmin()


class QProfileGUI:
    def __init__(self, master):
        self.master = master
        master.title("Q Profile Plotter")

        self.equilibria = []
        self.transp_runs = {}

        self.eq_pulse_entries = []
        self.transp_pulse_entries = []
        self.runid_entries = []
        self.dda_entries = []
        self.uid_entries = []
        self.seq_entries = []

        self.time_entry = None
        self.plot_type_var = tk.StringVar()

        self.create_equilibrium_inputs()
        self.create_transp_inputs()
        self.create_plot_type_input()
        self.create_time_input()
        self.create_plot_button()

    def create_equilibrium_inputs(self):
        """Create widgets for equilibrium inputs."""
        eq_frame = tk.Frame(self.master)
        eq_frame.pack()

        eq_label = tk.Label(eq_frame, text="Equilibria")
        eq_label.pack()

        num_eq_label = tk.Label(eq_frame, text="Number of Equilibria:")
        num_eq_label.pack(side=tk.LEFT)

        self.num_eq_entry = tk.Entry(eq_frame)
        self.num_eq_entry.pack(side=tk.LEFT)
        self.num_eq_entry.insert(0, "0")

        num_eq_button = tk.Button(eq_frame, text="Set", command=self.update_equilibrium_inputs)
        num_eq_button.pack(side=tk.LEFT)

    def create_time_input(self):
        """Create widgets for time input."""
        time_frame = tk.Frame(self.master)
        time_frame.pack()

        time_label = tk.Label(time_frame, text="Time:")
        time_label.pack(side=tk.LEFT)

        self.time_entry = tk.Entry(time_frame)
        self.time_entry.pack(side=tk.LEFT)

        self.plot_type_var.trace('w', lambda *args: self.update_time_label(time_label))

    def update_time_label(self, time_label):
        plot_type = self.plot_type_var.get()

        if plot_type == 'Profile':
            time_label.config(text='Times:')
        elif plot_type == 'Time Trace':
            time_label.config(text='Positions:')

    def create_plot_type_input(self):
        """Create widgets for plot type input."""
        plot_type_frame = tk.Frame(self.master)
        plot_type_frame.pack()

        plot_type_label = tk.Label(plot_type_frame, text="Plot Type:")
        plot_type_label.pack(side=tk.LEFT)

        plot_type_options = ['Profile', 'Time Trace']

        for option in plot_type_options:
            plot_type_radio = tk.Radiobutton(plot_type_frame, text=option, variable=self.plot_type_var, value=option)
            plot_type_radio.pack(side=tk.LEFT)

        self.plot_type_var.set('Profile')  # Set the default plot type as 'Profile'

    def create_plot_button(self):
        """Create the plot button."""
        plot_button = tk.Button(self.master, text="Plot", command=self.plot_q_profiles)
        plot_button.pack()

    def update_equilibrium_inputs(self):
        """Update widgets for equilibrium inputs based on user input."""
        num_eq = int(self.num_eq_entry.get())

        for entry in self.eq_pulse_entries + self.dda_entries + self.uid_entries + self.seq_entries:
            entry.destroy()

        self.eq_pulse_entries.clear()
        self.dda_entries.clear()
        self.uid_entries.clear()
        self.seq_entries.clear()

        for i in range(num_eq):
            eq_input_frame = tk.Frame(self.master)
            eq_input_frame.pack()

            pulse_label, pulse_entry = self.create_label_entry_pair(eq_input_frame, f"Pulse {i+1}:")
            self.eq_pulse_entries.append(pulse_entry)

            dda_label, dda_entry = self.create_label_entry_pair(eq_input_frame, "DDA:")
            self.dda_entries.append(dda_entry)

            uid_label, uid_entry = self.create_label_entry_pair(eq_input_frame, "UID:")
            self.uid_entries.append(uid_entry)

            seq_label, seq_entry = self.create_label_entry_pair(eq_input_frame, "Sequence:")
            self.seq_entries.append(seq_entry)

    def create_label_entry_pair(self, parent, label_text):
        """Create a pair of label and entry, and pack them into the parent frame."""
        label = tk.Label(parent, text=label_text)
        label.pack(side=tk.LEFT)
        entry = tk.Entry(parent)
        entry.pack(side=tk.LEFT)
        return label, entry

    def create_transp_inputs(self):
        """Create widgets for TRANSP inputs."""
        transp_frame = tk.Frame(self.master)
        transp_frame.pack()

        transp_label = tk.Label(transp_frame, text="TRANSP Runs")
        transp_label.pack()

        num_transp_label = tk.Label(transp_frame, text="TRANSP: Number of plotted pulses:")
        num_transp_label.pack(side=tk.LEFT)

        self.num_transp_entry = tk.Entry(transp_frame)
        self.num_transp_entry.pack(side=tk.LEFT)
        self.num_transp_entry.insert(0, "0")

        num_transp_button = tk.Button(transp_frame, text="Set", command=self.update_transp_inputs)
        num_transp_button.pack(side=tk.LEFT)

    def update_transp_inputs(self):
        """Update widgets for TRANSP inputs based on user input."""
        num_transp = int(self.num_transp_entry.get())

        for entry in self.transp_pulse_entries + self.runid_entries:
            entry.destroy()

        self.transp_pulse_entries.clear()
        self.runid_entries.clear()

        for i in range(num_transp):
            transp_input_frame = tk.Frame(self.master)
            transp_input_frame.pack()

            pulse_label, pulse_entry = self.create_label_entry_pair(transp_input_frame, f"Pulse {i+1}:")
            self.transp_pulse_entries.append(pulse_entry)

            runid_label, runid_entry = self.create_label_entry_pair(transp_input_frame, "RunIds:")
            self.runid_entries.append(runid_entry)

    def plot_q_profiles(self):
        try:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_facecolor('white')
            ax.grid(True, which='both', color='gray', linewidth=0.5)
            ax.tick_params(direction='out')
            ax.set_frame_on(True)
            try:
                time_slices = [float(t) for t in self.time_entry.get().split(",")]
            except ValueError:
                print("Invalid time input. Please input time as comma-separated numbers.")
                return

            try:
                num_eq = int(self.num_eq_entry.get())
            except ValueError:
                print("Invalid number of Equilibria. Please input an integer.")
                return

            try:
                num_transp = int(self.num_transp_entry.get())
                print(f'num_transp: {num_transp}')
            except ValueError:
                print("Invalid number of TRANSP Runs. Please input an integer.")
                return

            plot_type = self.plot_type_var.get()

            if plot_type not in ['Profile', 'Time Trace']:
                print("Invalid plot type. Please select either 'Profile' or 'Time Trace'.")
                return

            for i in range(num_eq):
                pulse = int(self.eq_pulse_entries[i].get())
                dda = self.dda_entries[i].get()
                uid = self.uid_entries[i].get()
                seq = int(self.seq_entries[i].get())

                eq = pr.Eq(pulse=pulse, dda=dda, uid=uid, seq=seq)
                self.equilibria.append(eq)

                for time in time_slices:
                    ti = np.abs(eq.t - time).argmin() if plot_type == 'Profile' else np.abs(eq.rntf - time).argmin(axis=1)
                    if plot_type == 'Profile':
                        ax.plot(eq.rntf[ti, :], eq.Q()[ti, :], label=f"{pulse}/{dda}/{uid}/{seq}, t={time}")
                    elif plot_type == 'Time Trace':
                        q = np.array([eq.Q()[i, ti[i]] for i in range(len(ti))])
                        ax.plot(eq.t, q, label=f"{pulse}/{dda}/{uid}/{seq}, X={time}")

            #linestyles = [ (0, ()), (0, (1, 1)),(0, (5, 10))), (0, (5, 5))), (0, (3, 5, 1, 5))), (0, (5, 1)))]

            linestyles_tuples = [
                ('solid', (0,())),
                ('densely dashed', (0, (5, 1))),
                ('dotted', (0, (1, 1))),
                ('dashdotted', (0, (3, 5, 1, 5))),
                ('densely dashdotted', (0, (3, 1, 1, 1))),
                ('long dash with offset', (5, (8, 3))),
                ('loosely dashed', (0, (5, 10))),


                ('loosely dashdotted', (0, (3, 10, 1, 10))),



                ('dashdotdotted', (0, (3, 5, 1, 5, 1, 5))),
                ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
                ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))]
            linestyles = [t[1] for t in linestyles_tuples]
            lnstyle=iter(linestyles)
            for i in range(num_transp):
                pulse = int(self.transp_pulse_entries[i].get())
                runids = self.runid_entries[i].get().split(',')
                for runid in runids:
                    runid = runid.strip()
                    self.transp_runs[runid] = pr.Transp(pulse, runid)
                    self.transp_runs[runid].add_data('Q')

                    for time in time_slices:
                        ln = next(lnstyle)
                        ti = get_index(self.transp_runs[runid].t, time-40.) if plot_type == 'Profile' else np.abs(self.transp_runs[runid].x - time).argmin(axis=1)
                        if plot_type == 'Profile':
                            ax.plot(
                                self.transp_runs[runid].x[ti, :],
                                self.transp_runs[runid].transp('Q')[ti, :],
                                label=f"{pulse}{runid}, t={time}",
                                color = 'k',
                                linestyle = ln
                            )
                            if num_eq==0 and i==0:
                                ax.plot([0, 1], [1, 1], linestyle='dashed', color='black', linewidth=0.5)
                        elif plot_type == 'Time Trace':
                            q = np.array([self.transp_runs[runid].transp('Q')[i, ti[i]] for i in range(len(ti))])
                            ax.plot(self.transp_runs[runid].t+40., q, color = 'k', linestyle=ln, label=f"{pulse}{runid}, x={time}")
                            if num_eq==0 and i==0:
                                ax.plot([self.transp_runs[runid].t[0]+40., self.transp_runs[runid].t[-1]+40.], [1, 1], linestyle='dashed', color='black', linewidth=0.5)

            if plot_type == 'Profile':
                ax.plot([0, 1], [1, 1], linestyle='dashed', color='black', linewidth=0.5)
            else:
                if num_eq>0:
                    ax.plot([eq.t[0], eq.t[-1]], [1, 1], linestyle='dashed', color='black', linewidth=0.5)

            ax.set_xlabel('X' if plot_type == 'Profile' else 'Time')
            ax.set_ylabel('Q')
            ax.set_ylim(bottom=0.0)
            ax.legend()

            plt.show()

            self.equilibria.clear()
            self.transp_runs.clear()
        except Exception as e:
            print(f"An error occurred: {e}")

root = tk.Tk()
gui = QProfileGUI(root)
root.mainloop()

