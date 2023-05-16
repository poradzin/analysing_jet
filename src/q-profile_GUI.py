import tkinter as tk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import profiles as pr
import numpy as np

def gi(x, y):
    return (np.abs(x - y)).argmin()

class QProfileGUI:
    def __init__(self, master):
        self.master = master
        master.title("Q Profile Plotter")
        
        self.equilibria = []
        
        self.pulse_entries = []
        self.dda_entries = []
        self.uid_entries = []
        self.seq_entries = []
        
        self.time_entry = None
        self.plot_type_var = tk.StringVar()  # Variable to store the selected plot type
        
        self.create_equilibrium_inputs()
        self.create_plot_type_input()
        self.create_time_input()
        self.create_plot_button()
        
    def create_equilibrium_inputs(self):
        eq_frame = tk.Frame(self.master)
        eq_frame.pack()
        
        eq_label = tk.Label(eq_frame, text="Equilibria")
        eq_label.pack()
        
        num_eq_label = tk.Label(eq_frame, text="Number of Equilibria:")
        num_eq_label.pack(side=tk.LEFT)
        
        self.num_eq_entry = tk.Entry(eq_frame)
        self.num_eq_entry.pack(side=tk.LEFT)
        self.num_eq_entry.insert(0, "1")
        
        num_eq_button = tk.Button(eq_frame, text="Set", command=self.update_equilibrium_inputs)
        num_eq_button.pack(side=tk.LEFT)
        
    def create_time_input(self):
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
        plot_button = tk.Button(self.master, text="Plot", command=self.plot_q_profiles)
        plot_button.pack()
        
    def update_equilibrium_inputs(self):
        num_eq = int(self.num_eq_entry.get())
        
        for entry in self.pulse_entries:
            entry.destroy()
        self.pulse_entries.clear()
        
        for entry in self.dda_entries:
            entry.destroy()
        self.dda_entries.clear()
        
        for entry in self.uid_entries:
            entry.destroy()
        self.uid_entries.clear()
        
        for entry in self.seq_entries:
            entry.destroy()
        self.seq_entries.clear()
        
        for i in range(num_eq):
            eq_input_frame = tk.Frame(self.master)
            eq_input_frame.pack()
            
            pulse_label = tk.Label(eq_input_frame, text=f"Pulse {i+1}:")
            pulse_label.pack(side=tk.LEFT)
            pulse_entry = tk.Entry(eq_input_frame)
            pulse_entry.pack(side=tk.LEFT)
            self.pulse_entries.append(pulse_entry)

            dda_label = tk.Label(eq_input_frame, text="DDA:")
            dda_label.pack(side=tk.LEFT)
            dda_entry = tk.Entry(eq_input_frame)
            dda_entry.pack(side=tk.LEFT)
            self.dda_entries.append(dda_entry)
            
            uid_label = tk.Label(eq_input_frame, text="UID:")
            uid_label.pack(side=tk.LEFT)
            uid_entry = tk.Entry(eq_input_frame)
            uid_entry.pack(side=tk.LEFT)
            self.uid_entries.append(uid_entry)
            
            seq_label = tk.Label(eq_input_frame, text="Sequence:")
            seq_label.pack(side=tk.LEFT)
            seq_entry = tk.Entry(eq_input_frame)
            seq_entry.pack(side=tk.LEFT)
            self.seq_entries.append(seq_entry)
            
    def plot_q_profiles(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        time_slices = self.time_entry.get().split(",")
        num_eq = int(self.num_eq_entry.get())
        plot_type = self.plot_type_var.get()  # Get the selected plot type
        
        for i in range(num_eq):
            pulse = int(self.pulse_entries[i].get())
            dda = self.dda_entries[i].get()
            uid = self.uid_entries[i].get()
            seq = int(self.seq_entries[i].get())
            
            eq = pr.Eq(pulse=pulse, dda=dda, uid=uid, seq=seq)
            self.equilibria.append(eq)
            
            for time in time_slices:
                time = float(time.strip())
                ti = np.abs(eq.t - time).argmin() if plot_type == 'Profile' else np.abs(eq.rntf - time).argmin(axis=1)
                
                if plot_type == 'Profile':
                    ax.plot(eq.rntf[ti, :], eq.Q()[ti, :], label=f"{pulse}/{dda}/{uid}/{seq}, t={time}")
                elif plot_type == 'Time Trace':
                    q = np.array([eq.Q()[i, ti[i]] for i in range(len(ti))])
                    ax.plot(eq.t, q, label=f"{pulse}/{dda}/{uid}/{seq}, rntf={time}")
        
        if plot_type == 'Profile':
            ax.plot([0, 1], [1, 1], linestyle='dashed', color='black', linewidth=1)
        else:
            ax.plot([eq.t[0], eq.t[-1]], [1, 1], linestyle='dashed', color='black', linewidth=1)

        ax.set_xlabel('X' if plot_type == 'Profile' else 'Time')
        ax.set_ylabel('Q')
        ax.set_ylim(bottom=0.0)
        ax.legend()
        
        plt.show()  # Show the plot in a pop-up window
        
        self.equilibria.clear()

root = tk.Tk()
gui = QProfileGUI(root)
root.mainloop()

