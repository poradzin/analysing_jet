import argparse
import ppf
import scipy
import numpy as np
from scipy.interpolate import interp1d, interp2d
import profiles as pr
import plotWindow as pw
import matplotlib.pyplot as plt

def gi(x, y):
    return (np.abs(x - y)).argmin()

class Equil(pr.Eq):
    def __init__(self, pulse, dda='EFTP', uid='jetppf', seq=0):
        super().__init__(pulse, dda=dda, uid=uid, seq=seq)
  
    def Q(self):
        return self._Q.reshape((len(self._t), len(self._x)))
  
    @property
    def rntf(self):
        '''
        Square root of normalized toroidal flux
        '''   
        return self._sqrt_ftor_norm


def plot_Q_profiles(win, equilibria, t):
    fig = plt.figure() 
    fig.suptitle(f'Q profiles at t = {t:.2f}', fontsize=13) 
    ax = fig.add_subplot(111)

    cmap = plt.cm.get_cmap('tab10')
    color = iter(cmap.colors)

    for eq in equilibria:
        ti = gi(t, eq.t)
        c = next(color)
        ax.plot(
            eq.rntf[ti, :],
            eq.Q()[ti, :],
            color=c,
            linewidth=1
        )

    ax.plot([0, 1], [1, 1], linestyle='dashed', color='black')
    ax.set_xlabel('X')
    ax.set_xlim(0.0, 1.0)

    win.addPlot('Q profiles', fig)

def main(eq_configs, time=None):
    win = pw.plotWindow()
    equilibria = []

    for config in eq_configs:
        eq = Equil(*config)
        equilibria.append(eq)

    if time is not None:
        plot_Q_profiles(win, equilibria, time)

    win.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Q-Profile Script")
    parser.add_argument(
        "-eq",
        "--equilibria",
        nargs="+",
        metavar=("pulse", "dda", "uid", "seq"),
        help="Equilibrium configurations: pulse dda uid seq"
    )
    parser.add_argument("-t", "--time", type=float, default=None, help="Requested time slice (default: None)")

    args = parser.parse_args()

    eq_configs = [list(map(int, eq.split())) for eq in args.equilibria]
    main(eq_configs, args.time)


