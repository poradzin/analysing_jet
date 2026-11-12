#!/usr/bin/env python3
"""
analyse_JET_pulse.py -- check which diagnostics are available for a JET pulse
before fitting profiles for TRANSP input.

For a given pulse the script
  1. finds when the heating (NBI, ICRH) is active -> heating window,
  2. sets the analysis window (heating window by default, can be overridden
     or widened from the command line),
  3. for every selected diagnostic (PPF signal) reports whether it exists,
     its full time range, and inside the window: first/last time slice,
     number of slices, average interval and frequency,
  4. optionally plots the heating, a timeline of all measurement times and
     each diagnostic.

Examples
  analyse_JET_pulse.py 103073                      # heating + all diagnostics
  analyse_JET_pulse.py --list                      # show known diagnostics
  analyse_JET_pulse.py 103073 -d HRTS              # all HRTS signals
  analyse_JET_pulse.py 103073 -d HRTS/NE CXD6/TI   # chosen signals
  analyse_JET_pulse.py 103073 -d KK3/TE30 CXD6/TI/cxsbatch/2   # ad-hoc signals
  analyse_JET_pulse.py 103073 -g ne Ti             # by group
  analyse_JET_pulse.py 103073 --tstart 44.98 --tend 46.38
  analyse_JET_pulse.py 103073 --before 0.2 --after 0.1  # widen heating window
  analyse_JET_pulse.py 103073 -d HRTS/NE --plot --x 3.2
  analyse_JET_pulse.py 103073 --save figs          # write PNGs (no display)

To add a diagnostic permanently, append a Signal(...) line to DIAGNOSTICS
below. Any signal can also be given ad hoc with -d DDA/DTYPE[/UID[/SEQ]].
"""
import argparse
import os
import sys
from dataclasses import dataclass, field

import numpy as np

try:
    import ppf
except ImportError:  # allows --list / --help outside JET
    ppf = None


# --------------------------------------------------------------------------
# Signal definitions
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Signal:
    dda: str
    dtype: str
    uid: str = 'jetppf'
    seq: int = 0
    group: str = ''

    @property
    def name(self):
        return f'{self.dda}/{self.dtype}/{self.uid}/{self.seq}'

    @classmethod
    def parse(cls, text, group='user'):
        """Build a Signal from 'DDA/DTYPE[/UID[/SEQ]]'."""
        parts = text.split('/')
        if not 2 <= len(parts) <= 4 or not all(parts[:2]):
            raise ValueError(f"cannot parse signal '{text}', "
                             "expected DDA/DTYPE[/UID[/SEQ]]")
        uid = parts[2].lower() if len(parts) > 2 and parts[2] else 'jetppf'
        seq = int(parts[3]) if len(parts) > 3 and parts[3] else 0
        return cls(parts[0].upper(), parts[1].upper(), uid, seq, group)

    def matches(self, token):
        """True if token ('HRTS', 'CXD6/TI', 'CXD6/TI/cxsbatch', ...)
        matches this signal component by component (case-insensitive)."""
        mine = [self.dda, self.dtype, self.uid, str(self.seq)]
        parts = [p for p in token.split('/')]
        return all(p == '' or p.lower() == m.lower()
                   for p, m in zip(parts, mine))


# Diagnostics checked by default. Add new ones here.
# group is used by --group and for ordering of the output.
DIAGNOSTICS = [
    # electron density
    Signal('HRTS', 'NE', group='ne'),
    Signal('LIDR', 'NE', group='ne'),
    # electron temperature
    Signal('HRTS', 'TE', group='Te'),
    Signal('KK3', 'TE21', group='Te'),
    Signal('KK1', 'TE21', group='Te'),
    # ion temperature
    Signal('CXD6', 'TI', group='Ti'),
    Signal('CXG6', 'TI', group='Ti'),
    Signal('CXD6', 'TIFS', group='Ti'),
    Signal('CXG6', 'TIFS', group='Ti'),
    Signal('CXDM', 'TIFS', group='Ti'),
    Signal('CXD6', 'TICR', group='Ti'),
    Signal('CXG6', 'TICR', group='Ti'),
    Signal('CXD6', 'TI', 'cxsbatch', group='Ti'),
    Signal('CXG6', 'TI', 'cxsbatch', group='Ti'),
    Signal('CXD6', 'TICR', 'cxsbatch', group='Ti'),
    Signal('CXG6', 'TICR', 'cxsbatch', group='Ti'),
    Signal('CX7C', 'TI', group='Ti'),
    Signal('CX7D', 'TI', group='Ti'),
    Signal('TION', 'TI', group='Ti'),
    Signal('CXKM', 'TI', group='Ti'),
    Signal('CXFM', 'TICR', group='Ti'),
    # rotation
    Signal('CXD6', 'AFCR', group='rotation'),
    Signal('CXG6', 'AFCR', group='rotation'),
]

# Heating systems used to define the heating window.
HEATING = [
    Signal('NBI', 'PTOT', group='heating'),
    Signal('NBI4', 'PTOT', group='heating'),
    Signal('NBI8', 'PTOT', group='heating'),
    Signal('ICRH', 'PTOT', group='heating'),
]
NBI_OCTANTS = ('NBI4', 'NBI8')   # GAS dtype is read for the beam species
HEATING_COLOURS = {'NBI': 'blue', 'NBI4': 'green', 'NBI8': 'orange',
                   'ICRH': 'red'}


# --------------------------------------------------------------------------
# Reading data
# --------------------------------------------------------------------------
@dataclass
class SignalData:
    signal: Signal
    t: np.ndarray            # (nt,)
    x: np.ndarray            # (nx,) or empty for a pure time trace
    data: np.ndarray         # (nt, max(nx,1))
    dunits: str = ''
    xunits: str = ''
    desc: str = ''
    seq_read: int = 0
    stats: dict = field(default_factory=dict)


def read_signal(pulse, sig):
    """Return (SignalData or None, status string)."""
    try:
        out = ppf.ppfdata(pulse, sig.dda, sig.dtype, seq=sig.seq, uid=sig.uid)
    except Exception as exc:  # ppf library errors
        return None, f'read error: {exc}'
    value, x, t, nd, nx, nt, dunits, xunits, tunits, desc, comm, seq, ier = out
    if ier != 0 or len(t) == 0:
        return None, 'not available' + (f' (ier={ier})' if ier else '')

    t = np.asarray(t, dtype=float)
    x = np.asarray(x, dtype=float) if len(x) > 0 else np.array([])
    value = np.asarray(value, dtype=float)
    ncol = max(len(x), 1)
    if value.size != len(t) * ncol:
        return None, f'unexpected size {value.size} != {len(t)}x{ncol}'
    return SignalData(sig, t, x, value.reshape(len(t), ncol),
                      str(dunits).strip(), str(xunits).strip(),
                      str(desc).strip(), int(seq)), 'OK'


# --------------------------------------------------------------------------
# Heating
# --------------------------------------------------------------------------
def check_heating(pulse, cutoff=0.01, verbose=True):
    """Find when each heating system is above `cutoff` [W].

    Returns (heating dict, (t_on, t_off) or None)."""
    heating = {}
    for sig in HEATING:
        sd, status = read_signal(pulse, sig)
        if sd is None:
            if verbose:
                print(f'  {sig.dda:<5} {status}')
            continue
        p = sd.data[:, 0]
        on = p > cutoff
        if not on.any():
            if verbose:
                print(f'  {sig.dda:<5} present but never above cutoff')
            continue
        h = dict(t=sd.t, p=p, t_on=sd.t[on][0], t_off=sd.t[on][-1],
                 pmax=p.max(), species='')
        if sig.dda in NBI_OCTANTS:
            try:
                gas = ppf.ppfdata(pulse, sig.dda, 'GAS', seq=0, uid='jetppf')
                if len(gas[0]) > 0 and gas[9].split():
                    h['species'] = gas[9].split()[-1]
            except Exception:
                pass
        heating[sig.dda] = h
        if verbose:
            sp = f'  species: {h["species"]}' if h['species'] else ''
            print(f'  {sig.dda:<5} on {h["t_on"]:8.4f} - {h["t_off"]:8.4f} s'
                  f'   Pmax {h["pmax"] / 1e6:6.2f} MW{sp}')

    if not heating:
        return heating, None
    t_on = min(h['t_on'] for h in heating.values())
    t_off = max(h['t_off'] for h in heating.values())
    return heating, (t_on, t_off)


# --------------------------------------------------------------------------
# Statistics
# --------------------------------------------------------------------------
def window_stats(sd, t0, t1):
    t = sd.t
    mask = (t >= t0) & (t <= t1)
    tw = t[mask]
    dw = sd.data[mask]
    n = len(tw)
    s = dict(n=n, t_first=np.nan, t_last=np.nan, dt=np.nan, freq=np.nan,
             dt_min=np.nan, dt_max=np.nan, dt_std=np.nan,
             n_valid=int(np.sum(np.any(np.isfinite(dw) & (dw != 0), axis=1))),
             before=t[t < t0][-1] if np.any(t < t0) else np.nan,
             after=t[t > t1][0] if np.any(t > t1) else np.nan)
    if n:
        s['t_first'], s['t_last'] = tw[0], tw[-1]
    if n > 1:
        dt = np.diff(tw)
        s.update(dt=dt.mean(), freq=1.0 / dt.mean(), dt_min=dt.min(),
                 dt_max=dt.max(), dt_std=dt.std())
    return s


def print_table(rows, t0, t1):
    w = max([len(r[0].name) for r in rows] + [6]) + 2
    print(f'\nDiagnostics in window {t0:.4f} - {t1:.4f} s '
          f'({t1 - t0:.4f} s)\n')
    head = (f'{"Signal":<{w}}{"nx":>4}  {"full range [s]":>17}  '
            f'{"first [s]":>9} {"last [s]":>9} {"N":>5} {"valid":>5} '
            f'{"<dt> [ms]":>9} {"f [Hz]":>8}')
    print(head)
    print('-' * len(head))
    group = None
    for sig, sd, status in rows:
        if sig.group != group:
            group = sig.group
            if group:
                print(f'[{group}]')
        if sd is None:
            print(f'{sig.name:<{w}}  -- {status}')
            continue
        s = sd.stats
        full = f'{sd.t[0]:7.3f} - {sd.t[-1]:7.3f}'
        if s['n'] == 0:
            print(f'{sig.name:<{w}}{len(sd.x):>4}  {full:>17}  '
                  f'no slices in window')
            continue
        print(f'{sig.name:<{w}}{len(sd.x):>4}  {full:>17}  '
              f'{s["t_first"]:9.4f} {s["t_last"]:9.4f} {s["n"]:5d} '
              f'{s["n_valid"]:5d} {s["dt"] * 1e3:9.2f} {s["freq"]:8.1f}')
    print()


def print_details(rows):
    print('Detailed interval statistics')
    for sig, sd, _ in rows:
        if sd is None or sd.stats['n'] == 0:
            continue
        s = sd.stats
        info = ', '.join(i for i in (f'seq read: {sd.seq_read}', sd.desc,
                                     f'[{sd.dunits}]') if i)
        print(f'  {sig.name} ({info})')
        if s['n'] > 1:
            print(f'    dt min/max/std: {s["dt_min"] * 1e3:.3f} / '
                  f'{s["dt_max"] * 1e3:.3f} / {s["dt_std"] * 1e3:.3f} ms')
        print(f'    nearest slice before window: {s["before"]:.4f} s, '
              f'after window: {s["after"]:.4f} s')
    print()


# --------------------------------------------------------------------------
# Plotting
# --------------------------------------------------------------------------
class FigureSink:
    """Collect figures: show them in tabs (plotWindow) and/or save PNGs."""

    def __init__(self, show, save_dir, pulse):
        self.show_figs = show
        self.save_dir = save_dir
        self.pulse = pulse
        self.win = None
        import matplotlib
        if show:
            try:
                import plotWindow  # sets the qt5agg backend
                self.win = plotWindow.plotWindow()
                self.win.MainWindow.setWindowTitle(f'JET pulse {pulse}')
            except Exception as exc:
                print(f'plotWindow unavailable ({exc}); '
                      'using plain matplotlib windows')
        else:
            matplotlib.use('Agg')
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)

    def add(self, title, fig):
        if self.save_dir:
            fname = f'{self.pulse}_{title.replace("/", "_")}.png'
            fig.savefig(os.path.join(self.save_dir, fname), dpi=120)
        if self.win is not None:
            self.win.addPlot(title, fig)

    def show(self):
        if self.win is not None:
            self.win.show()
        elif self.show_figs:
            import matplotlib.pyplot as plt
            plt.show()


def shade(ax, window, heat_window):
    if heat_window:
        ax.axvspan(*heat_window, color='orange', alpha=0.12,
                   label='heating')
    ax.axvline(window[0], color='k', ls='--', lw=1)
    ax.axvline(window[1], color='k', ls='--', lw=1)


def plot_heating(heating, window, heat_window, pulse):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    for dda, h in heating.items():
        label = f'{dda}: {h["species"]}' if h['species'] else dda
        ax.plot(h['t'], h['p'] / 1e6, color=HEATING_COLOURS.get(dda),
                ls='--' if dda == 'NBI' else '-', label=label)
    shade(ax, window, heat_window)
    ax.set_xlabel('time [s]')
    ax.set_ylabel('P [MW]')
    ax.set_title(f'Heating {pulse}')
    ax.legend()
    return fig


def plot_timeline(rows, window, heat_window, pulse):
    import matplotlib.pyplot as plt
    avail = [(sig, sd) for sig, sd, _ in rows if sd is not None]
    fig, ax = plt.subplots(figsize=(8, 1.0 + 0.35 * len(avail)))
    for i, (sig, sd) in enumerate(avail):
        inside = (sd.t >= window[0]) & (sd.t <= window[1])
        ax.eventplot(sd.t[~inside], lineoffsets=i, linelengths=0.6,
                     colors='lightgray')
        ax.eventplot(sd.t[inside], lineoffsets=i, linelengths=0.6,
                     colors=f'C{i % 10}')
    shade(ax, window, heat_window)
    ax.set_yticks(range(len(avail)))
    ax.set_yticklabels([f'{sig.name} ({sd.stats["n"]})'
                        for sig, sd in avail], fontsize=8)
    ax.invert_yaxis()
    pad = 0.25 * (window[1] - window[0])
    ax.set_xlim(window[0] - pad, window[1] + pad)
    ax.set_xlabel('time [s]')
    ax.set_title(f'Measurement times {pulse} (N in window)')
    fig.tight_layout()
    return fig


def plot_signal(sd, window, heat_window, channel=None, xpos=None):
    import matplotlib.pyplot as plt
    sig, s = sd.signal, sd.stats
    if xpos is not None and len(sd.x) > 1:
        ich = int(np.nanargmin(np.abs(sd.x - xpos)))
    else:
        ich = min(channel or 0, sd.data.shape[1] - 1)
    pad = 0.25 * (window[1] - window[0])
    xlim = (window[0] - pad, window[1] + pad)
    profile = len(sd.x) > 1

    if profile:
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True,
                                       figsize=(7, 7))
        view = (sd.t >= xlim[0]) & (sd.t <= xlim[1])
        if view.sum() > 1:
            pc = ax1.pcolormesh(sd.t[view], sd.x, sd.data[view].T,
                                shading='nearest')
            fig.colorbar(pc, ax=ax1, label=sd.dunits)
            # invisible twin colorbar keeps the two time axes aligned
            fig.colorbar(pc, ax=ax2).ax.set_visible(False)
        ax1.axhline(sd.x[ich], color='w', lw=0.8, ls=':')
        ax1.axvline(window[0], color='w', ls='--', lw=1)
        ax1.axvline(window[1], color='w', ls='--', lw=1)
        ax1.set_ylabel(f'x [{sd.xunits}]')
    else:
        fig, ax2 = plt.subplots(figsize=(7, 4.5))

    y = sd.data[:, ich]
    inside = (sd.t >= window[0]) & (sd.t <= window[1])
    ax2.plot(sd.t, y, color='gray', lw=0.8)
    ax2.plot(sd.t[inside], y[inside], 'o', ms=3, color='C0')
    shade(ax2, window, heat_window)
    ax2.set_xlim(xlim)
    ax2.set_xlabel('time [s]')
    xlab = f', x={sd.x[ich]:.3f} {sd.xunits}' if profile else ''
    ax2.set_ylabel(f'{sig.dtype} [{sd.dunits}] (ch {ich}{xlab})')

    if s['n'] > 0:
        txt = '\n'.join([
            f'Window: ({s["t_first"]:.3f}, {s["t_last"]:.3f}) s',
            f'No. of slices: {s["n"]}',
        ] + ([f'f = {s["freq"]:.1f} Hz',
              r'$\Delta t$ = ' + f'{s["dt"] * 1e3:.2f} ms']
             if s['n'] > 1 else []))
        ax2.text(0.98, 0.95, txt, transform=ax2.transAxes, fontsize=9,
                 va='top', ha='right',
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    fig.suptitle(f'{sig.name} (seq {sd.seq_read})')
    fig.tight_layout()
    return fig


# --------------------------------------------------------------------------
# Command line
# --------------------------------------------------------------------------
def select_signals(tokens, groups):
    selected = list(DIAGNOSTICS)
    if groups:
        g = {x.lower() for x in groups}
        selected = [s for s in selected if s.group.lower() in g]
    if tokens:
        chosen = []
        for tok in tokens:
            hits = [s for s in selected if s.matches(tok)]
            if not hits:  # unknown signal -> ad-hoc
                hits = [Signal.parse(tok)]
            chosen += [h for h in hits if h not in chosen]
        selected = chosen
    return selected


def list_diagnostics():
    print('Known diagnostics (edit DIAGNOSTICS in this script to add more):')
    group = None
    for s in DIAGNOSTICS:
        if s.group != group:
            group = s.group
            print(f'\n  [{group}]')
        print(f'    {s.name}')
    print('\nHeating: ' + ', '.join(s.name for s in HEATING))


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description='Check heating and diagnostic availability / time '
                    'slices for a JET pulse.',
        epilog='Run with -h and see the module docstring for examples.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('pulse', type=int, nargs='?', help='JET pulse number')
    p.add_argument('-l', '--list', action='store_true',
                   help='list known diagnostics and exit')
    p.add_argument('-d', '--diag', nargs='+', metavar='SIG',
                   help='diagnostics to check: DDA, DDA/DTYPE or '
                        'DDA/DTYPE/UID/SEQ; unknown ones are read ad hoc '
                        '(default: all known)')
    p.add_argument('-g', '--group', nargs='+',
                   help='only these groups, e.g. ne Te Ti rotation')
    w = p.add_argument_group('time window (default: heating window)')
    w.add_argument('--tstart', type=float, help='window start [s]')
    w.add_argument('--tend', type=float, help='window end [s]')
    w.add_argument('--before', type=float, default=0.0,
                   help='extend window start earlier by this many s')
    w.add_argument('--after', type=float, default=0.0,
                   help='extend window end later by this many s')
    w.add_argument('--power-cutoff', type=float, default=0.01,
                   help='heating considered on above this power [W]')
    w.add_argument('--no-heating', action='store_true',
                   help='skip heating check (window from --tstart/--tend '
                        'or full signal range)')
    o = p.add_argument_group('output')
    o.add_argument('-p', '--plot', action='store_true',
                   help='plot heating, timeline and each diagnostic')
    o.add_argument('--save', metavar='DIR',
                   help='save the plots as PNG files into DIR')
    o.add_argument('--channel', type=int, default=0,
                   help='profile channel index plotted as time trace')
    o.add_argument('--x', type=float,
                   help='plot time trace at channel nearest to this x '
                        '(overrides --channel)')
    o.add_argument('-s', '--stats', action='store_true',
                   help='print detailed interval statistics')
    o.add_argument('-q', '--quiet', action='store_true',
                   help='hide unavailable signals in the table')
    args = p.parse_args(argv)
    if args.pulse is None and not args.list:
        p.error('pulse number is required')
    return args


def main(argv=None):
    args = parse_args(argv)
    if args.list:
        list_diagnostics()
        return 0
    if ppf is None:
        sys.exit('Error: ppf module not found (run on a JET machine).')

    pulse = args.pulse
    print(f'PULSE {pulse}\n')
    ier = ppf.ppfgo(pulse, 0)
    if ier != 0:
        sys.exit(f'Error initialising PPF for pulse {pulse} (ier={ier}).')

    try:
        signals = select_signals(args.diag, args.group)
    except ValueError as exc:
        sys.exit(f'Error: {exc}')
    if not signals:
        sys.exit('No diagnostics selected.')

    heating, heat_window = {}, None
    if not args.no_heating:
        print('Heating')
        heating, heat_window = check_heating(pulse, args.power_cutoff)
        if heat_window:
            print(f'Heating window: {heat_window[0]:.4f} - '
                  f'{heat_window[1]:.4f} s\n')
        else:
            print('No active heating found.\n')

    rows = []
    for sig in signals:
        sd, status = read_signal(pulse, sig)
        rows.append((sig, sd, status))
    available = [sd for _, sd, _ in rows if sd is not None]

    # analysis window: explicit > heating > full data range
    if heat_window:
        t0, t1 = heat_window
    elif available:
        t0 = min(sd.t[0] for sd in available)
        t1 = max(sd.t[-1] for sd in available)
    else:
        t0 = t1 = None
    t0 = args.tstart if args.tstart is not None else t0
    t1 = args.tend if args.tend is not None else t1
    if t0 is None or t1 is None:
        print('No data available for the selected diagnostics.')
        return 1
    t0, t1 = t0 - args.before, t1 + args.after
    if t1 <= t0:
        sys.exit(f'Error: empty window ({t0}, {t1}).')

    for sd in available:
        sd.stats = window_stats(sd, t0, t1)

    shown = [r for r in rows if r[1] is not None] if args.quiet else rows
    print_table(shown, t0, t1)
    if args.stats:
        print_details(rows)

    if args.plot or args.save:
        sink = FigureSink(args.plot, args.save, pulse)
        window = (t0, t1)
        if heating:
            sink.add('heating', plot_heating(heating, window, heat_window,
                                             pulse))
        if available:
            sink.add('timeline', plot_timeline(rows, window, heat_window,
                                               pulse))
        for sd in available:
            sink.add(sd.signal.name,
                     plot_signal(sd, window, heat_window, args.channel,
                                 args.x))
        if args.save:
            print(f'Plots saved in {args.save}/')
        sink.show()
    return 0


if __name__ == '__main__':
    sys.exit(main())
