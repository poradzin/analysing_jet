import ppf
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse

parser = argparse.ArgumentParser(description='Read KT5P tritium data from PPF.')
parser.add_argument('--save', action='store_true', help='Save results to file')
parser.add_argument('--method', choices=['max', 'median', 'p95','sym_p95'], default='max',
                    help='Method for summarising each pulse (default: max)')
args = parser.parse_args()

p_start = 95000
p_stop  = 105900

pulse = []
tttd  = []
dthd  = []
hthd  = []

for pul in range(p_start, p_stop):
    PPFseq = 0
    user = 'JETPPF'
    ier = ppf.ppfgo(pul, PPFseq)
    if ier != 0:
        raise Exception(f'ier = {ier}. Error initialising PPF routines. Aborting.')
    ier = ppf.ppfuid(user, rw="R")

    data_t = ppf.ppfdata(pul, 'KT5P', 'TTTD')[0]
    data_d = ppf.ppfdata(pul, 'KT5P', 'DTHD')[0]
    data_h = ppf.ppfdata(pul, 'KT5P', 'HTHD')[0]

    print(f'Pulse {pul}, np.size(data_t): {np.size(data_t)}')

    # Step 1: mask unphysical values
    data = np.stack([data_t, data_d, data_h])
    mask = np.all((data >= 0) & (data <= 1), axis=0)
    data_t, data_d, data_h = data[:, mask]

    if np.size(data_t) > 0:
        # Step 2/4: summarise according to chosen method
        if args.method == 'max':
            idx = np.argmax(data_t)
        elif args.method == 'median':
            idx = np.argmin(np.abs(data_t - np.median(data_t)))
        elif args.method == 'p95':
            threshold = np.percentile(data_t, 95)
            idx = np.argmin(np.abs(data_t - threshold))
        elif args.method == 'sym_p95':
            magnitude = data_t**2 + data_d**2 + data_h**2
            threshold = np.percentile(magnitude, 95)
            idx = np.argmin(np.abs(magnitude - threshold))

        pulse.append(pul)
        tttd.append(data_t[idx])
        dthd.append(data_d[idx])
        hthd.append(data_h[idx])
    else:
        print(f'Pulse {pul}: skipping')

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.plot(pulse, tttd, marker='o')
ax1.set_ylabel('Max tritium')
ax1.grid(True)
ax2.plot(pulse, hthd, marker='o', color='orange')
ax2.set_ylabel('Hydrogen')
ax2.set_xlabel('Pulse')
ax2.grid(True)
plt.suptitle(f'Method: {args.method}')
plt.tight_layout()
plt.show()

if args.save:
    os.makedirs("outputs", exist_ok=True)
    np.savetxt("outputs/HDT_isotope_levels.txt",
               np.column_stack([pulse, tttd, dthd, hthd]),
               fmt='%d %.6e %.6e %.6e',
               header=f'pulse  max_tritium  deuterium  hydrogen  [method: {args.method}]')
    print("Saved to outputs/H_isotope_levels.txt")
