# -*- coding: utf-8 -*-
"""Performance-versus-NFE statistics of the equal-NFE section of the paper, from the 15,000-evaluation
traces: (i) evaluations after which EGSO's median best-so-far comes within
1% of its final median; (ii) at 2,000 evaluations, against how many
competitors EGSO's median is already shorter than the competitor's final
median. Usage: analyze_perf_nfe.py RESULT_DIR [RESULT_DIR ...]"""
import sys
from pathlib import Path
import numpy as np
import scipy.io as sio

for d in sys.argv[1:]:
    d = Path(d)
    print('==', d)
    for m in range(1, 6):
        tr = []
        for f in sorted(d.glob(f'EqNFE_map{m}_fe15000_r*.mat')):
            S = sio.loadmat(f, squeeze_me=True)
            algos = [str(a) for a in S['Algos']]
            T = np.asarray(S['Trace'], dtype=float)
            if T.ndim == 2:
                T = T[None]
            tr.append(T)
            ck = int(S['ckpt'])
        TR = np.concatenate(tr, axis=0)                 # runs x algos x checkpoints
        med = np.nanmedian(TR, axis=0)                  # algos x checkpoints
        e = algos.index('EGSO')
        final = med[e, -1]
        k1 = int(np.argmax(med[e] <= 1.01 * final))
        k2000 = 2000 // ck - 1
        comp = [a for a in range(len(algos)) if a != e]
        beaten = sum(med[e, k2000] < med[a, -1] for a in comp)
        print(f'Map {m}: runs {TR.shape[0]}, final median {final:.2f}, within 1% after {(k1 + 1) * ck} NFE, '
              f'median at 2,000 = {med[e, k2000]:.2f} shorter than the final median of {beaten}/{len(comp)} competitors')
