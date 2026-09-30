"""EME taper transfer for Fig. 2(b,c), check V9 and make_numbers.py.

Runs waveguide.eme_taper() with the arguments of the archived run
(Zenodo 10.5281/zenodo.22758086, data/taper_eme.json and eme_log2.txt):
TFLN bus width 0.6 um, diamond taper 0.35 -> 0.05 um in 30 slices,
4 modes per slice, taper lengths 1-20 um.  Writes results/taper_eme.json
(written by eme_taper) and results/eme_log2.txt (the printed log, which
make_figures.py reads for the supermode branches of Fig. 2(b))."""
import sys

import waveguide as wg

W_LN = 0.6                                   # TFLN bus width (um)
LENGTHS = [1, 2, 3, 4, 6, 8, 10, 12, 16, 20]  # taper lengths (um)
NUM_MODES = 4                                # modes kept per slice


class Tee:
    """Print to the console and to the log file at the same time."""
    def __init__(self, *streams):
        self.streams = streams

    def write(self, s):
        for st in self.streams:
            st.write(s)

    def flush(self):
        for st in self.streams:
            st.flush()


with open(wg.RESULTS / "eme_log2.txt", "w") as log:
    out, sys.stdout = sys.stdout, Tee(sys.stdout, log)
    try:
        wg.eme_taper(W_LN, w0=0.35, w1=0.05, lengths=LENGTHS, n_seg=30,
                     num_modes=NUM_MODES)
        print("EME DONE")
    finally:
        sys.stdout = out
