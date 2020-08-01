import argparse
import logging
import os
import shutil
import warnings

from mmgpbsa.amber_mmgpbsa import run_amber
from mmgpbsa.systemloader import SystemLoader


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--amber_path', type=str, default=None, required=False,
                        help='path to amber installation folder (ex /home/austin/amber20)')
    parser.add_argument('--pdb', type=str,
                        help='pdb input file with ligand and protein in complex. Only one chain for now please!')
    parser.add_argument('--platform', type=str, choices=['CPU', 'CUDA', 'OpenCL'], default=None)
    ##simulation options
    parser.add_argument('--ps', type=float, default=None, required=False, help='picoseconds to run simulation')
    parser.add_argument('--equil_ps', type=float, default=None, required=False,
                        help='number of ps to run equil (not used for gb/pbsa calc)')
    parser.add_argument('--calcFrames', type=int, default=None, required=False,
                        help='number of frames averaged over trajectory for calculation')
    parser.add_argument('--mbar', type=int, required=False, default=None, help='use pymbar to subsample from this number of frames (should be > 25 at least)')
    parser.add_argument('--method', type=str, choices=['gbsa', 'pbsa'], default='gbsa', help='use pbsa or gbsa')

    ## logging options
    parser.add_argument('--odir', type=str, default=None, required=False, help='directory for intermediate files '
                                                                               'and outputs')
    parser.add_argument('-v', type=int, choices=[0, 1, 2], default=1, required=False, help='verbose level')

    return parser.parse_args()

def setup_folder(args):
    if args.odir is None:
        args.odir = f'{os.getcwd()}/{args.pdb.split("/")[-1].split(".")[0]}'

    try:
        os.mkdir(args.odir)
    except OSError:
        if not os.path.isdir(args.odir):
            print("Creation of the directory %s failed" % args.odir)
            exit()
        else:
            shutil.rmtree(args.odir)
            os.mkdir(args.odir)
    else:
        if args.v >= 1:
            print("Successfully created the directory %s " % args.odir)

if __name__ == '__main__':
    args = get_args()

    setup_folder(args)

    if args.amber_path is None:
        args.amber_path = os.environ['AMBERHOME']

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.ERROR)
    logging.getLogger('simtk').setLevel(logging.ERROR)
    logging.getLogger('openforcefield').setLevel(logging.ERROR)
    logging.getLogger('openmmtools').setLevel(logging.ERROR)
    warnings.filterwarnings("ignore")

    sl = SystemLoader(dirpath=args.odir,
                      verbose=args.v,
                      input_pdb=args.pdb,
                      ps=args.ps,
                      calcFrames=args.calcFrames,
                      equil_ps=args.equil_ps,
                      platform_name=args.platform,
                      mbar=args.mbar
                      )
    sl.prepare_simulation()

    deltag, std = sl.run_amber(args.method, args.amber_path)
    print(f"{deltag} ± {std} (kcal/mol)")
