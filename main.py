import argparse
import logging
import os
import warnings

from amber_mmgpbsa import run_amber
from systemloader import SystemLoader

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--amber_path', type=str, default=None, required=True, help='path to amber installation folder (ex /home/austin/amber20)')
    parser.add_argument('--pdb', type=str, help='pdb input file with ligand and protein in complex. Only one chain for now please!')
    parser.add_argument('--use_gpu', action='store_true')
    ##simulation options
    parser.add_argument('--ps', type=int, default=None, required=False, help='picoseconds to run simulation')
    parser.add_argument('--equil_steps', type=int, default=None, required=False, help='number of steps to run equil (not used for gb/pbsa calc)')
    parser.add_argument('--traj_frames', type=int, default=None, required=False, help='number of frames averaged over trajectory for calculation')
    parser.add_argument('--method', type=str, choices=['gbsa', 'pbsa'], default='gbsa', help='use pbsa or gbsa')

    ## logging options
    parser.add_argument('--odir', type=str, default=None, required=False, help='directory for intermediate files '
                                                                                 'and outputs')
    parser.add_argument('-v', type=int, choices=[0, 1, 2], default=1, required=False, help='verbose level')

    return parser.parse_args()

if __name__ == '__main__':
    args = get_args()
    args.v = 0


    if args.odir is None:
        args.odir = f'{os.getcwd()}/{args.pdb.split("/")[-1].split(".")[0]}'

    try:
        os.mkdir(args.odir)
    except OSError:
        if not os.path.isdir(args.odir):
            print("Creation of the directory %s failed" % args.odir)
            exit()
    else:
        if args.v >= 1:
            print("Successfully created the directory %s " % args.odir)

    if args.method == 'pbsa':
        print("not ready yet.")
        exit()


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
                      calcFrames=args.traj_frames,
                      equil_steps=args.equil_steps,
                      use_gpu=args.use_gpu
                      )
    sl.prepare_simulation()

    print(run_amber(args.amber_path, args.odir, verbose=args.v, pbsa=(args.method == 'pbsa')))


