import subprocess

from mmgpbsa.utils import working_directory, parse_final


def get_mmgbsa_input(use_pbsa=False, igb=5):
    if use_pbsa:
        str = (f"""&general
startframe=1,
verbose=2,
keep_files=0,
/
&pb
        istrng=0.15, fillratio=4.0
/
""")
    else:
        str = (f"""&general
startframe=1,
verbose=2,
keep_files=0,
/
&gb
        igb={igb}, saltcon=0.150,
/
""")
    with open("mmpbsa.in", 'w') as f:
        f.write(str)


def run_amber(amber_path, dir_path, use_pbsa=False, igb=5, verbose=1):
    with working_directory(dir_path):
        get_mmgbsa_input(use_pbsa, igb)

        args = {}
        if verbose < 2:
            args = {'stdout': subprocess.DEVNULL, 'stderr': subprocess.DEVNULL}

        with open("runamber.sh", 'w') as f:
            f.write(f"source {amber_path}/amber.sh\n")
            f.write(f"MMPBSA.py -O -i mmpbsa.in -cp com.prmtop -rp apo.prmtop -lp lig.prmtop -y trajectory.dcd\n")
        subprocess.run(['bash', 'runamber.sh'], **args)

        if verbose >= 2:
            with open("FINAL_RESULTS_MMPBSA.dat", 'r') as f:
                print(f.read())
        return parse_final('FINAL_RESULTS_MMPBSA.dat')
