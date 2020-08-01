import subprocess
from sys import stdout

import mdtraj as md
import numpy as np
import simtk.openmm as mm
from openeye import oechem, oequacpac
from pymbar import timeseries
from simtk import unit
from simtk.openmm import app

from mmgpbsa.amber_mmgpbsa import run_amber
from mmgpbsa.utils import make_message_writer, working_directory


def subsample(enthalpies):
    """
    Subsamples the enthalpies using John Chodera's code.
    This is probably better than the simple cutoff we normally use.
    No output -- it modifies the lists directly
    """
    # Use automatic equilibration detection and pymbar.timeseries to subsample
    [t0, g, Neff_max] = timeseries.detectEquilibration(enthalpies)
    enthalpies = enthalpies[t0:]
    return timeseries.subsampleCorrelatedData(enthalpies, g=g)


class Config:
    defaults_config = {

        ## System Params
        'nonbondedMethod': app.CutoffNonPeriodic,
        'rigid_water': True,
        'removeCMMotion': True,

        "nonbondedCutoff": 1.0 * unit.nanometer,
        "constraints": app.HBonds,
        "implicitSolvent": app.GBn2,
        "soluteDielectric": 1.0,
        "solventDielectric": 78.5,
        "igb": 5,

        ## Integrator Params
        "temperature": 310.15 * unit.kelvin,
        "frictionCoeff": 1.0 / unit.picoseconds,
        "stepSize": 2.0 * unit.femtoseconds,
        "constraintTolerance": 0.00001,

        'platform_name': 'CPU',
        'platform_properties': {},

        'reportInterval': 1250,
        'equil_ps': 10,
        'ps': 10,
        'calcFrames': None,
        'mbar': 50
    }

    def __init__(self, **kwargs):

        self.__dict__.update(self.defaults_config)

        for k, v in kwargs.items():
            if k in self.__dict__:
                if v is not None:
                    self.__dict__[k] = v

        ## Simulation Params
        if self.platform_name in ['OpenCL', 'CUDA']:
            self.platform_properties = {'Precision': 'mixed'}

        if not isinstance(self.ps, unit.Quantity):
            self.ps = self.ps * unit.picosecond
        if not isinstance(self.equil_ps, unit.Quantity):
            self.equil_ps = self.equil_ps * unit.picosecond

        if self.mbar is not None and self.mbar > 0:
            self.calcFrames = self.mbar

        self.total_ps = self.equil_ps + self.ps
        self.equil_steps = int(self.equil_ps / self.stepSize)
        self.production_steps = int(self.ps / self.stepSize)
        self.total_steps = int(self.total_ps / self.stepSize)
        self.trajInterval = int(self.production_steps / self.calcFrames)


class SystemLoader:

    def __init__(self, dirpath, verbose, input_pdb, **kwargs):
        self.verbose = verbose
        self.logger = make_message_writer(self.verbose >= 1, self.__class__.__name__, verbose >= 2)

        with self.logger("__init__") as logger:
            self.dirpath = dirpath
            self.input_pdb = input_pdb
            self.config = Config(**kwargs)

    def split_complex_from_system(self):
        with self.logger("split_complex_from_system") as logger:
            pdb = oechem.OEMol()
            prot = oechem.OEMol()
            lig = oechem.OEMol()
            wat = oechem.OEGraphMol()
            other = oechem.OEGraphMol()
            ifs = oechem.oemolistream()
            ifs.SetFlavor(oechem.OEFormat_PDB,
                          oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_DATA | oechem.OEIFlavor_PDB_ALTLOC)  # noqa
            if ifs.open(self.input_pdb):
                oechem.OEReadMolecule(ifs, pdb)
            else:
                logger.error("Could not open file!")
            ifs.close()
            logger.log(f"Reading input PDB {self.input_pdb}")
            if not oechem.OESplitMolComplex(lig, prot, wat, other, pdb):
                logger.failure("could not split complex. exiting", exit_all=True)
            else:
                logger.log(
                    f"Split complex. atom sizes-- lig: {len(list(lig.GetAtoms()))}, prot: {len(list(prot.GetAtoms()))}, water: {len(list(wat.GetAtoms()))}, other: {len(list(other.GetAtoms()))}")

            return prot, lig

    def prepare_ligand(self, lig):
        with self.logger("prepare_ligand") as logger:
            ## prepare ligand
            oemol = oechem.OEMol(lig)
            oemol.SetTitle("UNL")
            oechem.OEAddExplicitHydrogens(oemol)

            ofs = oechem.oemolostream(f'{self.dirpath}/lig.mol2')
            oechem.OEWriteMolecule(ofs, oemol)
            ofs.close()

            ofs = oechem.oemolostream(f'{self.dirpath}/lig.pdb')
            oechem.OEWriteMolecule(ofs, oemol)
            ofs.close()

            self.lig = oechem.OEMol(oemol)

            oequacpac.OEAssignCharges(oemol, oequacpac.OEAM1BCCCharges())

            ofs = oechem.oemolostream(f'{self.dirpath}/charged.mol2')
            oechem.OEWriteMolecule(ofs, oemol)
            ofs.close()

            self.charged_lig = oechem.OEMol(oemol)

    def prepare_protein(self, protein):
        with self.logger("prepare_protein") as logger:
            ofs = oechem.oemolostream()
            oemol = oechem.OEMol(protein)

            if ofs.open(f'{self.dirpath}/apo.pdb'):
                oechem.OEWriteMolecule(ofs, oemol)
            ofs.close()

            self.apo = oemol

    def __setup_system_im(self):
        with self.logger("__setup_system_im") as logger:
            try:
                protein, lig = self.split_complex_from_system()
                self.prepare_ligand(lig)
                self.prepare_protein(protein)

                with working_directory(self.dirpath):
                    subprocess.run(
                        f'antechamber -i lig.mol2 -fi mol2 -o lig.mol2 -fo mol2 -pf y -an y -a charged.mol2 -fa mol2 -ao crg'.split(
                            " "), check=True, capture_output=True)
                    subprocess.run(f'parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod'.split(" "), check=True,
                                   capture_output=True)
                    try:
                        subprocess.run('pdb4amber -i apo.pdb -o apo_new.pdb --reduce --dry'.split(" "), check=True,
                                       capture_output=True)
                    except subprocess.CalledProcessError as e:
                        logger.error("Known bug, pdb4amber returns error when there was no error", e.stdout,
                                     e.stderr)
                        pass

                    # Wrap tleap
                    with open('leap.in', 'w+') as leap:
                        leap.write("source leaprc.protein.ff14SBonlysc\n")
                        # leap.write("source leaprc.phosaa10\n")
                        leap.write("source leaprc.gaff2\n")
                        leap.write("set default PBRadii mbondi3\n")
                        leap.write("rec = loadPDB apo_new.pdb\n")  # May need full filepath?
                        leap.write("saveAmberParm rec apo.prmtop apo.inpcrd\n")
                        leap.write("lig = loadmol2 lig.mol2\n")
                        leap.write("loadAmberParams lig.frcmod\n")
                        leap.write("saveAmberParm lig lig.prmtop lig.inpcrd\n")
                        leap.write("com = combine {rec lig}\n")
                        leap.write("saveAmberParm com com.prmtop com.inpcrd\n")
                        leap.write("quit\n")
                    try:
                        subprocess.run('tleap -f leap.in'.split(" "), check=True, capture_output=True)
                    except subprocess.CalledProcessError as e:
                        logger.error("tleap error", e.output.decode("UTF-8"))
                        exit()

                    prmtop = app.AmberPrmtopFile(f'com.prmtop')
                    inpcrd = app.AmberInpcrdFile(f'com.inpcrd')
                    with open('com.pdb', 'w') as f:
                        app.pdbfile.PDBFile.writeFile(prmtop.topology, inpcrd.positions, file=f)

                system = prmtop.createSystem(nonbondedMethod=self.config.nonbondedMethod,
                                             nonbondedCutoff=self.config.nonbondedCutoff,
                                             rigidWater=self.config.rigid_water,
                                             constraints=self.config.constraints,
                                             implicitSolvent=self.config.implicitSolvent,
                                             soluteDielectric=self.config.soluteDielectric,
                                             solventDielectric=self.config.solventDielectric,
                                             removeCMMotion=self.config.removeCMMotion)
                self.topology, self.positions = prmtop.topology, inpcrd.positions
                return system
            except Exception as e:
                logger.error("EXCEPTION CAUGHT BAD SPOT", e)

    def run_amber(self, method: str, amber_path):
        with self.logger('run_amber') as logger:
            logger.log("Calculating mmgb/pbsa value...may take awhile.")
            return run_amber(amber_path, self.dirpath, verbose=self.verbose, igb=self.config.igb,
                             use_pbsa=(method == 'pbsa'))

    def prepare_simulation(self):
        with self.logger("prepare_simulation") as logger:

            system = self.__setup_system_im()
            integrator = mm.LangevinIntegrator(self.config.temperature, self.config.frictionCoeff,
                                               self.config.stepSize)
            integrator.setConstraintTolerance(self.config.constraintTolerance)

            platform = mm.Platform.getPlatformByName(self.config.platform_name)
            simulation = app.Simulation(self.topology, system, integrator, platform, self.config.platform_properties)
            logger.log(
                f"Built simulation using platform {self.config.platform_name} with properties {self.config.platform_properties}")
            simulation.context.setPositions(self.positions)

            logger.log(f'Minimizing and setting velocities to {self.config.temperature}')
            simulation.minimizeEnergy()
            simulation.context.setVelocitiesToTemperature(self.config.temperature)

            if self.verbose >= 1:
                simulation.reporters.append(
                    app.StateDataReporter(stdout, max(self.config.reportInterval - 1, 1), step=True,
                                          potentialEnergy=True, temperature=True,
                                          progress=True,
                                          remainingTime=True,
                                          speed=True, totalSteps=self.config.total_steps,
                                          separator='\t'))

            logger.log(f'Equilibrating for {self.config.equil_steps} steps, or {self.config.equil_ps}')
            simulation.step(self.config.equil_steps)

            simulation.reporters.append(
                app.DCDReporter(f'{self.dirpath}/trajectory.dcd',
                                max(self.config.trajInterval - 1 if self.config.mbar == 0 else self.config.mbar, 1)))

            logger.log(f'Running Production for {self.config.production_steps} steps, or {self.config.ps}')
            if self.config.mbar > 0:
                enthalpies = np.zeros((self.config.mbar))
                for i in range(self.config.mbar):
                    simulation.step(int(self.config.production_steps / self.config.mbar))
                    state = simulation.context.getState(getEnergy=True)
                    potential_energy = state.getPotentialEnergy()
                    enthalpies[i] = potential_energy.value_in_unit(unit.kilojoules_per_mole)
                logger.log(f"Simulation Done. Running MBAR on {self.config.mbar} snapshots.")

                idx = np.array(subsample(enthalpies))
                traj = md.load(f'{self.dirpath}/trajectory.dcd', top=f'{self.dirpath}/com.pdb')
                traj = traj[idx]
                traj.save_dcd(f'{self.dirpath}/trajectory.dcd', force_overwrite=True)

                logger.log(f"Done. Subsampled {len(idx)} from {self.config.mbar} snapshots.")

                return enthalpies
            else:
                simulation.step(self.config.production_steps)
                logger.log(f"Simulation Done")
