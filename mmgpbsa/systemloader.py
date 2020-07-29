import subprocess
from sys import stdout

import simtk.openmm as mm
from openeye import oechem, oequacpac
from simtk import unit
from simtk.openmm import app

from mmgpbsa.amber_mmgpbsa import run_amber
from mmgpbsa.utils import make_message_writer, working_directory


class Config:
    def __init__(self, calcFrames = None, temperature=None, frictionCoeff=None, stepSize=None, constraintTolerance=None,
                 use_gpu=None, reportInterval=None, ps=None, nonbondedCutoff=None, constraints=None,
                 implicitSolvent=None, equil_steps=None, solventDielectric=None, soluteDielectric=None):


        ## System Params
        self.nonbondedMethod = app.CutoffNonPeriodic
        self.rigid_water = True
        self.removeCMMotion = True

        self.nonbondedCutoff = nonbondedCutoff or (1.0 * unit.nanometer)
        self.constraints = constraints or app.HBonds
        self.implicitSolvent = implicitSolvent or app.GBn2
        self.soluteDielectric = soluteDielectric or 1.0
        self.solventDielectric = solventDielectric or 78.5

        ## Integrator Params
        self.temperature = temperature or (300 * unit.kelvin)
        self.frictionCoeff = frictionCoeff or (1.0 / unit.picoseconds)
        self.stepSize = stepSize or (2.0 * unit.femtoseconds)
        self.constraintTolerance = constraintTolerance or 0.00001

        ## Simulation Params
        if use_gpu:
            self.platform_name = 'CUDA'
            self.platform_properties = {'Precision': 'mixed'}
        else:
            self.platform_name = 'OpenCL'
            self.platform_properties = {'Precision': 'mixed'}

            # self.platform_name = 'CPU'
            # self.platform_properties = {}

        ## Log Parms
        self.reportInterval = reportInterval or 5000 # 10ps if stepsize is 2fs
        self.equil_steps = equil_steps or 100
        self.ps = ps or 20 #
        self.ps = self.ps * unit.picosecond
        self.step_size = self.stepSize
        self.total_steps = int(self.ps / self.step_size)

        self.calcFrames = calcFrames or 10
        self.trajInterval = int(self.total_steps / self.calcFrames)


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

    def prepare_protein(self, protein, verbose=True):
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

                system = prmtop.createSystem(nonbondedMethod=self.config.nonbondedMethod,
                                             nonbondedCutoff=self.config.nonbondedCutoff, rigidWater=self.config.rigid_water,
                                             constraints=self.config.constraints, implicitSolvent=self.config.implicitSolvent,
                                             soluteDielectric=self.config.soluteDielectric, solventDielectric=self.config.solventDielectric, removeCMMotion=self.config.removeCMMotion)
                self.topology, self.positions = prmtop.topology, inpcrd.positions
                return system
            except Exception as e:
                logger.error("EXCEPTION CAUGHT BAD SPOT", e)


    def run_amber(self, method, amber_path):
        with self.logger('run_amber') as logger:
            logger.log("Calculating mmgb/pbsa value...may take awhile.")
            return run_amber(amber_path, self.dirpath, verbose=self.verbose, pbsa=(method == 'pbsa'))

    def prepare_simulation(self):
        with self.logger("prepare_simulation") as logger:

            system = self.__setup_system_im()
            integrator = mm.LangevinIntegrator(self.config.temperature, self.config.frictionCoeff,
                                               self.config.stepSize)
            integrator.setConstraintTolerance(self.config.constraintTolerance)

            platform = mm.Platform.getPlatformByName(self.config.platform_name)
            simulation = app.Simulation(self.topology, system, integrator, platform, self.config.platform_properties)
            logger.log(f"Built simulation using platform {self.config.platform_name} with properties {self.config.platform_properties}")
            simulation.context.setPositions(self.positions)

            logger.log(f'Minimizing and setting velocities to {self.config.temperature}')
            simulation.minimizeEnergy()
            simulation.context.setVelocitiesToTemperature(self.config.temperature)
            logger.log(f'Equilibrating for {self.config.equil_steps * self.config.step_size} steps, or {(self.config.equil_steps * self.config.step_size).in_units_of(unit.picosecond)}')
            simulation.step(self.config.equil_steps)

            simulation.reporters.append(
                    app.DCDReporter(f'{self.dirpath}/trajectory.dcd', self.config.trajInterval))

            if self.verbose >= 1:
                simulation.reporters.append(app.StateDataReporter(stdout, self.config.reportInterval, step=True,
                                                                  potentialEnergy=True, temperature=True,
                                                                  progress=True,
                                                                  remainingTime=True,
                                                                  speed=True, totalSteps=self.config.total_steps + self.config.equil_steps,
                                                                  separator='\t'))

            logger.log(f'Running Production for {self.config.total_steps} steps, or {(self.config.total_steps * self.config.step_size).in_units_of(unit.nanosecond)}')
            simulation.step(self.config.total_steps)
            logger.log(f"")