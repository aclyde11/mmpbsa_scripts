import argparse

from openeye import oechem, oedocking


def split_complex_from_system(complex_file):
    pdb = oechem.OEMol()
    prot = oechem.OEMol()
    lig = oechem.OEMol()
    wat = oechem.OEGraphMol()
    other = oechem.OEGraphMol()
    ifs = oechem.oemolistream()
    ifs.SetFlavor(oechem.OEFormat_PDB,
                  oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_DATA | oechem.OEIFlavor_PDB_ALTLOC)  # noqa
    if ifs.open(complex_file):
        oechem.OEReadMolecule(ifs, pdb)
    else:
        print("Could not open file!")
        exit()
    ifs.close()
    if not oechem.OESplitMolComplex(lig, prot, wat, other, pdb):
        print("could not split complex. exiting")
        exit()

    return prot, lig


def get_receptor(complex_file):
    prot, lig = split_complex_from_system(complex_file)
    receptor = oechem.OEGraphMol()
    oedocking.OEMakeReceptor(receptor, prot, lig)
    return receptor, lig


def score_pose(receptor, ligand):
    dock = oedocking.OEDock(oedocking.OEDockMethod_Chemgauss4, oedocking.OESearchResolution_Default)
    dock.Initialize(receptor)

    score = dock.ScoreLigand(ligand)

    comps = {'Chemgauss4 Score': score}
    for component in dock.GetComponentNames():
        comps[component] = dock.ScoreLigandComponent(ligand, component)

    return comps


def score_and_build(complex_pdb):
    receptor, lig = get_receptor(complex_pdb)
    scores = score_pose(receptor, lig)
    return scores


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--complex', type=str, default=None, required=False)
    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()
    scores = score_and_build(args.complex)
    print(scores)
