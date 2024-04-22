import numpy as np
from openff.qcsubmit.results import (
    BasicResultCollection,
    OptimizationResultCollection,
)
from openff.qcsubmit.results.filters import LowestEnergyFilter
from openff.toolkit import ForceField, Molecule
from qubekit.utils import constants
from tqdm import tqdm


def write_hess():
    """grab lowest-energy conformers and keep only those with hessians, storing
    the resulting `BasicResultCollection` to file"""
    opt = (
        OptimizationResultCollection.parse_file(
            "../02_curate-data/datasets/combined-opt.json"
        )
        .filter(LowestEnergyFilter())
        .to_basic_result_collection(driver="hessian")
    )

    with open("hess.json", "w") as f:
        f.write(opt.json(indent=2))


def get_hess():
    "Write the first molecule and its hessian to files for testing"
    ds = BasicResultCollection.parse_file("hess.json")
    for rec, mol in tqdm(ds.to_records(), desc="Processing records"):
        with open("test.mol", "w") as out:
            print(mol.to_json(indent=2), file=out)
        np.savetxt("test.hess", rec.return_result)
        break


def load_mol(filename):
    with open(filename) as inp:
        return Molecule.from_json(inp.read())


def label_molecule(ff: ForceField, mol: Molecule, kind="ProperTorsions"):
    return ff.label_molecules(mol.to_topology())[0][kind]


def unit(v):
    return v / np.linalg.norm(v)


def mag(v):
    return np.linalg.norm(v)


def mag2(v):
    t = np.linalg.norm(v)
    return t * t


def partial_hess(hess, ai, bi):
    """returns the eigenvalues and eigenvectors of the partial hessian for
    atomic indices `ai` and `bi`"""
    h = hess[(ai * 3) : ((ai + 1) * 3), (bi * 3) : ((bi + 1) * 3)]
    return np.linalg.eig(h)


def msm_torsion(hess, mol, indices):
    # based on the qubekit code, the hessian is in atomic units, convert to
    # kcal/mol/A^2
    conversion = constants.HA_TO_KCAL_P_MOL / (constants.BOHR_TO_ANGS**2)
    _hess = hess * conversion

    ai, bi, ci, di = indices

    assert len(mol.conformers) == 1

    geom = mol.conformers[0]

    # these are quantities in angstroms
    a, b, c, d = (
        geom[0].magnitude,
        geom[1].magnitude,
        geom[4].magnitude,
        geom[3].magnitude,
    )

    u_ab = unit(b - a)
    u_bc = unit(c - b)
    u_cb = -u_bc
    u_cd = unit(d - c)
    u_dc = -u_cd

    r_ba = mag(a - b)
    r_cd = mag(d - c)

    u_n_abc = unit(np.cross(u_cb, u_ab))
    u_n_bcd = unit(np.cross(u_dc, u_bc))

    lam_ab, nu_ab = partial_hess(_hess, ai, bi)
    lam_dc, nu_dc = partial_hess(_hess, di, ci)

    # seminario paper mentions "including the minus sign" so this is probably
    # off by a sign at least. qubekit also multiplies bonds by -1/2, which
    # aligns (kinda) with the msm paper saying angles are overcounted by a
    # factor of 2
    inv_k = 1.0 / (
        r_ba
        * r_ba
        * mag2(np.cross(u_ab, u_bc))
        * sum(lam_ab[i] * mag(np.dot(u_n_abc, nu_ab[i])) for i in range(3))
    ) + 1.0 / (
        r_cd
        * r_cd
        * mag2(np.cross(u_bc, u_cd))
        * sum(lam_dc[i] * mag(np.dot(u_n_bcd, nu_dc[i])) for i in range(3))
    )

    return 1.0 / inv_k


ff = ForceField("openff-2.1.0.offxml")


# first torsion is (0, 1, 4, 3) and assigned t44 with k1 = 3.26 kcal/mol
hess = np.loadtxt("test.hess")
mol = load_mol("test.mol")

for env, p in label_molecule(ff, mol).items():
    k = msm_torsion(hess, mol, env)
    print(env, p.id, [k.magnitude for k in p.k], k)
