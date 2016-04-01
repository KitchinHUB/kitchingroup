from ase.db import connect
# from ase.visualize import view
import numpy as np

db = connect('example-db1.db')

ids, images = [], []
# A search performed only on key 'Pd!=0' to exclude pure oxygen.
for d in db.select(['Pd!=0']):
    # Calculate the difference between NN predicted energies.

    # I do not store these because each key can dramatically increase
    # the size of a database.
    delta = abs(d.data['nn1'] - d.data['nn2']) / d.natoms

    if delta > 0.2:
        atoms = db.get_atoms(d.id)
        # Uncommenting this will produce 30 instances of ase-gui
        # view(atoms)
        images += [atoms]
        ids += [d.id]

print('{} poorly fit calculations'.format(len(images)))

# Not all poorly fit calculations are necessarily bad.
# This loop determines the distances between oxygen atoms in each 
# poorly fit structure.
bad_ids = []
for i, atoms in enumerate(images):
    n = len(atoms)
    dis = []

    for j, atom in enumerate(atoms):
        if atom.symbol == 'O':
            dis += [atoms.get_distances(j, range(n)[36:])]

    dis = np.matrix(dis)
    off_diag = np.extract(1 -  np.eye(len(dis)), dis)
    if (off_diag < 0.9).any():
        bad_ids += [ids[i]]
        # view(atoms)

print('{} infeasible configurations'.format(len(bad_ids)))

# This will delete the bad_ids from the next iteration of the database.
# db0 = connect('networks/db2/data.db')
# db0.delete(bad_ids)
