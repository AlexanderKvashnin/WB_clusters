# CI-NEB reaction-path dataset

This directory contains the reported climbing-image nudged elastic band
(CI-NEB) calculations for the first N-H bond cleavage of molecular ammonia on
isolated W-B nanoclusters,

`NH3* -> NH2* + H*`.

## Scope of the sampling

The dataset is a targeted sampling of local, chemically plausible first-cleavage
pathways. It is not an exhaustive enumeration of all ammonia adsorption states,
all three N-H bonds, all possible H-accepting atoms, reconstructed products, or
subsequent H-migration pathways.

Initial states (IS) were selected from relaxed, molecular, W-bound NH3
configurations. The selection was intended to represent distinct local W
environments and configurations for which an N-H bond was oriented toward a
nearby W- or B-containing environment. Candidate final states (FS) were then
constructed by transferring one H atom to a nearby local acceptor chosen by
geometric inspection. When more than one H orientation or acceptor environment
was judged plausible, multiple pathways were calculated from the same IS.
The resulting NH2*+H* endpoints were relaxed before the NEB band was built.
No universal distance cutoff or exhaustive nearest-neighbour enumeration was
used; consequently, the number of sampled pathways differs among clusters.

Every pathway connects separately relaxed IS and FS structures through eight
intermediate images, giving ten archived images in total. The highest-energy
image in each archived profile is reported as the CI-NEB transition-state
estimate (TS). The reported pathway means and standard deviations are
descriptive statistics of this explicitly sampled set, not exhaustive
composition-specific kinetic averages.

## Files and identifiers

- `neb_trajectories_xyz/` contains one ten-frame XYZ file per pathway.
- A filename of the form
  `W_MM_B_NN_<initial-state>_<pathway>_movie_contcar.xyz` records the cluster
  composition, the relaxed molecular-NH3 initial-state identifier, and the
  pathway identifier within that initial state.
- Each XYZ comment line stores the electronic energy relative to image 0.
- `NEB_wb.csv` is the canonical channel-level manifest. The absolute IS energy
  is retained as the reference from the original calculation table; FS and TS
  energies, barriers, reaction energies, and endpoint labels are regenerated
  from the archived trajectory associated with that exact pathway identifier.
- `NEB_sampling_by_cluster.csv` contains the corresponding cluster-level counts
  and descriptive statistics.
- `build_neb_manifest.py` regenerates both CSV files and validates the pathway
  identifiers and the energy identities.

## Endpoint labels in the manifest

Atom identities are derived without manual relabelling:

- `Transferred_H` is the H atom with the largest endpoint displacement after
  subtracting the median displacement of the cluster atoms.
- `Initial_NH3_W_anchor` is the W atom closest to N in image 0.
- `Final_H_anchor` is the W or B atom closest to the transferred H in image 9.
- Element-local indices are one-based; for example, `B9` is the ninth B atom in
  the trajectory atom order.
- `Final_H_anchor_rank_initial_excluding_binding_W` ranks that final anchor by
  its image-0 distance from the transferred H after excluding the W atom that
  binds NH3. A value of zero means that H ends on the NH3-binding W itself.

The endpoint labels describe the relaxed archived structures. They do not imply
that the final nearest atom was necessarily the only atom considered during the
manual construction of the unrelaxed product guess, because endpoint relaxation
can reconstruct the local environment.

## Energy definitions

For every row,

`Delta_E = E_FS - E_IS`

`E_act_fwd = E_TS - E_IS`

`E_act_rev = E_TS - E_FS`.

The build script derives these quantities from a single trajectory energy
profile and checks the three identities before writing the manifest. Energies
are electronic energies in eV; zero-point, entropic, and finite-temperature
corrections are not included.
