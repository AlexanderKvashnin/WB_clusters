# DFT Dataset for W-B Nanoclusters

This folder contains processed DFT data for W-B nanoclusters. The dataset includes reference energies of bare nanoclusters, NH3 adsorption data, and atomic H adsorption data.

All energies are given in eV, distances in Å, and angles in degrees.

## Files

### `WB_nanoclusters_bare.csv`

This table contains total DFT energies of optimized bare W-B nanoclusters. These energies are used as reference values for adsorption-energy calculations.

Columns:

- `cluster`: W-B nanocluster composition.
- `E_total_eV`: total energy of the optimized bare nanocluster.

---

### `WB_nanoclusters_NH3_adsorption_clean.csv`

This table contains optimized NH3 adsorption configurations on W-B nanoclusters.

Columns:

- `cluster`: W-B nanocluster composition.
- `site_id`: NH3 adsorption-site index. Numbering starts from `0` for each cluster.
- `E_total_eV`: total energy of the optimized NH3-adsorbed system.
- `E_ads_NH3_eV`: NH3 adsorption energy.
- `d_N_W_min_A`: minimum distance between the N atom of NH3 and the nearest W atom.
- `d_N_H_mean_A`: mean N-H bond length in the adsorbed NH3 molecule.
- `angle_HNH_mean_deg`: mean H-N-H angle in the adsorbed NH3 molecule.

The NH3 adsorption energy was calculated as:
E_ads(NH3) = E(cluster + NH3) - E(cluster) - E(NH3)
More negative values correspond to stronger NH3 adsorption.
---

### `WB_nanoclusters_H_adsorption.csv`

This table contains atomic H adsorption energies and selected geometric descriptors for H adsorption on W-B nanoclusters.

Repeated or nearly equivalent initial H adsorption sites were grouped into representative adsorption states.

Columns:

- `cluster`: W-B nanocluster composition.
- `site_label`: processed adsorption-site label, for example `B1`, `B2`, or `W1`.
- `anchor_element`: host atom associated with H adsorption, either `B` or `W`.
- `E_ads_H_eV`: atomic H adsorption energy.
- `d_H_host_A`: distance between the adsorbed H atom and the host atom.
- `E_total_H_adsorbed_eV`: total energy of the optimized H-adsorbed system.
- `representative_source_site`: original source site selected as the representative structure.
- `n_equivalent_sites`: number of original sites grouped into this representative state.
- `equivalent_source_sites`: original source sites grouped together, separated by semicolons.

The atomic H adsorption energy was calculated as:
E_ads(H) = E(cluster + H) - E(cluster) - E(H)
More negative values correspond to stronger H adsorption.

## Notes

The file `WB_nanoclusters_bare.csv` provides the reference energies for adsorption-energy calculations.
For NH3 adsorption, site numbering starts from `0` for each cluster.
For H adsorption, labels such as `B1` and `W1` indicate whether H is associated with a B or W host atom.