# DFT adsorption dataset for W-B nanoclusters

This directory contains the processed DFT data used for the W-B nanocluster
adsorption analysis. Energies are given in eV, distances in angstrom, and angles
in degrees.

## Files

### `WB_nanoclusters_bare.csv`

Reference total energies for 21 optimized bare W-B nanoclusters.

- `cluster`: W-B nanocluster composition.
- `E_total_eV`: total energy of the optimized bare cluster.

### `WB_nanoclusters_NH3_adsorption.csv`

Detailed data for 113 optimized molecular-NH3 configurations on W sites. This
table retains total energies and intramolecular geometric descriptors.

- `cluster`: W-B nanocluster composition.
- `site_id`: zero-based adsorption-configuration identifier within a cluster.
- `E_total_eV`: total energy of the optimized NH3-adsorbed system.
- `E_ads_NH3_eV`: molecular-NH3 adsorption energy.
- `d_N_W_min_A`: minimum N-W distance.
- `d_N_H_mean_A`: mean N-H bond length.
- `angle_HNH_mean_deg`: mean H-N-H angle.

### `WB_nanoclusters_NH3_adsorption_W_B_sites.csv`

Canonical consolidated molecular-NH3 adsorption table used for the W/B site
comparison. It contains 286 configurations: 113 W-site configurations and 173
accessible B-site configurations. The B-site geometric descriptors were
calculated from the final optimized `CONTCAR` structures.

- `cluster`: W-B nanocluster composition.
- `site_label`: element-local adsorption-site label, for example `W1` or `B1`.
- `anchor_element`: adsorption-site element, `W` or `B`.
- `E_total_eV`: total energy of the optimized NH3-adsorbed system.
- `E_ads_NH3_eV`: molecular-NH3 adsorption energy.
- `d_NH3_host_A`: optimized distance between N and the host atom.
- `d_N_H_mean_A`: mean of the three N-H bond lengths.
- `angle_HNH_mean_deg`: mean of the three H-N-H angles.

### `WB_nanoclusters_H_adsorption.csv`

Current table of 304 optimized atomic-H adsorption configurations.

- `cluster`: W-B nanocluster composition.
- `site_label`: element-local adsorption-site label, for example `W1` or `B1`.
- `anchor_element`: nearest host element after optimization, `W` or `B`.
- `E_ads_H_eV`: atomic-H adsorption energy.
- `d_H_host_A`: optimized distance between H and the host atom.

## Energy definitions

The molecular-NH3 adsorption energy is

`E_ads(NH3) = E(cluster + NH3) - E(cluster) - E(NH3)`.

The atomic-H adsorption energy is referenced to molecular hydrogen,

`E_ads(H) = E(cluster + H) - E(cluster) - 1/2 E(H2)`.

Negative adsorption energies correspond to exothermic adsorption.
