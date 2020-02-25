# GNB2_analysis
In silico analysis of the variants for the GNB2 paper

## Templates

In the folder templates are 6 models.

All use the [Swissmodel GNB2 model](https://swissmodel.expasy.org/repository/5e5058728fd6f9e51e1ef5bc.pdb), 
which is based on the 90% identical GNB1 structure (template: chain B of PDB:3CIK)
and are energy minimised with Rosetta using 15 cycles of FastRelax.

See [templates/template_generation.py](templates/template_generation.py) for code.

The models are:

* **unpaired**: β-subunit (GNB2) and γ-subunit (GNG2) of G-protein
* **with Alpha**: β-subunit (GNB2) and γ-subunit (GNG2) with the α-subunit (GNAI1) based upon PDB:6CRK
* **with GRK2**: β-subunit (GNB2) and γ-subunit (GNG2) with β-adrenergic receptor kinase 1 (GRK2) based upon PDB:3CIK
* **with KCTD12**: complex of multiple β- and γ-subunit with KCTD12 based upon PDB:6M8S
* **with PREX1**: β-subunit (GNB2) and γ-subunit (GNG2)  with PREX1 based upon PDB:6PCV (latter left bovine)
* **phosphorylated**: post-translational modifications taken from Phosphosite Plus

RASD2 (Rhes) has strong homology to solved homologues, except in the GNB2 binding region (C-terminus),
consequently a model could not be made. Literature suggests it is were the alpha helix of the α-subunit.
See [RASD2_notes](RASD2/RASD2_notes.md).
