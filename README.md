# GNB2_analysis
In silico analysis of the variants for the GNB2 paper

## Templates

In the folder templates are 6 models.

All use the [Swissmodel GNB2 model](https://swissmodel.expasy.org/repository/5e5058728fd6f9e51e1ef5bc.pdb), which is based on the 90% identical GNB1 structure (template: chain B of PDB:3CIK)
and are energy minimised with Rosetta using 15 cycles of FastRelax.

See [templates/template_generation.py](templates/template_generation.py) for code.

The models are:

* **unpaired**: β-subunit (GNB2) and γ-subunit (GNG2) of G-protein
* **with Alpha**: β-subunit (GNB2) and γ-subunit (GNG2) with the α-subunit (GNAI1) based upon PDB:6CRK
* **with GRK2**: β-subunit (GNB2) and γ-subunit (GNG2)  with β-adrenergic receptor kinase 1 (GRK2) based upon 3CIK
* **with KCTD12**:
* **with PREX1**:
* **phosphorylated**: post-translational modifications taken from Phosphosite Plus

RASD2 (Rhes) has strong homology to solved homologues, except in the GNB2 binding region (C-terminus), consequently a model could not be made.