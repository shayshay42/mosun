This folder is a complete SimBiology export bundle for Julia reconstruction.

Key files:
- model_summary.txt: counts and active configset summary
- compartments/species/parameters/reactions/rules/events/observables TSVs
- variants.tsv and doses.tsv: flattened variant content and dose metadata
- stoich_sparse.tsv: sparse stoichiometry entries with indices and names
- configset_active_dump.txt: active SimBiology configset dump
- TDBr26_6_paper.sbml: SBML export (if sbmlexport succeeds)
- equations/*.txt + equation_manifest.tsv: getequations text per scenario

All dose amounts in phase1 representative scenarios are converted mg->ug/kg using BW=70 kg.
