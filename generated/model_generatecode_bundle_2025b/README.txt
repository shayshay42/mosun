This folder contains MATLAB code generated from the SimBiology project
using SimBiology.Model/generateCode in MATLAB R2025b.

Files:
- models/model_<scenario>.m : full model reconstruction code (components + defaults + configset)
- scenario_manifest.tsv : active variant/dose context used for each generated model file

Notes:
- generateCode captures all model components (species/reactions/rules/variants/doses).
- Active flags in each generated model file match the scenario context at export time.
- This code is intended as a lossless intermediate for Julia translation.
