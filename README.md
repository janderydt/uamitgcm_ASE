# Repository in support of De Rydt and Naughten (2024)

The repository contains scripts to

* generate model diagnostics from the raw Úa-MITgcm output:
  + mean melt rates and cavity transfer coefficients (Eq. (1) and Figs. 5, S3 and S4 in De Rydt and Naugthen (2024))
  + streamfunction amplitudes (Figs. 7 and S5 in De Rydt and Naughten (2024))
  + ice-shelf mass balance components (Figs. 9 and S6 in De Rydt and Naughten (2024))
    
  This requires the Úa-MITgcm toolbox, which is linked as a submodule. To clone the toolbox, use `git clone --recursive` instead of just `git clone`.
  
* reproduce all figures in De Rydt and Naughten (2024)

The results in De Rydt and Naughten (2024) are based on a regional Amundsen Sea Embayement configuration of the coupled ice-ocean model Úa-MITgcm. The Úa and MITgcm input files and the Úa-MITgcm config files can be found [here](https://github.com/knaughten/UaMITgcm/tree/archer2/example/ASE_999).

Unprocessed Úa-MITgcm output at 10-year intervals for the _ref\_melt_, _hi\_melt_, _av\_melt_ and _var\_melt_ experiments in De Rydt and Naughten (2024) is available [here](https://zenodo.org/records/10778033).

## References
De Rydt, J. and Naughten, K.: Geometric amplification and suppression of ice-shelf basal melt in West Antarctica, The Cryosphere, 18, 1863–1888, https://doi.org/10.5194/tc-18-1863-2024, 2024.


