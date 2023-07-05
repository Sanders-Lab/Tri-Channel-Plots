# Strand-seq Channel plots
R functions for combinatorial data channel plotting of scTRIP single cell Strand-seq data.

# Roadmap

# 📕 Technical features
- [X] Make yaml file for conda env to launch into RStudio server on the cluster
- [X] Option to split H1 and H2 in seperate plots
- [ ] Option to plot sv call track along the plotting region
- [X] General option to chose combination of channel plots through numeric ID (e.g. 3=classical trichannel, 4= same but split haps etc.)
- [ ] Update ggarrange composition based on channels ID after including merged haps
- [X] Region of interest highlighting
- [ ] Add Counts channel (hisogram plot)
- [ ] loop through all cells of samples for given ranges for manual curation of hotspots

# 🛑 Small issues

- [X] Scaling
- [ ] Input parameter checks especially for channels beeing either 1,2,3 or 4 for now!
- [ ] Formatting of breaks in smaller regions (scale_mb function)
- [ ] subsetting is slow
- [ ] parallelization for haplotagger

Example output plot:

![trichannelplot](example_cell.png)

# Authors
- [Patrick Weidner](https://github.com/pweidner)
- [Suharto Bannerjee](https://github.com/suhartobanerjee)
- David Porubsky wrote some of the functions were sourcing!

Please contact us with any problems or submit them as an issue in this Github repository.
