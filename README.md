# :sunflower: Strand-seq Channel plots
R functions for combinatorial data channel plotting of scTRIP single cell Strand-seq data.

# Usage
1. Clone the repository

```bash
git clone https://github.com/pweidner/tri-channel-plots
```
2. Open the run_script.R and run the plot_channels function

Input files:
- d <- scTRIP counts output file (e.g. AGLCD.txt.gz)
- d.hap <- run haplotagger on bam files in the haplotag folder or import such table

Function Parameters:
- cell_ID (as in your table, e.g. "i504" or "P1530_i504")
- chromosome (e.g. "chr2")
- channels (1=State, 2=Depth, 3=Haps (merged), 4 = Haps (split))
- plot_range (e.g. c(40000000,120000000))
- count = FALSE (e.g. show count plot instead of W:C ratio)
- roi (e.g. c(7500000,85000000)

# Roadmap

# 📕 Technical features
- [X] Make yaml file for conda env to launch into RStudio server on the cluster
- [X] Option to split H1 and H2 in seperate plots
- [ ] Option to plot sv call track along the plotting region
- [X] General option to chose combination of channel plots through numeric ID (e.g. 3=classical trichannel, 4= same but split haps etc.)
- [X] Update ggarrange composition based on channels ID after including merged haps
- [X] Region of interest highlighting
- [X] Add Counts channel (hisogram plot)
- [ ] Unified title including set params
- [ ] coord_flip = T option to plot horizontal version
- [ ] gene = "CD3" param to zoom in on genomic ranges of interest
- [ ] margin = 2 param exclusive to plot_range that set a margin to plot around roi/gene
- [ ] loop through all cells of samples for given ranges for manual curation of hotspots (sample param)
- [ ] parallelization for haplotagger

# 🛑 Small issues
- [X] filtering for whole chromosomes when plot_range=NULL corrently not working
- [X] Scaling
- [X] Input parameter checks especially for channels beeing either 1,2,3 or 4 for now!
- [X] Formatting of breaks in smaller regions (scale_mb function)
- [X] subsetting is slow


Example output plot:

![trichannelplot](example_cell.png)

# Authors
- [Patrick Weidner](https://github.com/pweidner)
- [Suharto Bannerjee](https://github.com/suhartobanerjee)
- David Porubsky wrote the bamregion2ranges.R function! :)

:mailbox_with_mail:  Please contact us with any problems or submit them as an issue in this Github repository.
