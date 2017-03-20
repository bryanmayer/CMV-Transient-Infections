# Description of scripts

1. setup_data.R - process the data to identify infection episodes (transient infections) with duration.

2. CMV_transient_descriptives.Rmd -> Basic descriptives for CMV, compares transient infections and primary infection

3. CMV_blip_comparison.Rmd  -> Basic descriptives all viruses, compares percent positive, frequencies. Also compares blips to demographic variables.

4. analyze_blip_prob.Rmd-> MLE model of probability of blips

#Other Documents

*questionable_blips.Rmd - old file to compare questionable blips between me and Elizabeth, currently obsolete with updated defintions.

# Scripts in old_scripts/

*all_virus_blip_analysis.R - old version of CMV_blip_comparison.Rmd

*blip_duration_analysis.R - something old for analyzing duration_data, not sure what is does

*cluster_clip_analysis.R - Analyzed how positive swabs clustered: ++, -+, --, +- probabilities. Old and probably would need an update to be useful.

*CMV_blip_descriptives.R - old messy version of CMV_blip_descriptives.Rmd.

*transmission_analysis.R - looked at correlation of demofactors with blip risk, including consecutive blip probabilities. Old and better done in CMV_blip_comparison.Rmd
