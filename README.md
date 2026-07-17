# Code for, "Prey age modifies risk-taking behavior in a multi-predator environment."

_Authors_:  

  - Brian J. Smith <a itemprop="sameAs" content="https://orcid.org/0000-0002-0531-0492" href="https://orcid.org/0000-0002-0531-0492" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" style="width:1em;margin-right:.5em;"/></a>
  - Tal Avgar <a itemprop="sameAs" content="https://orcid.org/0000-0002-8764-6976" href="https://orcid.org/0000-0002-8764-6976" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" style="width:1em;margin-right:.5em;"/></a>
  - Scott D. Peacor <a itemprop="sameAs" content="https://orcid.org/0000-0002-5334-7775" href="https://orcid.org/0000-0002-5334-7775" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" style="width:1em;margin-right:.5em;"/></a>
  - Daniel R. Stahler <a itemprop="sameAs" content="https://orcid.org/0000-0002-8740-6075" href="https://orcid.org/0000-0002-8740-6075" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" style="width:1em;margin-right:.5em;"/></a>
  - Matthew C. Metz <a itemprop="sameAs" content="https://orcid.org/0000-0002-7037-9891" href="https://orcid.org/0000-0002-7037-9891" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" style="width:1em;margin-right:.5em;"/></a>
  - Jack W. Rabe <a itemprop="sameAs" content="https://orcid.org/0000-0003-3227-2484" href="https://orcid.org/0000-0003-3227-2484" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" style="width:1em;margin-right:.5em;"/></a>
  - Wesley Binder <a itemprop="sameAs" content="https://orcid.org/0009-0002-4393-3647" href="https://orcid.org/0009-0002-4393-3647" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" style="width:1em;margin-right:.5em;"/></a>
  - Daniel R. MacNulty <a itemprop="sameAs" content="https://orcid.org/0000-0002-9173-8910" href="https://orcid.org/0000-0002-9173-8910" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" alt="ORCID iD icon" style="width:1em;margin-right:.5em;"/></a>


## Manuscript status
Published to bioRxiv on July XX, 2026.  https://doi.org/

Submitted for review to *Ecology Letters* on July XX, 2026.

## Code Repository

This GitHub repository contains the code needed to fit the models, generate figures, and create summary statistics described in the manuscript. The data and covariates are stored separately (see below).

### Versions 

The following DOI will always resolve to the latest version:

[![DOI](https://zenodo.org/badge/DOI/XX.svg)](https://doi.org/XX)

Other links are version-specific.

#### Repository version 0.1 -- Prior to peer review

This release (v0.1) was created before publishing preprint and before peer-review. The repository is archived here:

[![DOI](https://zenodo.org/badge/DOI/XX.svg)](https://doi.org/XX)


## Data Availability

The data are stored in a Zenodo archive ([10.5281/zenodo.21402472](https://doi.org/10.5281/zenodo.21402472)). The code assumes data are stored in a parent directory up one level from the scripts (`../elk_RT_data/`).

Outputs (e.g., intermediate processed data files) are also stored in the Zenodo archive, again, assumed to be located in a directory one level above (`../out/`).

Note there are also silhouettes from PhyloPic.org stored one level above (`../PhyloPic/`).

Thus, to recreate the workflow, users should have their directory structured like this:

```
<Parent Directory>
|-- elk_RT/
|-- elk_RT_data/
|-- out/
|-- PhyloPic/
```

where `elk_RT/` is this GitHub repo, and the other three files come from the other Zenodo archive.

### Scripts

- `00_import_data.R` -- imports data from particular directories on BJS's computer; should not be run by other users
- `01_clean_data.R` -- uses `amt` cleaning workflow to clean raw GPS data
- `02_subset_data.R` -- subsets cleaned GPS data to the seasons of interest and ensures constant (1-h) fix rates
- `03_steps.R` -- creates steps from GPS points and samples available steps
- `04_covariates.R` -- attaches all covariates to used and available steps
- `05_final_prep.R` -- calculates mean and SD for all covariates and scales/centers data (z-transform)
- `06_iSSA.R` -- first stage analysis; runs integrated step selection analysis (iSSA)
- `07_UHC.R` -- evaluates fitted models using graphical used-habitat calibration plots (UHC plots)
- `08_visualize_iSSA.R` -- creates figures from fitted iSSA
- `09_risk_model.R` -- second stage analysis; estimates risk taking as a function of covariates
- `10_visualize_risk_model.R` -- creates figures from risk-taking analysis
- `11_age_structure.R` -- estimates population-level risk taking from age structure reconstruction
- `12_graphical_model.R` -- creates graphical model shown in Box 1
- `99_fun.R` -- contains various helper functions used throughout other scripts
