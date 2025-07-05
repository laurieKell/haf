# haf: FLR-based Implementation of ICES Harvest Control Rules

[![R-CMD-check](https://github.com/yourusername/haf/workflows/R-CMD-check/badge.svg)](https://github.com/yourusername/haf/actions)
[![License: GPL-3](https://img.shields.io/badge/License-GPL%203-green.svg)](https://opensource.org/licenses/GPL-3.0)

A comprehensive R package implementing ICES (International Council for the Exploration of the Sea) harvest control rules within the FLR (Fisheries Library for R) framework. The package provides standardized implementations of various ICES advice rules including hockey-stick HCRs, RFB rules for data-limited stocks, and Nephrops-specific rules.

## Installation

```r
# Install from GitHub (requires devtools)
devtools::install_github("yourusername/haf")

# Or install from local source
install.packages("path/to/haf", repos = NULL, type = "source")
```

## Dependencies

The package requires the following R packages:
- **FLCore** (>= 2.6.0) - Core FLR functionality
- **FLBRP** (>= 2.5.0) - FLR reference points
- **ggplot2** (>= 3.4.0) - Plotting functionality

## Quick Start

```r
library(haf)

# Load example data
data(ple4)

# Set up HCR parameters for hockey-stick rule
params <- FLPar(
  blim = 100000,   # Biomass limit reference point
  btrig = 150000,  # Biomass trigger reference point  
  bmin = 100000,   # Minimum biomass
  ftar = 0.3,      # Target fishing mortality
  fmin = 0.05      # Minimum fishing mortality
)

# Run HCR simulation
result <- hcrICESAR(ple4, eql, sr_deviances, params, 
                  start = 2010, end = 2030)

# Create HCR visualization
hcrPlot2(btarget = 1.0, blim = 0.1, ftarget = 0.8)
```

## Available Harvest Control Rules

### 1. Hockey-Stick HCR (Category 1 & 2 stocks)

The standard ICES approach for well-assessed stocks:

```r
# Standard hockey-stick HCR
result <- hcrICESAR(stock, eql, sr_deviances, params)
```

**Features:**
- Biomass-based trigger system
- Linear reduction between Blim and Btrigger
- Optional TAC bounds and implementation error
- Support for multiple iterations

### 2. RFB Rule (Category 3 stocks)

For data-limited stocks using the r × f × b × m approach:

```r
# Set up RFB parameters
cntrl <- rfbParams(linf = 35, lc = 25, k = 0.25, index = historical_index)

# Apply RFB rule
advice <- rfbRule(iYr = 2024, indx = biomass_index, 
                  FIndx = length_index, tac = previous_tac, cntrl = cntrl)
```

**Components:**
- **r**: Biomass trend ratio (recent vs historical)
- **f**: Length-based fishing pressure proxy
- **b**: Biomass safeguard
- **m**: Precautionary multiplier

### 3. Survey-Based HCR

Generic survey-based management for stocks managed through survey indices:

```r
# Generic survey-based HCR
cntrl <- FLPar(fmsy = 0.16, btrig = 5000, bbuf = 2500)
advice <- hcrSurvey(iYr = 2024, index = survey_index, 
                    catch = historical_catch, cntrl = cntrl)

# Nephrops-specific HCR (wrapper function)
advice <- hcrNephrops(iYr = 2024, index = survey_index, 
                      catch = historical_catch, cntrl = cntrl)
```

**Zones:**
- **Green**: B ≥ Btrigger → F = FMSY
- **Yellow**: Bbuffer ≤ B < Btrigger → F = FMSY × reduction factor
- **Red**: B < Bbuffer → F = 0.8 × previous catch

**Suitable for:**
- Nephrops and other crustacean stocks
- Stocks with limited analytical assessments
- Species managed primarily through survey indices
- Data-limited stocks requiring precautionary management

## Visualization

The package includes comprehensive plotting functions:

```r
# Standard HCR plot
hcrPlot(btarget = 1.0, blim = 0.1, ftarget = 0.8)

# Enhanced HCR plot with management zones
hcrPlot2(btarget = 1.0, blim = 0.1, ftarget = 0.8, fmsy = 1.0)

# Survey-based HCR visualization
plotSurveyHCR(cntrl, index, catch)

# Nephrops HCR visualization
plotNephropsHCR(cntrl, index, catch)
```

## Documentation

- **Package vignette**: `vignette("haf")`
- **Function documentation**: `?hcrICESAR`, `?rfbRule`, `?hcrNephrops`
- **Examples**: See individual function help pages

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## Citation

If you use haf in your research, please cite:

```r
citation("haf")
```

## License

This package is licensed under the GNU General Public License v3.0. See the [LICENSE](LICENSE) file for details.

## References

- ICES (2025). Technical Guidelines for Data-limited Stock Assessments. ICES CM 2025/ACOM:68
- ICES (2023). ICES Technical Guidelines for Category 1 and 2 stocks. ICES CM 2023/ACOM:XX
- Kell, L.T., et al. (2023). FLR: Fisheries Library for R. Journal of Statistical Software, 58(12), 1-34.

## Support

For questions, bug reports, or feature requests, please:
- Open an issue on GitHub
- Contact the maintainer at your.email@example.com 