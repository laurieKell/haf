# FLRICES Installation Guide

## Prerequisites

Before installing FLRICES, you need to have the following R packages installed:

### Required Packages

1. **FLCore** (>= 2.6.0) - Core FLR functionality
   ```r
   install.packages("FLCore", repos = "https://flr-project.org/R")
   ```

2. **FLBRP** (>= 2.5.0) - FLR reference points
   ```r
   install.packages("FLBRP", repos = "https://flr-project.org/R")
   ```

3. **ggplot2** (>= 3.4.0) - Plotting functionality
   ```r
   install.packages("ggplot2")
   ```

### Optional Packages (for development)

- **testthat** (>= 3.0.0) - For running tests
- **knitr** - For building vignettes
- **rmarkdown** - For building documentation

## Installation Methods

### Method 1: Install from GitHub (Recommended)

```r
# Install devtools if not already installed
if (!require(devtools)) install.packages("devtools")

# Install FLRICES from GitHub
devtools::install_github("yourusername/FLRICES")
```

### Method 2: Install from Local Source

1. Download or clone the repository
2. Open R in the package directory
3. Run:
   ```r
   install.packages(".", repos = NULL, type = "source")
   ```

### Method 3: Install using R CMD INSTALL

From the command line in the package directory:
```bash
R CMD INSTALL .
```

## Verification

After installation, verify that the package works correctly:

```r
# Load the package
library(FLRICES)

# Check package version
packageVersion("FLRICES")

# Run basic tests
testthat::test_package("FLRICES")
```

## Troubleshooting

### Common Issues

1. **FLCore not found**: Make sure you have the FLR repository added:
   ```r
   options(repos = c(CRAN = "https://cran.r-project.org/",
                    FLR = "https://flr-project.org/R"))
   ```

2. **Dependencies missing**: Install missing dependencies:
   ```r
   install.packages(c("FLCore", "FLBRP", "ggplot2"))
   ```

3. **Build tools missing**: On Windows, you may need Rtools:
   - Download from: https://cran.r-project.org/bin/windows/Rtools/
   - Add to PATH environment variable

### Platform-Specific Notes

#### Windows
- Install Rtools for building packages from source
- Ensure R is in your PATH

#### macOS
- Install Xcode command line tools: `xcode-select --install`
- May need to install additional libraries via Homebrew

#### Linux
- Install build essentials: `sudo apt-get install build-essential` (Ubuntu/Debian)
- Install R development headers: `sudo apt-get install r-base-dev`

## Development Installation

For development work, install in development mode:

```r
# Install in development mode
devtools::install(".", reload = TRUE)

# Or use load_all for faster development cycle
devtools::load_all()
```

## Building Documentation

To build the package documentation:

```r
# Build documentation
devtools::document()

# Build vignettes
devtools::build_vignettes()

# Check package
devtools::check()
```

## Support

If you encounter installation issues:

1. Check the troubleshooting section above
2. Ensure all dependencies are properly installed
3. Try installing from a clean R session
4. Report issues on the GitHub repository

For additional help, contact the maintainer at your.email@example.com 