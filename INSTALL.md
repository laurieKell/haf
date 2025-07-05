# haf Installation Guide

## Prerequisites

Before installing haf, you need to have the following R packages installed:

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

- **devtools** - For installing from GitHub
- **testthat** - For running tests
- **knitr** - For building vignettes
- **rmarkdown** - For documentation

## Installation Methods

### Method 1: Install from GitHub (Recommended)

```r
# Install devtools if not already installed
if (!require(devtools)) install.packages("devtools")

# Install haf from GitHub
devtools::install_github("yourusername/haf")
```

### Method 2: Install from Local Source

If you have downloaded the source code:

```r
# Navigate to the package directory
setwd("path/to/haf")

# Install the package
install.packages(".", repos = NULL, type = "source")
```

### Method 3: Install Dependencies First

If you encounter dependency issues:

```r
# Install FLR packages first
install.packages(c("FLCore", "FLBRP"), repos = "https://flr-project.org/R")

# Then install haf
devtools::install_github("yourusername/haf")
```

## Verification

After installation, verify that the package works correctly:

```r
# Load the package
library(haf)

# Check version
packageVersion("haf")

# Run tests
testthat::test_package("haf")
```

## Troubleshooting

### Common Issues

1. **FLR packages not found**: Make sure you have the FLR repository added:
   ```r
   install.packages("FLCore", repos = "https://flr-project.org/R")
   ```

2. **Build tools required**: On Windows, you may need Rtools:
   ```r
   install.packages("devtools")
   devtools::find_rtools()
   ```

3. **Permission errors**: Try running R as administrator or use:
   ```r
   devtools::install_github("yourusername/haf", force = TRUE)
   ```

## Development Installation

For development work:

```r
# Clone the repository
git clone https://github.com/yourusername/haf.git

# Set working directory
setwd("haf")

# Install in development mode
devtools::install()
```

## Updating

To update to the latest version:

```r
# Remove old version
remove.packages("haf")

# Install new version
devtools::install_github("yourusername/haf")
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