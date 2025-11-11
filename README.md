# Advanced Bacterial Dynamics Model

## 🔬 Model Overview

This advanced computational model simulates bacterial population dynamics, integrating demographic, birth, and contact data to provide a comprehensive understanding of bacterial transmission and resistance.

## 🌟 Key Features

### Comprehensive Data Integration
- Population demographics
- Birth rate data
- Contact matrix integration
- Multi-compartment bacterial population tracking

### Detailed Modeling Aspects
- Age-stratified population dynamics
- Sensitive and resistant bacterial strain interactions
- Co-colonization mechanisms
- Transmission cost of resistance

## 📊 Model Compartments

1. **Uninfected Population** (`X`)
2. **Sensitive Bacterial Strain** (`S`)
3. **Resistant Bacterial Strain** (`R`)
4. **Treated Sensitive Strain** (`Sr`)
5. **Treated Resistant Strain** (`Rs`)

## 🛠 Technical Specifications

### Transmission Dynamics
- Strain-specific transmission rates
- Age-dependent contact patterns
- Co-colonization efficiency modeling

### Key Parameters
- Transmission rates
- Clearance rates
- Resistance transmission cost
- Co-colonization efficiency

## 📦 Dependencies

- R (version 4.x recommended)
- Packages:
  - deSolve
  - viridis
  - ggplot2
  - tidyverse
  - data.table
  - patchwork
  - here

## 🚀 Installation

### R Package Installation
\`\`\`r
# Install pacman if not already installed
install.packages("pacman")

# Install required packages
pacman::p_load(
  deSolve, 
  viridis, 
  ggplot2, 
  tidyverse, 
  data.table,
  patchwork,
  here
)
\`\`\`

## 🔧 Usage

1. Ensure data files are in the same directory:
   - \`Population_emro_2023_1yearage.csv\`
   - \`3_U_1_Birth_1year_emro.csv\`
   - \`3_U_1_contact_Pakistan_1y.csv\`

2. Run the script:
\`\`\`bash
Rscript davies_model_advanced.R
\`\`\`

## 📈 Outputs

- Bacterial population dynamics visualization
- CSV file with model results
- PNG plot of population dynamics

## 🧪 Customization

Modify the following in the script:
- \`config$country\`: Change target country
- Model parameters in \`prepare_model_parameters()\`
- Initial state conditions

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Commit your changes
4. Push to the branch
5. Create a Pull Request


[Your contact information]
