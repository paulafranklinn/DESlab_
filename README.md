# DESlab_

Welcome to the **DESlab_** project! This repository contains all the necessary scripts and instructions to run the DESlab_ library. Please follow the instructions below to set up your environment and run the library.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
  - [Installing Miniconda](#installing-miniconda)
  - [Setting up Environment](#setting-up-environment)
  - [Platform-specific Configuration](#platform-specific-configuration)
  - [Installing DESlab](#installing-deslab)
- [Testing](#testing)
- [Troubleshooting](#troubleshooting)

## Prerequisites

- Git (for cloning the repository)
- Terminal/Command Prompt access
- Internet connection for downloading dependencies

## Installation

### Installing Miniconda

Choose the appropriate installation method for your operating system:

#### Windows

Open the terminal and run:
```sh
curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe --output .\Downloads\Miniconda3-latest-Windows-x86_64.exe
```

#### macOS

##### Apple Silicon (M1/M2/M3 Macs)
```sh
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh    
bash ~/Miniconda3-latest-MacOSX-arm64.sh
```

##### Intel Macs
```sh
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh    
bash ~/Miniconda3-latest-MacOSX-x86_64.sh
```

#### Linux
```sh
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash ~/Miniconda3-latest-Linux-x86_64.sh
```

> 💡 **Note**: For detailed installation instructions, visit the [official Miniconda documentation](https://www.anaconda.com/docs/getting-started/miniconda/install#linux-terminal-installer).

### Setting up Environment

1. **Clone the repository:**
   ```sh
   git clone https://github.com/paulafranklinn/DESlab_
   cd DESlab_
   ```

2. **Create the conda environment:**
   ```sh
   conda env create -f environment.yml
   ```

3. **Activate the environment:**
   ```sh
   conda activate DESlab
   ```

### Platform-specific Configuration

#### Windows Configuration

**Setting up Environment Variables:**

1. Open System Properties:
   - Go to Control Panel → System → Advanced System Settings
   - In the Advanced tab, click "Environment Variables"

2. Edit the PATH variable in System Variables:
   - Click "Edit" and look for paths separated by semicolons (`;`)
   - Remove any paths for older Python versions
   - Add the path for Python 3.12 (e.g., `C:\...\Python\Python312\`)

   **Examples:**
   - ✅ Correct: `C:\...\Python\Python312\`
   - ❌ Incorrect: `C:\...\Python\Python312\python.exe`

3. **Install DESlab dependencies:**
   - Run the `Install.bat` script located in the DESlab folder
   - This installs:
     - Graphviz (version 2.28)
     - FaDo (version 2.2.0)
     - Future (version 1.0.0)
     - Networkx (version 2.8.8)
     - Pyparsing (version 3.1.4)
     - Pydot (version 3.0.1)
     - DESlab

4. **Verify Graphviz Path:**
   - Ensure the Graphviz path (e.g., `...\Graphviz\bin`) is added to Windows Environment Variables

#### macOS/Linux Configuration

1. **Install Graphviz:**
   ```sh
   pip install graphviz
   ```

2. **Install LaTeX (required for PDF generation):**

   **For Linux:**
   ```sh
   sudo apt install texlive-full
   ```

   **For macOS:**
   ```sh
   brew install --cask mactex
   ```

3. **Installing DESlab**

  Navigate to the DESlab directory and install:
  ```sh
  cd DESlab_
  pip install setup.py
  ```

## Testing

To verify that everything is working correctly, run the test file:
```sh
python testando123.py
```

If the test runs without errors, your installation is complete!

## Troubleshooting

### Common Issues

- **Environment activation fails**: Ensure Miniconda is properly installed and added to your PATH
- **Graphviz not found**: Verify that Graphviz is installed and its binary path is in your system's PATH
- **LaTeX errors**: Make sure you have a complete LaTeX installation (texlive-full or MacTeX)
- **Permission errors**: On Unix systems, you may need to use `sudo` for system-wide installations

### Getting Help

If you encounter issues:
1. Check that all prerequisites are installed
2. Verify your Python version is compatible
3. Consult the project's issue tracker on GitHub
