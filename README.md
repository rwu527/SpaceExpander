# SpaceExpander

<p align="center">
  <a href="https://www.lilab-ecust.cn/markushclaim" target="_blank">
    <img src="https://github.com/rwu527/SpaceExpander/raw/main/src/logo.png" alt="SpaceExpander Logo" width="300">
  </a>
</p>

<p align="center">
  <b>https://www.lilab-ecust.cn/markushclaim</b>
</p>
<p align="center">
  <b>An innovative solution for automated drafting and expansion of Markush claims.</b>
</p>

---

SpaceExpander is a web-based platform that allows users to input molecular structures and automatically generate Markush claims.  
It supports intelligent substructure analysis and systematic enumeration of chemical space for patent drafting and exploration.

**Offline Usage**: Standalone executables are available for Windows, macOS, and Linux platforms. These provide full functionality without dependencies or installation.  
Download here:  
https://drive.google.com/drive/folders/1DNeaFxnVlkqT0lW7YxIK2sOQMCxz6SkM?usp=drive_link

---

## Features

- Input molecular structures via SMILES or supported file formats
- Automatically generate Markush claims based on molecular analysis
- Intelligent substructure recognition and classification
- Offline usage supported with standalone executables for major platforms

## Installation

```bash
git clone git@github.com:rwu527/SpaceExpander.git
cd SpaceExpander
conda env create -f environment.yml
conda activate SpaceExpander
```

## Usage

Run the application with:

```javascript
python main.py
```

You can use any file from the test_dataset folder as an example.
The generated results are stored in the directory under the input path.
