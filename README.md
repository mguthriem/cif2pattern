This is a script intended to be run in [mantidworkbench](https://www.mantidproject.org/) in order to quantitaively simulate measured diffraction data from samples measured on the SNAP TOF neutron diffractometer. The simulation attempts to accurately replicate signal to background and experimental noise using a reference data set. At present, only a single experimental configuration is supported: 

* A powder sample measured a diamond anvil cell with 90 degree detector setting.

In this setting, background levels, noise levels and peak widths were empirically fitted to a reference dataset of nickel powder measured in a DAC.

To run, open `cif2pattern.py` in an instance of mantidworkbench, configure as needed (following instructions in script), then execute.
