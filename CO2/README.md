# CO2L
Component library in EcosimPro for dynamic simulations of CO2-based cooling systems used at CERN for Silicon detector cooling.

Auto-generated documentation may be found [here](https://co2l.web.cern.ch)

![Example model made with CO2L](example.svg)

## Dependencies
Requires the following additional EcosimPro libraries to run:
* CRYOLIB
* MATH
* PORTS_LIB
* PLC
* THERMAL
* THERMO_TABLE_INTERP

### Additional dependencies
Requires the following additional C++ libraries:
* moistairlib: External C++ library for moist air calculations based on ASHRAE equations. This static library needs to be compiled for the C++ compiler you have on your machine. For this reason, the static library files are not included in this repo. Please get in touch with me at viren.bhanot@cern.ch to discuss compiling these files. You will need these to run the plate fin heat exchanger models.
* refpropecosimpro: call the refprop C++ wrapper from ecosimpro. The full refprop suite has NOT been implemented. Instead, only the functions necessary have been implemented, and only for CO2 as a refrigerant.

## Files Tracked:
  
| Element          | Associated File                                     |
|------------------|-----------------------------------------------------|
| Library          | LIB_PATH\.lsp.xml                                   |
|                  | LIB_PATH\etc\.ler.xml                               |
| EL source code   | LIB_PATH\sources.el                                 |
| Schematic        | LIB_PATH\schematics\.eds                            |
|                  | LIB_PATH\schematics\.edx                            |
|                  | LIB_PATH\schematics\.prop.xml                       |
|                  | LIB_PATH\schematics\component_images\.png           |
|                  | LIB_PATH\schematics\component_csv\.csv              |
|                  | LIB_PATH\schematics\component_xml\.xml              |
| Filter file list | LIB_PATH\filters\.flt                               |
| Symbol           | LIB_PATH\palette\.sym                               |
|                  | LIB_PATH\etc\.lsc.xml                               |
|                  | LIB_PATH\palette\symbol_images\.png                 |
| Partition        | LIB_PATH\experiments\PAR_FOLDER\.par.xml            |
|                  | LIB_PATH\experiments\PAR_FOLDER\prop.xml            |
| Experiment       | LIB_PATH\experiments\PAR_FOLDER\EXP_FOLDER\.exp     |
|                  | LIB_PATH\experiments\PAR_FOLDER\EXP_FOLDER\.cfg.xml |
