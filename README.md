# PpkAssocLoc

Package for processing raw continuous waveform. procudures include: (1) phase picking, (2) associate picks to events, and (3) locate events and estimate magnitude. Each of the three procedures are implemented in seperate classes, i.e. the 'pickers.py', 'associators.py' and 'locators.py'. An example for combining these processes are shown in 'mkctlg.py', which aims to get an earthquake catalog directly from raw waveforms. 'parallel.py' are also provided for parallel computing.
