# PpkAssocLoc

Package for processing raw continuous waveform. <br>
<br>
procudures include: <br>
(1) phase picking <br>
(2) associate picks to events <br>
(3) locate events and estimate magnitude. <br>
<br>
Each of the three procedures are implemented in seperate scripts, i.e. the 'pickers.py', 'associators.py' and 'locators.py'. An example for combining these processes are shown in 'mkctlg.py', which aims to get an earthquake catalog directly from raw waveforms. 'parallel.py' are also provided for parallel computing.
