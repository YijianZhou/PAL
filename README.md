# PpkAssocLoc

Package for processing raw continuous waveform. <br>
<br>
procudures include: <br>
(1) phase picking <br>
(2) associate picks to events <br>
(3) locate events and estimate magnitude. <br>
<br>
Each of the three procedures are implemented in seperate scripts, i.e. the 'pickers.py', 'associators.py' and 'locators.py'. An example for combining these processes are shown in 'mkctlg.py', which aims to get an earthquake catalog directly from raw waveforms. 'parallel.py' are also provided for parallel computing.
<br>
* phase pickers
    pickers.py defines various picking algorithms as picker classes. 
```python
# use picker
# 1. waveform --> picks
import pickers
picker = pivkers.Trad_PS()
picks = picker.pick(stream) # input obspy.stream
```
<br>
    associators.py defines various methods to associate phase picks to picks of different events.
```python
# use associator
# 2. associate: picks --> events
event_picks = associator.pick2event(picks)
# write pahse file
associator.write(event_picks, out_pha)
```
<br>
    locators.py defines various earthquake locate methods.
```python
# use locator
# 3. locate evnets
event_loc = locator.locate(event_pick)
# 4. estimate magnitude
event_loc_mag = locator.calc_mag(event_pick, event_loc)
# write catalog
locator.write(event_loc_mag, out_file)
```
