# PAL

PAL is an earthquake detection architecture. <br>
<br>
The workflow can be described as: <br>
(1) phase **P**icking <br>
(2) phase **A**ssociation <br>
(3) earthquake **L**ocation <br>
<br>
Both procedures are implemented in seperate scripts: the *pickers.py* and *associators.py*. Main function for processing raw waveform is *run_pick_assoc.py*; for parallel computing, use *parallel_pick_assoc.py*. For tunning association parameters, use *run_assoc.py* or *parallel_assoc.py*
<br>
  
* phase pickers  
*pickers.py* defines various phase picking algorithms. 
```python
# use picker
# 1. waveform --> picks
import pickers
picker = pickers.STA_LTA_PCA()
picks = picker.pick(stream, out_pick) # input obspy.stream
```
  
* phase associators  
*associators.py* defines various phase associate methods.
```python
# use associator
# 2. picks --> events
import associators
associator = associators.TS_Assoc(sta_dict)
associator.associate(picks, out_ctlg, out_pha)
```
