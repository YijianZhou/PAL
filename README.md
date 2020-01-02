# PAD

Earthquake detection from raw continuous waveform. <br>
<br>
procudures include: <br>
(1) phase picking <br>
(2) phase association <br>
<br>
Both procedures are implemented in seperate scripts, i.e. the 'pickers.py' and 'associators.py'. An example for combining these two processes for earthquake detection are shown in 'run_ppk_assoc.py'. 'parallel_ppk_assoc.py' are also provided for parallel computing.
<br>
  
* phase pickers  
*pickers.py* defines various phase picking algorithms. 
```python
# use picker
# 1. waveform --> picks
import pickers
picker = pickers.Trad_PS()
picks = picker.pick(stream, out_ppk) # input obspy.stream
```
  
* phase associators  
*associators.py* defines various phase associate methods.
```python
# use associator
# 2. picks --> events
import associators
associator = associators.TS_Assoc(sta_dict, resp_dict)
associator.associate(picks, out_ctlg, out_pha)
```
