# PAL

PAL is an earthquake detection architecture. <br>
<br>
The workflow can be described as: <br>
(1) phase **P**icking <br>
(2) phase **A**ssociation <br>
(3) earthquake **L**ocation <br>
<br>
* phase **P**ickers  
*pickers.py* defines phase picking algorithms. 
```python
# use picker
# 1. waveform --> picks
import pickers
picker = pickers.STA_LTA_PCA()
picks = picker.pick(stream, out_pick) # input obspy.stream
```
  
* phase **A**ssociators  
*associators.py* defines phase associate methods.
```python
# use associator
# 2. picks --> events
import associators
associator = associators.TS_Assoc(sta_dict)
associator.associate(picks, out_ctlg, out_pha)
```

* earthquake **L**ocation <br>
HypoInverse and HypoDD interfaces are provided for location purpose
