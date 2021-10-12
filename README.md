# PAL

PAL is an earthquake detection and location architecture. <br>
<br>
The workflow can be described as: <br>
(1) phase **P**icking <br>
(2) phase **A**ssociation <br>
(3) event **L**ocation <br>
<br>
* phase **P**ickers  
*picker_pal.py* defines the default PAL phase picker. 
```python
# use picker
# 1. waveform --> picks
import picker_pal
picker = picker_pal.STA_LTA_Kurtosis()
picks = picker.pick(stream, out_pick) # input obspy.stream
```
  
* phase **A**ssociators  
*associator_pal.py* defines the default PAL phase associator.
```python
# use associator
# 2. picks --> events
import associator_pal
associator = associator_pal.TS_Assoc(sta_dict)
associator.associate(picks, out_ctlg, out_pha)
```

* event **L**ocation <br>
HypoInverse and HypoDD interfaces are provided for location purpose. <br>

***
## Installation <br>
PAL is a set of codes. All you need is to setup proper Python environment. This can be accomplished easily by installing [Anaconda](https://www.anaconda.com/products/individual#Downloads) and [Obspy](https://github.com/obspy/obspy/wiki/Installation-via-Anaconda) sequentially. Or you can use the *env/pal.yml* file with conda. 
