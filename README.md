# PAD

Package for detecting earthquakes from raw continuous waveform. <br>
<br>
procudures include: <br>
(1) phase picking <br>
(2) phase association <br>
<br>
Each of the three procedures are implemented in seperate scripts, i.e. the 'pickers.py' and 'detectors.py'. An example for combining these processes are shown in 'mk_ctlg.py', which aims to get an earthquake catalog directly from raw waveforms. 'parallel.py' are also provided for parallel computing.
<br>
  
* phase pickers  
*pickers.py* defines various picking algorithms as picker classes. 
```python
# use picker
# 1. waveform --> picks
import pickers
picker = pickers.Trad_PS()
picks = picker.pick(stream) # input obspy.stream
```
  
* earthquake detectors  
*detectors.py* defines various earthquake detection methods.
```python
# use detector
# 2. associate by original time (ot) cluster: picks --> events
event_picks = detector.pick2event(picks)
```

```python
# 3. associate by earthquake location (P travel time cluster)
event_loc, event_pick = detector.locate(event_pick)
# 4. estimate magnitude
event_loc_mag = detector.calc_mag(event_pick, event_loc)
```
