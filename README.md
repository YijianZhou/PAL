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

## Installation

PAL is a set of codes. All you need is to setup proper Python environment. This can be accomplished easily by installing [Anaconda](https://www.anaconda.com/products/individual#Downloads) and [Obspy](https://github.com/obspy/obspy/wiki/Installation-via-Anaconda) sequentially. Or you can use the *env/pal.yml* file with conda. <br>

## References

- **Zhou, Y.**, H. Yue, S. Zhou, L. Fang, Y. Zhou, L. Xu, Z. Liu, T. Wang, L. Zhao, & A. Ghosh (2022). Microseismicity along Xiaojiang Fault Zone (Southeastern Tibetan Plateau) and the Characterization of Interseismic Fault Behavior. *Tectonophysics*; 833: 229364. doi: [10.1016/j.tecto.2022.229364](https://doi.org/10.1016/j.tecto.2022.229364)  

- **Zhou, Y.**, H. Yue, L. Fang, S. Zhou, L. Zhao, & A. Ghosh (2021). An Earthquake Detection and Location Architecture for Continuous Seismograms: Phase Picking, Association, Location, and Matched Filter (PALM). *Seismological Research Letters*; 93(1): 413–425. doi: [10.1785/0220210111](https://doi.org/10.1785/0220210111)  

- **Zhou, Y.**, A. Ghosh, L. Fang, H. Yue, S. Zhou, & Y. Su (2021). A High-Resolution Seismic Catalog for the 2021 MS6.4/Mw6.1 YangBi Earthquake Sequence, Yunnan, China: Application of AI picker and Matched Filter. *Earthquake Science*; 34(5): 390-398.doi: [10.29382/eqs-2021-0031](https://doi.org/10.29382/eqs-2021-0031)  

- Lu, W., **Y. Zhou**, Z. Zhao, H. Yue, & S. Zhou (2021). Aftershock sequence of the 2017 Mw 6.5 Jiuzhaigou, China earthquake monitored by an AsA network and its implication to fault structures and strength. *Geophysical Journal International*; 228(3): 1763-1779. doi: [10.1093/gji/ggab443](https://doi.org/10.1093/gji/ggab443)  
