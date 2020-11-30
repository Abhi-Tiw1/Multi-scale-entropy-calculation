# Multi-scale entropy calculation

This repository implements various multi-scale entropy algorithms. The algorithm usually has two main steps.
- **Scaling**: This involves scaling the time series before calculting the entropy values for each scale. Various different methods have been implemtned for scaling.
  - base: normal coarse grained scaling
  - mov_avg: moving average scaling
  - mom: momentum scaling
  - mavg_mom: moving average momentum scaling
  - comp: composite scaling
  - c2f: coarse to fine EMD based scaling
  - f2c: fine to coarse EMD based scaling
- **Entropy calculation**: Entropy calculation for each scale has been done using two different algorithms:
  - Sample entropy algorithm: Uses embedding dimension (m) and tolerance (r) as the main parameters
  - Permutation entropy: Calculates permutation entropy with main parameters being lag and degee of motifs. Additionally two different variants of PE have also been implemented. These are the modified permutation entropy (by setting mod_flag=1) and weighted PE (by setting typ='wt')
 
## Requirements

PyEMD and nolds library: 

    $ pip install EMD-signal
    $ pip install nolds
  
## Reference
Tiwari, Abhishek, Isabela Albuquerque, Mark Parent, Jean-François Gagnon, Daniel Lafond, Sébastien Tremblay, and Tiago H Falk. "Multi-Scale Heart Beat Entropy Measures for Mental Workload Assessment of Ambulant Users." Entropy 21, no. 8 (2019): 783.

