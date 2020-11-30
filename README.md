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
  
  
## Citation

@article{tiwari2019multi,
  title={Multi-Scale Heart Beat Entropy Measures for Mental Workload Assessment of Ambulant Users},
  author={Tiwari, Abhishek and Albuquerque, Isabela and Parent, Mark and Gagnon, Jean-Fran{\c{c}}ois and Lafond, Daniel and Tremblay, S{\'e}bastien and H Falk, Tiago},
  journal={Entropy},
  volume={21},
  number={8},
  pages={783},
  year={2019},
  publisher={Multidisciplinary Digital Publishing Institute}
}
