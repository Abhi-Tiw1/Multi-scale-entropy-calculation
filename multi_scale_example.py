from multi_scale import *

"""
The main function consists of two parts:
	Scaling algorithm for time series:
		base: normal coarse grained scaling
		mov_avg: moving average scaling
		mom: momentum scaling
		mavg_mom: moving average momentum scaling
		comp: composite scaling
		c2f: coarse to fine EMD based scaling
		f2c: fine to coarse EMD based scaling
	
	function:
		get_scale_series(xt,scale=5,cg_type='base')
	Inputs:
		xt: time series
		scale: number of scales for entropy calculation
		cg_type: scaling options, mentioned above
	Output:
		x_cg: scaled series returned as list of lists
	
	Entropy type calculation:
		mspe: multiscale permutation entropy
		msse: multi-scale sample entropy
		comp_mspe: composite multiscale permutation entropy
		comp_msse: composite multiscale sample entropy
	
	function:
		calc_mse(x_cg,ent_typ='pe',typ='',mod_flag=0,deg=3,lag=1,m=2,r=0.15)
	Inputs:
		x_cg: scaled series passed as list of lists
		ent_type: entropy type as mentioned above
		Permutation entropy paarameters
		typ: '' or 'wt' for normal or weighted permutation entropy
		mod_flag: 0 or 1 for normal or modified permutation entropy
		deg: degree of motif
		lag: lag for motif
		Sample entrpoy parameters
		r: tolerance parameter
		m: embedding dimension 

 """

x_te=abs(np.random.rand(500))

#options for type of scaling:
# base, mov_avg, mom, mavg_mom, comp , f2c, c2f
x_cg=get_scale_series(x_te,scale=6,cg_type='base')

#multi-scale entropy options:
# mspe, msse, comp_mspe, comp_msse
print((calc_mse(x_cg,'msse')))

#generating names for the features:
print(get_mse_fname(cg_type='base',ent_type='msse',scale=6))