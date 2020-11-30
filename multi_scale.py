import numpy as np
from scipy.stats import rankdata
import itertools
from PyEMD import EMD
import nolds

"""Complete package for calculating any kind of multiscale entropy features """

"""Coarse graining methods
- Normal
- Moving average (X)
- Volatility series (X)
- Moving average volatility series (X) 
- EMD- Coarse to fine series (X) --> undetermind scale --> put limitation on it
- EMD - Fine to coarse series (X) --> undetermind scale --> put limitation on it
- Composite coarse graining (X) 

Entropy measurement methods
THe permutation methods mentioned are returned in the same way as mentioned above
- Permutation entropy (X)
- Modified permutation entropy (X)
- Weighted permutation entropy (X)
- Weighted modified PE (X)
- Sample entropy (X)
- Composite variation for above (X)"""


def add_perm(perm):
	"""Add extra permutations for modified PE case """
	perm.append((1,1,0))
	perm.append((0,1,1))
	perm.append((1,0,1))
	perm.append((0,1,0))
	perm.append((0,0,0))
	perm.append((0,0,1))
	perm.append((1,0,0))
	
	
	return perm
	
def get_perms(m,mod_flag=0):
	"""get all the permutation for entropy calculation """
	perm = (list(itertools.permutations(range(m))))
	#adding similar instances	
	if mod_flag==1:
		perm=add_perm(perm)
		
	perm=np.array(perm)
	return np.array(perm)

def get_1s_pe(x,m=3,lag=1,mod_flag=0,typ=''):
	"""All the combinations of permutation entropy for a single scale """
	mot_x, wt_x=make_mot_series(x,m,lag,mod_flag=0)
	n_p=len(get_perms(m,mod_flag))
	dist=get_mot_dist(mot_x,n_p,wt_x,typ=typ)
	pe=perm_ent(dist)
	
	return np.array(pe)
	
	


def make_mot_series(time_series,m=3,lag=1,mod_flag=0):
	"""Creates a motif series and returns their with the motif distribution
	Input: 
	- time_series
	- m: permutaiton degree, lag: permutation lag
	- mod_flag: flag to use modfied PE
	Output: 
	- motif time series, 
	- corrsponding weights 
	"""
	
	time_series=np.array(time_series).squeeze()
	n=len(time_series)
	mot_x, wt_x, mod_mot_x=[], [], []

	perms=get_perms(m,0)
	perms_mod=get_perms(m,1)
	
	for i in range(n - lag * (m - 1)):
		smp=time_series[i:i + lag * m:lag]
		wt=np.var(smp)
		#orginal dense ranking of data
		mot_array1 = np.array(rankdata(smp, method='dense')-1)
		val=np.where(np.all(perms==mot_array1,axis=1))[0]
		val_mod=val
		
		if val.shape[0]==0:
			mot_array = np.array(rankdata(smp, method='ordinal')-1)
			val=np.where(np.all(perms==mot_array,axis=1))[0]
			val_mod=np.where(np.all(perms_mod==mot_array1,axis=1))[0]
			
		mot_x.append(val[0])
		mod_mot_x.append(val_mod[0])
		wt_x.append(wt)
	
	if mod_flag==1:
		return np.array(mod_mot_x), np.array(wt_x)
	elif mod_flag==0:
		return np.array(mot_x), np.array(wt_x)


def get_mot_dist(mot_x,n_p,wt_x,typ=''):
	"""Create the distribution of motifs
	Input: 
	- mot_x: Motif time series, 
	- n_p: number of permutations,
	- wt_x: weight time series
	- typ: type of entropy, normal: '', or weighted: 'wt' 
	Output:  
	- motif distribution
	"""
	
	mot_dist = [0] * n_p
	for j in range(n_p):
		if typ=='wt':
			wts=wt_x[np.where(abs(mot_x-j)==0)]
			num_mots=np.ones(len(np.where(abs(mot_x-j)==0)[0]))
			mot_dist[j]=sum(np.multiply(num_mots,wts))
			
		else:
			mot_dist[j] = len(np.where(abs(mot_x-j)==0)[0])
	#removing non occring patterns as it breaks entropy
	if len(mot_x)==0:
		mot_dist=np.zeros(n_p)*np.nan
		
	return mot_dist
	
def perm_ent(mot_dist,m=3):
	"""Returns permutation entropy for the motif distribution given --> basic function for permutation entropy """ 
	c=mot_dist
	c = [element for element in c if element != 0]
	p = np.divide(np.array(c), float(sum(c)))
	pe = -sum(p * np.log(p))
	return pe#/np.log(factorial(m))


def get_mot_ent_dist(RRs,m,lag,typ='',mod_flag=0):
	
	""" 
	#RR series for all the scales (list of lists)
	Returns four kind of motif distributions
	--> normal motif distribution ('' + 0)
	--> modified motif distribution ('' + 1)
	--> weighted motif distribution ('wt' + 0)
	--> weighted modified motif distribution ('wt' + 1)
	"""
	
	dist=[]
	for rr in RRs:
		mot_x , wt_x=make_mot_series(rr,m,lag,mod_flag = mod_flag)
		n_p=len(get_perms(m,mod_flag))
		mot_dist=get_mot_dist(mot_x,n_p,wt_x,typ='')
		dist.append(mot_dist)   #Contains motif distribution for all the different scales

	d_all=[dist]
	
	return d_all
	
def ord_dist(mot_dist_x,mot_dist_y):
	"""Returns ordinal distance between two motif distributions
	Not used anywhere in the code """
	c_x=mot_dist_x
	c_y=mot_dist_y
	m=len(c_x)
	p_x=np.divide(np.array(c_x), float(sum(c_x)))
	p_y=np.divide(np.array(c_y), float(sum(c_y)))
	sq_diff=0
	for j in range(m):
		sq_diff=sq_diff+(p_x[j] -p_y[j])**2
		
	dm=np.sqrt(m/(m-1))*np.sqrt(sq_diff)
	return dm	
	

def get_com_mspe(distS,scale,mspe):
	"""Calculate center of mass entropy using ordinal distances as weights
	NOT USED ANYWHERE IN THE CODE
	"""
	
	distS=np.array(distS)
	dm_mat=np.zeros((scale,scale))
	
	for i  in range(0,scale-1):
		for j  in range(i,scale):
			dm=ord_dist(distS[i],distS[j])
			dm_mat[i,j]=dm
			dm_mat[j,i]=dm


	com_wts=np.zeros(scale)
	for i in range(0,scale):
		com_wts[i]=np.sum(dm_mat[i,:])/(scale-1)

	com_mspe=np.sum(np.multiply(com_wts,mspe))/np.sum(com_wts)
	
	return com_mspe

def calc_mspe(distS):
	"""Calculates the scaled permutation entropy and thier oridnal avg and normal average"""
	"""Takes an input which is a list of lists where distS[i] is motif dist with scale i """
	
	mspe=[]
	scale=len(distS)
	for s in range(0,scale):
		distm=distS[s]
		pe=perm_ent(distm)
		mspe.append(pe)

	mspe=np.array(mspe)
	#com_mspe=get_com_mspe(distS,scale,mspe)
	mspe_fin=np.hstack((mspe))
	
	return mspe_fin

def scale_series(x,scale,cg_typ):

	"""Get the different scales of the series based on specific scaling type 
	Except: composite and emd scaling types
	Input: 
	Time series (x), number of scales (scale), scale type (cg_type)
	 """
	x_scale=[]
	if cg_typ=='base':
		for i in range(0,len(x),scale):
			#not divided by scale
			if i+scale<=len(x):
				val=np.sum(x[i:i+scale])/len(x[i:i+scale])
				x_scale.append(val)
				
	elif cg_typ=='mov_avg':
		wts = np.ones(scale) / scale
		val=np.convolve(x, wts, mode='valid')
		x_scale.append(val)
			
	elif cg_typ=='mom':
		for i in range(0,len(x),scale):
			#not divided by scale
			if i+scale<=len(x):
				val=np.std(x[i:i+scale])
				x_scale.append(val)
	
	elif cg_typ=='mavg_mom':
		for i in range(0,len(x)):
			#not divided by scale
			if i+scale<=len(x):
				val=np.std(x[i:i+scale])
				x_scale.append(val)	
	
	
	return np.array(x_scale)

def scale_series_emd(x,cg_typ):

	x_scale=[]
	emd = EMD()
	imf = emd.emd(x)
	m=imf.shape[0]
	
	if cg_typ=='c2f':
		for scl in range(1,m+1):
			val=np.sum(imf[:scl,:],axis=0)
			x_scale.append(val)	
			
	elif cg_typ=='f2c':
		imf=np.flip(imf,0)
		for scl in range(1,m+1):
			val=np.sum(imf[:scl],axis=0)
			x_scale.append(val)	
	return np.array(x_scale)

def get_scale_series(xt,scale=5,cg_type='base'):
	"""Returns a list of lists with different scaled RR series """
	Xs=[]
	
	if cg_type!='f2c' and cg_type!='c2f':
		Xs.append(xt)
		if cg_type=='comp':
			for s in range(2,scale+1):
				x_comp=[]
				for j in range(0,s):
					x_comp.append(scale_series(xt[j:],s,cg_typ='base'))
				Xs.append(x_comp)
		else:
			for s in range(2,scale+1):
				Xs.append(scale_series(xt,s,cg_type))
	elif cg_type=='f2c' or cg_type=='c2f':
		Xs=scale_series_emd(xt,cg_type)
	

	return Xs



def calc_mse(x_cg,ent_typ='mspe',typ='',mod_flag=0,deg=3,lag=1,m=2,r=0.15):
	"""
	Recieves a scaled series and returns multiscale entropy from it
	Possible options:
	mspe: multi-scale permutation entropy
	 - with modifications of modified PE, and weighted PE
	msse: multi-scale sample entropy
	comp_msse: composite scaled sample entropy
	comp_mspe: composite scaled permutation entropy
	"""
	ent=[]
	if ent_typ=='mspe':
		d_all= get_mot_ent_dist(x_cg,deg,lag,typ=typ,mod_flag=mod_flag)
		for dist in d_all:
			#dist will have distribution from given entropy type for all scales
			# Squence of values - PE:  PE1, PE2,...PEs, 
			ent.append(calc_mspe(dist))
				
	elif ent_typ=='msse':
		for x in x_cg:
			x=np.array(x).squeeze()
			x_std=np.std(x)
			rval=r*x_std
			ent.append(nolds.sampen(x,emb_dim=m,tolerance=rval))
			
	elif ent_typ=='comp_msse':
		for scl,x in enumerate(x_cg):
			cmp=[]
			if type(x) is np.ndarray:
				x=np.array(x).squeeze()
				x_std=np.std(x)
				rval=r*x_std
				cmp.append(nolds.sampen(x,emb_dim=m,tolerance=rval)/(scl+1))
			else:
				for x_cp in x:
					x_cp=np.array(x_cp).squeeze()
					x_std=np.std(x_cp)
					rval=r*x_std
					cmp.append(nolds.sampen(x_cp,emb_dim=m,tolerance=rval)/(scl+1))
			#naking sure infinite doesnt bother the set-up
			cmp=np.array(cmp)
			cmp=cmp[np.isfinite(cmp)]
			ent.append(np.sum(np.array(cmp),axis=0))  #protect against infinities
		
	
	elif ent_typ=='comp_mspe':
		for scl,x in enumerate(x_cg):
			cmp=[]
			if type(x) is np.ndarray:
				cmp.append(get_1s_pe(x,deg,lag,typ=typ ,mod_flag=mod_flag)/(scl+1))
			else:
				for x_cp in x:
					cmp.append(get_1s_pe(x_cp,deg,lag,typ=typ ,mod_flag=mod_flag)/(scl+1))
			
			ent.append(np.sum(np.array(cmp),axis=0))  #protect against infinities	
		ent=np.array(ent).T
	

	return np.array(ent)

def get_mse_fname(cg_type='base',ent_type='mspe',scale=6):
	
	nm=[]
	
	if cg_type=='base':
		cg='_scl_b'
	elif cg_type=='mov_avg':
		cg='_scl_ma'
	elif cg_type=='mom':
		cg='_scl_mom'
	elif cg_type=='mavg_mom':
		cg='_scl_mm'
	elif cg_type=='f2c':
		cg='_scl_fc'
	elif cg_type=='c2f':
		cg='_scl_cf'
		
	
	if ent_type=='msse':
		ent=['msse']
	elif ent_type=='comp_se':
		ent=['comp_se']
	elif ent_type=='mspe':
		ent=['mspe']
	elif ent_type=='comp_pe':
		ent=['comp_pe']
		
	scl=[]
	for s in range(1,scale+1):
		tmp='_s'+str(s)+cg
		scl.append(tmp)
	
	for en in ent:
		for sca in scl:
			nm.append(en+sca)

	
	return nm	
		

	
	
