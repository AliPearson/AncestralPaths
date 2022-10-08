
import numpy as np
from chunk import chunk
from find_nearest import find_nearest
from scipy.optimize import curve_fit
from scipy.stats.distributions import t
np.seterr(divide='ignore', invalid='ignore')
import random

def tracts_prob(ints, sample, cleaned=False):
	if cleaned:
		chunks=chunk(ints, sample, cleaned=True)
	else:
		chunks=chunk(ints, sample)
	results_array=np.zeros((50,6,6), dtype=float)
	for dist in range(1,51):
		#start_point=float(25-(0.5*dist))
		start_point=random.uniform(0,1)
		for point in np.arange(start_point, (chunks[-1,3]-dist), 1):
			row_x=find_nearest(chunks[:,3], point)
			path_x=int(chunks[row_x,2])-1
			row_y=find_nearest(chunks[:,3],	point+dist)
			path_y=int(chunks[row_y, 2])-1
			results_array[(dist-1), (path_x), (path_y)]+=1
	return results_array

def bronze_comb(res_chrom, path1, path2):
	path1=path1-1
	path2=path2-1
	#res_chrom[res_chrom==0]=np.nan
	res_comb=res_chrom[:,path1,:,:]+res_chrom[:,path2,:,:]
	print(np.where((res_comb[:,path1,:]+res_comb[:,path2,:])/np.sum(res_comb,1)==np.nan))
	std=np.nanvar((res_comb[:,path1,:]+res_comb[:,path2,:])/np.sum(res_comb,1),1)
	res=np.sum(res_chrom, 3)
	res_comb=res[:,path1,:]+res[:,path2,:]
	y=(res_comb[:,path1]+res_comb[:,path2])/np.nansum(res_comb,1)
	return y, std

def comb(res_chrom, path):
	path=path-1
	#res_chrom[res_chrom==0]=np.nan
	std=np.nanvar((res_chrom[:,path,path,:]/np.sum(res_chrom[:,path,:,:], 1)), 1)
	res=np.sum(res_chrom, 3)
	y=(res[:,path,path])/np.nansum(res[:,path,:],1)
	return y, std

def exp_fit(y, std):
	def func(x, alpha, beta, theta):
		return alpha*np.exp(-beta*x)+theta
	x=np.linspace(0.01,0.5,50)
	pars, pcov = curve_fit(func, x, y, sigma=1/std, absolute_sigma=True)
	a=0.05
	n=len(y)
	p=len(pars)
	dof=max(0,n-p)
	tval=t.ppf((1-a/2), dof)
	vars=np.diag(pcov)
	sigma=vars**0.5
	for p, var in zip(pars, sigma):
		print('p{0}: [{1} {2}]'.format(p, p-sigma*tval, p+sigma*tval)) 
	


