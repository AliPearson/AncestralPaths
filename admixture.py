from chunk import chunk
from tracts_prob import tracts_prob
import tqdm, argparse, numpy as np, pandas as pd, math
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri

from rpy2.robjects import r
from rpy2.robjects import numpy2ri
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
np. set_printoptions(threshold=np. inf)

r=robjects.r
r['source']('exp_fit.R')
exp_fit_r=robjects.globalenv['exp_fit']
r['source']('exp_fit_bronze.R')
exp_fit_bronze_r=robjects.globalenv['exp_fit']

parser = argparse.ArgumentParser()
parser.add_argument("-out", help="Filename to write admixture times and fractions to", required=True)
parser.add_argument("-poplab",help="Population labels file", required=True)
parser.add_argument("-samples",help="File containing order of samples by name", required=True)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-ints",nargs="+",help="List of intervals file names to load")
group.add_argument("-intsf",help="Name of file containing intervals file names")
args = parser.parse_args()

ints_mes={}
if args.ints==None:
	count=0
	with open(args.intsf, 'r') as ints_file:
		for line in ints_file:
			data=np.load(str(line))
			print(line)
			ints_mes[count] = data['arr_0']
			count+=1
else:
	count=0
	for num in args.ints:
		data=np.load(str(num))
		print(num)
		ints_mes[count] = data['arr_0']
		count+=1

print("data loaded", len(ints_mes), sep=' ')

sam_order=[]
with open(args.samples, 'r') as sams_ordered:
	for line in sams_ordered:
		line=line.strip()
		sam_order.append(line)

samples={}
with open(args.poplab, 'r') as poplab:
	i=0
	for line in poplab:
		if (line.startswith("sample")):
			continue
		line=line.strip()
		field=line.split(' ')
		if field[0] not in sam_order:
			continue
		index=sam_order.index(field[0])
		ind_1=(index*2)
		ind_2=(index*2)+1
		samples[int(ind_1)]=str(field[1])
		samples[int(ind_2)]=str(field[1])
		i+=1
print(len(samples))

num_samples=len(sam_order)
results_times=np.zeros((num_samples, 2,2), dtype=float)
results_fractions=np.zeros((num_samples,2,2), dtype=float)

progress_bar=tqdm.tqdm(total=num_samples)
for hap in range(0,(num_samples*2),2):
	if samples[hap]=="WHG" or samples[hap]=="EHG" or samples[hap]=="CHG" or samples[hap]=="Ana" or samples[hap]=="BAA" or samples[hap]=="Armenia" or samples[hap]=="FIN" or samples[hap]=="CEU" or samples[hap]=="GBR" or samples[hap]=="TSI" or samples[hap]=="IBS":
		results_times[int(hap/2),:,:]=np.nan
		results_fractions[int(hap/2),:,:]=np.nan
		progress_bar.update()
		continue
	res_chrom=np.zeros((50,6,6,len(ints_mes)), dtype=float)
	count=0
	for i in range(len(ints_mes)):
		ints=ints_mes[i]
		res_chrom[:,:,:,count]+=tracts_prob(ints, hap)
		res_chrom[:,:,:,count]+=tracts_prob(ints, (hap+1))
		count+=1

	if samples[hap]=="Bronze_Age":
		for ind, paths in enumerate([(1,3), (2,4)]):
			got_it=False
			for cut in range(0,25,2):
				res_chrom_try=res_chrom[np.arange((50-cut)), :,:,:]
				with localconverter(ro.default_converter + numpy2ri.converter):
					data.df = ro.conversion.py2rpy(res_chrom_try)
				try:
					results_r=exp_fit_bronze_r(data.df, paths[0], paths[1])
					got_it=True
					break
				except:
					continue
			if got_it:
				results=np.array(results_r)
				results_times[int(hap/2),ind,:]=results[0:2]
				results_fractions[int(hap/2),ind,:]=results[2:]
			else:
				results_times[int(hap/2),:,:]=np.nan
				results_fractions[int(hap/2),:,:]=np.nan

	elif samples[hap]=="Yam":
		for ind, path in enumerate([2,4]):
			got_it=False
			for cut in range(0,25,2):
				res_chrom_try=res_chrom[np.arange((50-cut)), :,:,:]
				with localconverter(ro.default_converter + numpy2ri.converter):
					data.df = ro.conversion.py2rpy(res_chrom_try)
				try:
					results_r=exp_fit_r(data.df, path)
					got_it=True
					break
				except:
					continue
			if got_it:
				results=np.array(results_r)
				results_times[int(hap/2),ind,:]=results[0:2]
				results_fractions[int(hap/2),ind,:]=results[2:]
			else:
				results_times[int(hap/2),:,:]=np.nan
				results_fractions[int(hap/2),:,:]=np.nan

	elif samples[hap]=="Neo":
		for ind, path in enumerate([1,(3,4)]):
			got_it=False
			for cut in range(0,25,2):
				res_chrom_try=res_chrom[np.arange((50-cut)), :,:,:]
				with localconverter(ro.default_converter + numpy2ri.converter):
					data.df = ro.conversion.py2rpy(res_chrom_try)
				try:
					if type(path)==int:
						results_r=exp_fit_r(data.df, path)
					else:
						results_r=exp_fit_bronze_r(data.df, path[0], path[1])
					got_it=True
					break
				except:
					continue
			if got_it:
				results=np.array(results_r)
				results_times[int(hap/2),ind,:]=results[0:2]
				results_fractions[int(hap/2),ind,:]=results[2:]
			else:
				results_times[int(hap/2),:,:]=np.nan
				results_fractions[int(hap/2),:,:]=np.nan
	progress_bar.update()
progress_bar.close()
np.savez_compressed(args.out+".times.npz", results_times)
np.savez_compressed(args.out+".fractions.npz", results_fractions)
