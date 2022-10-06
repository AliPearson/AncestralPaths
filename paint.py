import msprime, sys, argparse, csv, itertools, random, numpy as np, threading, tskit, math, tqdm
from random import sample

import argparse, tensorflow as tf, keras, pandas as pd, numpy as np, sys
import keras.backend.tensorflow_backend as tfback
from keras import backend
from numpy import loadtxt
from keras.models import load_model

def _get_available_gpus():
	"""Get a list of available gpu devices (formatted as strings).

	# Returns
		A list of available GPU devices.
	"""
	#global _LOCAL_DEVICES
	if tfback._LOCAL_DEVICES is None:
		devices = tf.config.list_logical_devices()
		tfback._LOCAL_DEVICES = [x.name for x in devices]
	return [x for x in tfback._LOCAL_DEVICES if 'device:gpu' in x.lower()]

tfback._get_available_gpus = _get_available_gpus
tfback._get_available_gpus()
tf.config.list_logical_devices()

backend.set_image_data_format('channels_first')


parser = argparse.ArgumentParser()
parser.add_argument("-out", help="Filename to write training matrix to ", required=True)
parser.add_argument("-model",help="name of the NN model",required=True)
parser.add_argument("-ts", help="Tree sequence file to load", required=True)
parser.add_argument("-poplab", help="File of population labels", required=True)
parser.add_argument("-nn", help="Number of nodes to traverse up to", required=True)
args = parser.parse_args()

#############################3
model=load_model(args.model)
ts=msprime.load(args.ts)
num_trees=ts.get_num_trees()
num_sites=ts.num_sites
num_muts=ts.get_num_mutations()
print("Number of sites=",num_sites)
print("Number of mutations=",num_muts)
print("Number of trees in tree sequence =", ts.get_num_trees())

#############################
def all_nodes(tree, samples, results):
	counts=np.zeros((num_samples,(9*int(args.nn))), dtype=np.float32)
	counted=set()

	for sam in range(num_samples):
		if sam in counted:
			continue
		sib=tree.left_sib(sam)
		if sib==-1:
			sib=tree.right_sib(sam)
		if sib<num_samples:
			pop=samples[sam]
			pop_sib=samples[sib]
			if pop_sib==pop:
				sam_rel=np.array([sam, sib])
			else:
				sam_rel=np.array([sam])
		else:
			child_pops=[samples[leaf] for leaf in tree.leaves(sib)]
			if len(set(child_pops))==1 and child_pops[0]==8 and samples[sam]==8:				
				sam_rel=np.append(np.array(list(tree.leaves(sib))), sam)
			else:
				sam_rel=np.array([sam])
		u=tree.parent(sam)
		nodes=[]
		v=0
		leaves_previous=set()
		while u != -1 and v != int(args.nn):
			total={leaves for leaves in tree.leaves(u) if samples[leaves]!=8 and leaves not in leaves_previous}
#### Moving to next node if none of the samples are found without recording age: Different to previous and how classifier is trained
			if len(total)==0:
				u=tree.parent(u)
				continue

##Finding the label from population and age of nodes
			age=tree.time(u)
			age_norm=float(age)/1500
			counts[sam_rel,8+(9*v)]=age_norm

				#if len(total)==0:
				#	continue
##Finding the GNN distributions and the age of nodes
			add=float(1/len(total))
			leaves_previous.update(total)
			for leaf in total:
				counts[sam_rel,(9*v)+samples[leaf]]+=add
			v+=1
			u=tree.parent(u)

		if (v!=int(args.nn)):
			while (v!=int(args.nn)):
				counts[sam,(9*v):(9+(9*v))]=-15
				v+=1
		counted.update(sam_rel)
##################################
	counts_test=counts.reshape(counts.shape[0],1,int(args.nn),9)
	results[t,:, 0]=model.predict_classes(counts_test, batch_size = 128)+1
	prediction=model.predict(counts_test, batch_size=128)
	results[t,:,1]=np.max(prediction, axis=1)

###################################
sam_order=[]
with open("sams_order.txt", 'r') as sams_ordered:
	for line in sams_ordered:
		line=line.strip()
		sam_order.append(line)
print(sam_order)
	

samples={}
with open(args.poplab, 'r') as poplab:
	i=0
	for line in poplab:
		if (line.startswith("sample")):
			continue
		line=line.strip()
		field=line.split(' ')
		index=sam_order.index(field[0])
		ind_1=(index*2)
		ind_2=(index*2)+1
		if str(field[1])=="Bronze_Age":		
			samples[int(ind_1)]=0
			samples[int(ind_2)]=0
		elif str(field[1])=="BAA":
			samples[int(ind_1)]=1
			samples[int(ind_2)]=1
		elif str(field[1])=="Yam":
			samples[int(ind_1)]=2
			samples[int(ind_2)]=2
		elif str(field[1])=="Neo":
			samples[int(ind_1)]=3
			samples[int(ind_2)]=3
		elif str(field[1])=="WHG":
			samples[int(ind_1)]=4
			samples[int(ind_2)]=4
		elif str(field[1])=="EHG":
			samples[int(ind_1)]=5
			samples[int(ind_2)]=5
		elif str(field[1])=="Ana":
			samples[int(ind_1)]=6
			samples[int(ind_2)]=6
		elif str(field[1])=="CHG":
			samples[int(ind_1)]=7
			samples[int(ind_2)]=7
		else:
			samples[int(ind_1)]=8
			samples[int(ind_2)]=8
		i+=1

print(len(samples), file=sys.stderr)
num_samples=len(list(ts.samples()))
print(num_samples)
##############################
results=np.zeros((num_trees,num_samples,3),dtype=float)
progress_bar = tqdm.tqdm(total=num_trees)
for t, tree in enumerate(ts.trees()):
	results[t,:,2]=tree.interval.right
	all_nodes(tree, samples, results)
	progress_bar.update()
np.savez_compressed(str(args.out)+"_painted.npz", paths=results)
progress_bar.close()

