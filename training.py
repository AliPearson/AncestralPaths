import msprime, sys, argparse, csv, itertools, math, random, numpy as np, tqdm, tskit, tensorflow as tf, keras, sys
from random import sample

import keras.backend.tensorflow_backend as tfback
from keras import backend
import sklearn

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
from keras.models import Sequential, load_model
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Convolution2D, MaxPooling2D

parser = argparse.ArgumentParser()
parser.add_argument("-out", help="Filename to write model to", required=True)
parser.add_argument("-ts", help="Tree sequences prefix to load", required=True)
parser.add_argument("-poplab", help="Population labels file to load", required=True)
parser.add_argument("-nn", help="The numder of nodes towards the root to traverse up the tree", required=False, default=5)
args = parser.parse_args()

samples_pops={}
with open(args.poplab, 'r') as poplab:
	i=0
	for line in poplab:
		if (line.startswith("sample")):
			continue
		line=line.strip()
		field=line.split(' ')
		ind_1=(i*2)
		ind_2=ind_1+1
		if str(field[1])=="Bronze_Age":
			samples_pops[int(ind_1)]=0
			samples_pops[int(ind_2)]=0
		elif str(field[1])=="BAA":
			samples_pops[int(ind_1)]=1
			samples_pops[int(ind_2)]=1
		elif str(field[1])=="Yam":
			samples_pops[int(ind_1)]=2
			samples_pops[int(ind_2)]=2
		elif str(field[1])=="Neo":
			samples_pops[int(ind_1)]=3
			samples_pops[int(ind_2)]=3
		elif str(field[1])=="WHG":	
			samples_pops[int(ind_1)]=4
			samples_pops[int(ind_2)]=4
		elif str(field[1])=="EHG":
			samples_pops[int(ind_1)]=5
			samples_pops[int(ind_2)]=5
		elif str(field[1])=="Ana":
			samples_pops[int(ind_1)]=6
			samples_pops[int(ind_2)]=6
		elif str(field[1])=="CHG":
			samples_pops[int(ind_1)]=7
			samples_pops[int(ind_2)]=7
		else:
			samples_pops[int(ind_1)]=8
			samples_pops[int(ind_2)]=8
		i+=1

samples_inv={}
for k, v in samples_pops.items():
	samples_inv[v] = samples_inv.get(v, []) + [k]
	
num_seq=5
samples={}
samples[0]=["present_day", 0, max(samples_inv[8])+1, 100000*num_seq, 50000*num_seq]
samples[1]=["Bronze_Age", min(samples_inv[0]), max(samples_inv[0])+1,100000*num_seq, 50000*num_seq]
samples[2]=["BAA", min(samples_inv[1]),max(samples_inv[1])+1,25000*num_seq, 12500*num_seq]
samples[3]=["Neo", min(samples_inv[3]), max(samples_inv[3])+1,50000*num_seq, 25000*num_seq]
samples[4]=["Yam", min(samples_inv[2]), max(samples_inv[2])+1,25000*num_seq, 25000*num_seq]

counts=np.zeros((162500*num_seq,(9*int(args.nn))), dtype=float)
labels=np.zeros((162500*num_seq), dtype=int)

t_num=0
for tseq in range(5):
	ts_rel=tskit.load(str(args.ts)+"_relate_popsize_"+str(tseq)+".trees")
	ts_sim=tskit.load(str(args.ts)+"_"+str(tseq)+".trees")
	num_trees=ts_rel.get_num_trees()
	num_sites=ts_rel.num_sites
	print("Number of sites=",num_sites, file=sys.stderr)
	print("Number of trees in tree sequence =", ts_rel.get_num_trees(), file=sys.stderr)
	for samset in range(5):
		num_examples=(samples[samset][3]+10000)/num_seq #Number of places we need to visit in each tseq
		num_samples=samples[samset][2]-samples[samset][1] 
		split=num_examples/num_samples #Number of places within each sample
		step=math.floor(num_sites/split) #The step in sites visited

		print(samples[samset][0], file=sys.stderr)
		print(num_samples, file=sys.stderr)
		print(num_examples, file=sys.stderr)
		print(step, file=sys.stderr)
		progress_bar=tqdm.tqdm(total=round((samples[samset][4])/num_seq))
		counter_list=[0]*6

		while sum(counter_list)<round((samples[samset][4])/num_seq):	
			print("step=", step, file=sys.stderr)
			tree_sim=ts_sim.first()
			for tree_rel in ts_rel.trees():
				for site in tree_rel.sites():
					if site.id % step != 0:
						continue
					pos=site.position
					while tree_sim.interval.right <= pos:
						tree_sim.next()
					for sample in range(samples[samset][1], samples[samset][2]):
						############# Get label from sim
						parent_sim = tree_sim.parent(sample)
						while int(tree_sim.time(parent_sim))!=350:
							parent_sim=tree_sim.parent(parent_sim)
						pop=tree_sim.get_population(parent_sim)
						if samples[samset][0]=="BAA":
							if pop==8:
								path=5
							elif pop==7:
								path=6
						elif pop==3:
							path=4
						elif pop==4:
							path=3
						elif pop==7:
							path=2
						elif pop==8:
							path=1
						
						if samples[samset][0]=="BAA":
							if counter_list[path-1]>=(6250):
								continue
						elif counter_list[path-1]>=(12500):
							continue
				
						counter_list[path-1]+=1
						labels[t_num]=path
					
						############### Get relate GNNs
						parent_rel = tree_rel.parent(sample)
						v=0
						leaves_previous=set()
						while v != int(args.nn): #Traversing up to nn nodes 
							if (parent_rel==-1): ##Checking for root node
								break
					
							total={leaves for leaves in tree_rel.leaves(parent_rel) if leaves>=int(samples[0][2]) and leaves not in leaves_previous}
							if (len(total)==0): #If no leaves are from the ancestral groups then move to the next node
								parent_rel=tree_rel.parent(parent_rel)
								continue
	
							##Finding the GNN distributions and the age of nodes
							age=tree_rel.time(parent_rel)
							age_norm=int(age)/1500
							counts[t_num][8+(9*v)]+=age_norm

							add=float(1/len(total)) #Normalised so it is robust to sample size
							for leaf in total:
								leaves_previous.add(leaf)
								counts[t_num][(9*v)+samples_pops[leaf]]+=add
							parent_rel=tree_rel.parent(parent_rel) #Traverse up the tree to the next node
							v+=1

						if (parent_rel==-1 and v!=int(args.nn)):
							while (v!=int(args.nn)):
								counts[t_num][(9*v):(9+(9*v))]=-15
								v+=1
						
						########### Next trees
						t_num+=1
						progress_bar.update()
						if sum(counter_list)>=round((samples[samset][4])/num_seq) :
							break
					if sum(counter_list)>=round((samples[samset][4])/num_seq):
						break
				if sum(counter_list)>=round((samples[samset][4])/num_seq):
					print("breaking", file=sys.stderr)
					break
			step=step+(step/2)
		progress_bar.close()
		print("t_num=",t_num, file=sys.stderr)
		print("counter_list=", *counter_list, file=sys.stderr)
		print("sum counter list=", sum(counter_list), file=sys.stderr)


################### Classifier training
print("training classifier", file=sys.stderr)
X_train=counts
X_train = X_train.reshape(X_train.shape[0], 1,5,9)
print(X_train.shape, file=sys.stderr)
print(labels.shape, file=sys.stderr)

lab_train_target=keras.utils.to_categorical(labels-1, 6)
model = Sequential()
model.add(Convolution2D(32, (5,1), activation='relu', input_shape=(1,5,9)))
model.add(Flatten())
model.add(Dense(units=64, activation='relu'))
model.add(Dense(units=6, activation='softmax'))
model.summary()
model.compile(loss='categorical_crossentropy',optimizer='adam',metrics=['accuracy'])
model.fit(X_train, lab_train_target ,batch_size=30, nb_epoch=5, verbose=1)
name="model_"+str(args.out)+str(".h5")
model.save(name)

