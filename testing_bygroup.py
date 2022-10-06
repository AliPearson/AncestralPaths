import msprime, sys, argparse, csv, itertools, math, random, numpy as np, tqdm, tskit, tensorflow as tf, keras, sys
from random import sample

import keras.backend.tensorflow_backend as tfback
from keras import backend
#import eli5
#from eli5.sklearn import PermutationImportance
#from sklearn.metrics import scorer

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
parser.add_argument("-model", help="Classifier model file to load", required=True)
parser.add_argument("-nn", help="The numder of nodes towards the root to traverse up the tree", required=False, default=4)
parser.add_argument("-poplab", help="Population labels file to load", required=True)
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
samples[0]=["present_day", 0, max(samples_inv[8])+1, 40000, 10000]
samples[1]=["Bronze_Age", min(samples_inv[0]), max(samples_inv[0])+1, 40000, 10000]
samples[2]=["BAA", min(samples_inv[1]),max(samples_inv[1])+1,20000,10000]
samples[3]=["Neo", min(samples_inv[3]), max(samples_inv[3])+1,20000,10000]
samples[4]=["Yam", min(samples_inv[2]), max(samples_inv[2])+1,20000,10000]
print(samples)


model = load_model(args.model)
model.summary()

cvscores=np.zeros((5,5))
num_seq=5
X_test_pool=np.zeros((50000, 1,5,9,5))
labels_test_pool=np.zeros((50000,5))
for ts_num,tseq in enumerate(range(5,10)):
	print("TREE SEQUENCE", tseq)
	#t_num=0
	counts=np.zeros((10000,5,(9*int(args.nn))), dtype=float)
	labels=np.zeros((10000,5), dtype=int)
	
	ts_rel=tskit.load(str(args.ts)+"_relate_popsize_"+str(tseq)+".trees")
	ts_sim=tskit.load(str(args.ts)+"_"+str(tseq)+".trees")
	num_trees=ts_rel.get_num_trees()
	num_sites=ts_rel.num_sites
	print("Number of sites=",num_sites, file=sys.stderr)
	print("Number of trees in tree sequence =", ts_rel.get_num_trees(), file=sys.stderr)
	for samset in range(5):
		t_num=0
		num_examples=(samples[samset][4]+10000) #Number of places we need to visit in each tseq
		num_samples=samples[samset][2]-samples[samset][1] 
		split=num_examples/num_samples #Number of places within each sample
		step=math.floor(num_sites/split) #The step in sites visited

		print(samples[samset][0], file=sys.stderr)
		print(num_samples, file=sys.stderr)
		print(num_examples, file=sys.stderr)
		print(step, file=sys.stderr)
		progress_bar=tqdm.tqdm(total=(samples[samset][4]))
		counter_list=[0]*6

		while sum(counter_list)<round(samples[samset][4]):	
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
							if counter_list[path-1]>=(25000):
								continue
						elif counter_list[path-1]>=(50000):
							continue
				
						counter_list[path-1]+=1
						labels[t_num, samset]=path
					
						############### Get relate GNNs
						parent_rel = tree_rel.parent(sample)
						v=0
						leaves_previous=set()
						while v != int(args.nn): #Traversing up to nn nodes 
							if (parent_rel==-1): ##Checking for root node
								break
					
							total={leaves for leaves in tree_rel.leaves(parent_rel) if leaves>=samples[0][2] and leaves not in leaves_previous}
							if (len(total)==0): #If no leaves are from the ancestral groups then move to the next node
								parent_rel=tree_rel.parent(parent_rel)
								continue
	
							##Finding the GNN distributions and the age of nodes
							age=tree_rel.time(parent_rel)
							age_norm=int(age)/1500
							counts[t_num][samset][8+(9*v)]+=age_norm

							add=float(1/len(total)) #Normalised so it is robust to sample size
							for leaf in total:
								leaves_previous.add(leaf)
								counts[t_num][samset][(9*v)+samples_pops[leaf]]+=add
							parent_rel=tree_rel.parent(parent_rel) #Traverse up the tree to the next node
							v+=1

						if (parent_rel==-1 and v!=int(args.nn)):
							while (v!=int(args.nn)):
								counts[t_num][samset][(9*v):(9+(9*v))]=-15
								v+=1
						
						########### Next trees
						t_num+=1
						progress_bar.update()
						if sum(counter_list)==round((samples[samset][4])) :
							break
					if sum(counter_list)==round((samples[samset][4])):
						break
				if sum(counter_list)==round((samples[samset][4])):
					print("breaking", file=sys.stderr)
					break
			step=step+(step/2)
		progress_bar.close()
		print("t_num=",t_num, file=sys.stderr)
		print("counter_list=", *counter_list, file=sys.stderr)
		print("sum counter list=", sum(counter_list), file=sys.stderr)


	################### Classifier training
	for samset in range(5):
		print("testing classifier",samples[samset][0], file=sys.stderr)
		X_train=counts[:,samset,:]
		data = X_train.reshape(X_train.shape[0], 1,5,9)
		lab_target=keras.utils.to_categorical(labels[:,samset]-1, 6)
	
		#if tseq==0:
		#	X_test_pool[:,:,:,:,samset]=data
		#	labels_test_pool=(labels-1)
		#else:
		X_test_pool[(ts_num*10000):((ts_num*10000)+10000),:,:,:,samset]=data
		labels_test_pool[(ts_num*10000):((ts_num*10000)+10000),samset]=labels[:,samset]-1
		#np.vstack((X_test_pool, data))
		#np.append(labels_test_pool, (labels-1))

		scores = model.evaluate(data, lab_target, verbose=0)
		print("%s: %.2f%%" % (model.metrics_names[1], scores[1]*100))
		classes=model.predict_classes(data, batch_size = 128)
		print(np.where((labels-1)<0))
		print(tf.math.confusion_matrix(classes, labels[:,samset]-1))
		#np.savetxt(str(args.out)+str(samples[samset][0])+"_confusion.txt", tf.math.confusion_matrix(classes, labels[:,samset]-1))
		cvscores[ts_num,samset]=(scores[1] * 100)

print(cvscores)
print(np.mean(cvscores, axis=1)[0])
logfile=open(str(args.out)+".kfold.log", 'w')
classes_pool=np.array([])
labels_pool=np.array([])
for i in range(5):
	print(i)
	print((np.mean(cvscores, axis=1)[i], np.std(cvscores,axis=1)[i]), file=logfile)
	classes=model.predict_classes(X_test_pool[:,:,:,:,i], batch_size = 128)
	classes_pool=np.concatenate((classes_pool, classes))
	confusion=tf.math.confusion_matrix(classes, labels_test_pool[:,i])
	labels_pool=np.concatenate((labels_pool, labels_test_pool[:,i]))
	np.savetxt(str(args.out)+"_"+str(samples[i][0])+"_confusion.txt",confusion)
np.savetxt(str(args.out)+".kfold.txt", np.array(cvscores))

print(classes_pool.shape)
X_test_pool=X_test_pool.reshape(X_test_pool.shape[0]*5, 1, 5, 9)
labels_test_pool=labels_test_pool.reshape(labels_test_pool.shape[0]*5)
classes=model.predict_classes(X_test_pool, batch_size = 128)
print(np.unique(classes))
print(X_test_pool.shape)

print(np.array_equal(classes, classes_pool))
print(np.array_equal(labels_test_pool, labels_pool))
confusion=tf.math.confusion_matrix(classes_pool, labels_pool)
print(confusion)
np.savetxt(str(args.out)+"_confusion.txt", confusion)

