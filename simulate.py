import msprime, argparse, sys, math, random, demes, numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-out", help="Path and filename to save the tree sequences and vcfs", required=True)
parser.add_argument("-ages", help="File of population and ages for each diploid sample, tab separated", required=False)
args = parser.parse_args()

bronze_times=[]
baa_times=[]
yam_times=[]
neo_times=[]
whg_times=[]
ehg_times=[]
chg_times=[]
ana_times=[]
with open(args.ages, 'r') as sample_times:
	for line in sample_times:
		line=line.strip()
		field=line.split('\t')
		if field[0]=="Bronze_Age":
			bronze_times.append(int(field[1]))
		if field[0]=="BAA":
			baa_times.append(int(field[1]))
		if field[0]=="Yam":
			yam_times.append(int(field[1]))
		if field[0]=="Neo":
			neo_times.append(int(field[1]))
		if field[0]=="WHG":
			whg_times.append(int(field[1]))
		if field[0]=="EHG":
			ehg_times.append(int(field[1]))
		if field[0]=="CHG":
			chg_times.append(int(field[1]))
		if field[0]=="Ana":
			ana_times.append(int(field[1]))	

#initial population sizes:
N_bronze = 50000
N_Yam = 20000
N_baa = 20000
N_whg = 2000
N_ehg = 3000
N_ana = 5000
N_neo = 50000
N_chg = 3000
N_A = 15000 #Ancestor of WHG and EHG
N_B = 15000 #Ancestor of CHG and Neolithic farmers
N_present = 50000

#Time of events
T_bronze = max(bronze_times)+1
T_Yam = max(yam_times)+1
T_neo = max(neo_times)+1
T_baa = max(baa_times)+1
T_near_east = 800
T_europe = 600
T_basal = 1500

demography = msprime.Demography()
demography.add_population(name="present_bronze", initial_size=N_present) #0
demography.add_population(name="Yam", initial_size=N_Yam) #1
demography.add_population(name="Neo", initial_size=N_neo) #2
demography.add_population(name="EHG", initial_size=N_ehg) #3
demography.add_population(name="WHG", initial_size=N_whg) #4
demography.add_population(name="NE", initial_size=N_A) #5
demography.add_population(name="BAA", initial_size=N_baa) #6
demography.add_population(name="CHG", initial_size=N_chg) #7
demography.add_population(name="Ana", initial_size=N_ana) #8
demography.add_population(name="WA", initial_size=N_B) #9
demography.add_population(name="trunk", initial_size=N_B) #10

demography.add_admixture(time=T_bronze, derived="present_bronze", ancestral=["Neo", "Yam"], proportions=[0.5,0.5])
demography.add_admixture(time=T_neo, derived="Neo", ancestral=["WHG", "Ana"], proportions=[0.25,0.75])
demography.add_admixture(time=T_baa, derived="BAA", ancestral=["CHG", "Ana"], proportions=[0.25,0.75])
demography.add_admixture(time=T_Yam, derived="Yam", ancestral=["CHG", "EHG"], proportions=[0.5,0.5])
demography.add_population_split(time=T_europe, derived=["WHG", "EHG"], ancestral="NE")
demography.add_population_split(time=T_near_east, derived=["CHG", "Ana"], ancestral="WA")
demography.add_population_split(time=T_basal, derived=["WA", "NE"], ancestral="trunk")

demography.add_census(time=350)

present_sample = [msprime.SampleSet(91, time=0, population="present_bronze", ploidy=2)]
bronze_samples=[]
Baa_samples=[]
neolithic_samples=[]
Yam_samples=[]
CHG_samples=[]
Ana_samples=[]
WHG_samples=[]
EHG_samples=[]
for i in range(0,len(bronze_times)):
	time=bronze_times[i]
	bronze_samples.append(msprime.SampleSet(1, population="present_bronze", time=time, ploidy=2))
for time in baa_times:
	Baa_samples.append(msprime.SampleSet(1, time=time, population="BAA", ploidy=2))
for i in range(0, len(neo_times)):
	time=neo_times[i]
	neolithic_samples.append(msprime.SampleSet(1, time=time, population="Neo", ploidy=2))
for i in range(0,len(yam_times)):
	time=yam_times[i]
	Yam_samples.append(msprime.SampleSet(1, time=time, population="Yam", ploidy=2))
for i in range(0, len(chg_times)):
	time=chg_times[i]
	CHG_samples.append(msprime.SampleSet(1, time=time, population="CHG", ploidy=2))
for time in ana_times:
	Ana_samples.append(msprime.SampleSet(1, time=time, population="Ana", ploidy=2))
for time in whg_times:
	WHG_samples.append(msprime.SampleSet(1, time=time, population="WHG", ploidy=2))
for time in ehg_times:
	EHG_samples.append(msprime.SampleSet(1, time=time, population="EHG", ploidy=2))

samples = present_sample + bronze_samples + Baa_samples + neolithic_samples + Yam_samples + WHG_samples + EHG_samples + Ana_samples + CHG_samples  #Can multiply each by the number of samples needed for each population

demography.sort_events()
print(demography.debug())
graph=demography.to_demes()
demes.dump(graph, str(args.out)+".yaml")

#Simulate chromosome 3 only
tree_sequence_rep = msprime.sim_ancestry(num_replicates=10, recombination_rate=1e-8, demography=demography,samples=samples, sequence_length=2000000)


for j, tree_sequence in enumerate(tree_sequence_rep):
	ts=msprime.sim_mutations(tree_sequence, rate=1.25e-8)
	name=str(args.out)+"_"+str(j)+".trees"
	ts.dump(name)
 #Saving tree.sequence object to desired directory
	with open(str(args.out)+"_"+str(j)+".vcf", "w") as vcf_file:
		ts.write_vcf(vcf_file)

ages=tree_sequence.tables.nodes.time[0:tree_sequence.get_sample_size()]
np.savetxt(str(args.out)+".ages", ages)

poplab=open(str(args.out)+".poplabels", 'w')
nodes=tree_sequence.tables.nodes
pop_names=tree_sequence.tables.populations
print("sample", "population","group", "sex", sep=' ', file=poplab)
for i in range(int(tree_sequence.get_sample_size()/2)):
	if nodes.population[(i*2)]==0:
		if nodes.time[(i*2)]==0:
			print("tsk_"+str(i), "present_day", "present_day", "NA", sep=' ', file=poplab)
		else:
			print("tsk_"+str(i), "Bronze_Age", "Bronze_Age", "NA", sep=' ', file=poplab)
	else:
		print("tsk_"+str(i), pop_names[nodes.population[(i*2)]].metadata['name'], pop_names[nodes.population[(i*2)]].metadata['name'], "NA", sep=' ', file=poplab)
poplab.close()


