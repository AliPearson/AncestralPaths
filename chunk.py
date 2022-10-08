import numpy as np

def chunk(ints, sample, cleaned=False):
	if(cleaned==True):
		chunks=ints
	else:
		chunks=np.transpose(ints[sample,:,:])
	rows_keep=[]
	row=0
	while row<(chunks.shape[0]-1):
		while(chunks[row, 2]==chunks[(row+1),2]):
			row+=1
			if(row==(chunks.shape[0]-1)):
				break
		rows_keep.append(row)
		row+=1
	chunks_rem=chunks[np.array(rows_keep),:]
	return(chunks_rem)
