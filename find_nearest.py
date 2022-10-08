import numpy as np

def find_nearest(array,value):
	idx = np.searchsorted(array, value, side="left")
	if idx==array.shape[0]:
		return idx-1
	if (value - array[idx]) > 0:
		return idx+1
	else:
	     	return idx
