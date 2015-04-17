from itertools import combinations
import csv
from mltemp import get_annotations

def combine_files():
	annots = get_annotations()

	ihead = open('header.csv','r')
	hread = csv.reader(ihead, delimiter=',',quotechar='"')
	head = hread.next()

	ihand = open('measures_20150128.csv','r')
	reader = csv.reader(ihand, delimiter=',', quotechar='"')
	ohand = open('measures_20150206.tab','w')
	writer = csv.writer(ohand, delimiter='\t')

	temp = [head[0]]
	temp.append('activity')
	temp.append('chelix')
	temp.append('dfg')
	temp += head[1:]
	writer.writerow(temp)

	ftype = ['c' for x in temp]
	ftype[0] = 's'
	ftype[1] = 'd'
	ftype[2] = 'd'
	ftype[3] = 'd'
	writer.writerow(ftype)

	fclass = ['' for x in temp]
	fclass[1] = 'class'
	writer.writerow(fclass)
	
	for row in reader:
		res = row[0]
		temp = [row[0]]
		temp.append(annots[res]['active'])
		temp.append(annots[res]['chelix'])
		temp.append(annots[res]['dfg'])
		temp += row[1:]
		writer.writerow(temp)
	return	

def make_header():
	temp = ['pdb_id']

	want = combinations(range(1,242),2)

	for a,b in want:
		temp.append('dist_%d_CA_%d_CA' % (a,b))

	for a in range(1,239):
		temp.append('pd_CA_%d_%d_%d_%d' % (a,a+1,a+2,a+3))

	for a in range(1,242):
		temp.append('phi_%d' % (a,))
		temp.append('psi_%d' % (a,))
		temp.append('chi1_%d' % (a,))

	outfile = open('header.csv','w')
	writer = csv.writer(outfile, delimiter=',', quotechar='"')
	writer.writerow(temp)
	outfile.close()

if __name__ == '__main__':
	res = combine_files()
