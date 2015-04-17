import csv
from progress.bar import Bar
from copy import copy
from itertools import combinations
from glob import glob
import MDAnalysis as md
from second_look import get_pseudokinase_annotations

#temp, configurable
datadir = '/share/Data/UPDATE_PDBS/ESBG_KINASE_PDBS_02-03-15/'
#datadir = '/home/dim/ESBG_KINASE_PDBS_12-09-14/'

def get_annotations(infile='annotated_full_dataset.txt'):
	head = ['family','chelix','dfg','a_loop','active']
	reader = csv.reader(open(infile,'r'),delimiter='\t')
	res = {}
	for row in reader:
		res[row[0]] = { x:y for x,y in zip(head, row[1:6]) }
	return res

def get_mappings():
	#temporary hack..  should import alignment (hmm or cma)
	nums = datadir+'PDB_NUM/*.num'
	nums = glob(nums)
	mapp = {}
	tnums = [x.split('/')[-1] for x in nums]
	for i,item in enumerate(tnums):
		pdb = item.split('.')[0].lower()+'_'+item.split('_')[1][0].lower()
		mapp.update({pdb:map_one(nums[i])})
	return mapp

def map_one(numfile):
	reader = csv.reader(open(numfile,'r'),delimiter=' ')
	mapp = {}
	for row in reader:
		if len(row) < 2:
			pass
			print row
		elif row[1].isdigit():
			#mapp[int(row[1])] = int(row[5].split('|')[0])
			mapp[int(row[1])] = int(row[-1])
			name = row[0].lower()+'_'+row[3].lower()
	return mapp

def get_concordant_list(infile='concordant_structs.txt'):
	res = []
	for line in open(infile,'r'):
		if len(line) > 0:
			res.append(line.strip())
	return res


def process_chain(infile, mapp):
	univ = md.Universe(infile)
	want = combinations(range(1,242),2)
	vals = []
	#pairwise distances
	for a,b in want:
		temp = '?'
		if a in mapp and b in mapp:
			ta = univ.selectAtoms("resid %d and name %s" % (mapp[a],'CA'))
			tb = univ.selectAtoms("resid %d and name %s" % (mapp[b],'CA'))
			if len(ta) == 1 and len(tb) == 1:
				tac = ta.coordinates()
				tbc = tb.coordinates()
				dist = tbc-tac
				tdist = [pow(x,2) for x in dist[0]]
				dist = sum(tdist)
				dist = pow(dist,0.5)
				temp = dist
		vals.append(copy(temp))

	#pseudo-dihedrals
	for a in range(1,239):
		temp = '?'
		if a in mapp and a+1 in mapp and a+2 in mapp and a+3 in mapp:
			ta = univ.selectAtoms("resid %d and name %s" % (mapp[a],'CA'))
			ta += univ.selectAtoms("resid %d and name %s" % (mapp[a+1],'CA'))
			ta += univ.selectAtoms("resid %d and name %s" % (mapp[a+2],'CA'))
			ta += univ.selectAtoms("resid %d and name %s" % (mapp[a+3],'CA'))
			if len(ta) == 4:
				temp = ta.dihedral()
		vals.append(copy(temp))

	#phi,psi,chi1
	for a in range(1,242):
		temp = ['?','?','?']
		#find a
		res = univ.residues
		tres = [tr.id for tr in res]
		if a in mapp and mapp[a] in tres:
			pos = tres.index(mapp[a])
			ta = univ.residues[pos]
			#phi
			tphi = ta.phi_selection()
			if tphi is not None and len(tphi) == 4:
				temp[0] = tphi.dihedral()
			#psi
			tpsi = ta.psi_selection()
			if tpsi is not None and len(tpsi) == 4:
				temp[1] = tpsi.dihedral()
			#chi1
			tchi1 = ta.chi1_selection()
			if tchi1 is not None and len(tchi1) == 4:
				temp[2] = tchi1.dihedral()
		vals += copy(temp)
	return vals

def get_have_list(infile='measures_20150127.csv'):
	reader = csv.reader(open(infile,'r'),delimiter=',')
	have = set()
	trow = reader.next()
	tlen = len(trow)
	for row in reader:
		if len(row) == len(trow):
			have.add(row[0])
	return have

def get_all_measurements(clist = get_concordant_list(), outfile = 'measures.csv'):
	mapp = get_mappings()

	outfile = open(outfile,'w')
	writer = csv.writer(outfile,delimiter=',')
	bar = Bar("Measuring pdbs.", max = len(clist)+1)
	for item in clist:
		if item not in mapp:
			print item,' not in mapp'
                a,b = item.split('_')
                ifile = a+'.pdb_'+b.upper()+'.pdb'
		vals = process_chain(datadir+'PDBS/'+ifile,mapp[item])	
		temp = [item]+vals
		writer.writerow(temp)
		bar.next()
	bar.finish()
	outfile.close()
	return

def get_measurements(flist, slist, ofile):
            
        return

def get_non_concordant_list():
        clist = get_concordant_list()
        cfiles = set([x[:4]+'.pdb_'+x[-1].upper()+'.pdb' for x in clist])
        temp = glob(datadir+'/PDB_NUM/*.num')
        pdbs = set([x.split('/')[-1] for x in temp])
        want = pdbs - cfiles
        want = [x[:4]+'_'+x[9].lower() for x in want]
        return list(want)
        
if __name__ == '__main__':
	clist=get_non_concordant_list()
        #get_measurements(clist,outfile='non_concordant_measures_20150417.csv')
