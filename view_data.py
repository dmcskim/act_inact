#from sklearn.preprocessing import Imputer
from sklearn import decomposition, svm, cross_validation
from sklearn.metrics import roc_curve, auc
from scipy import interp
import matplotlib.pyplot as plt
import numpy as np
import pickle

def evaluate_features(data, annot, head):
	mean_tpr = 0.0
	mean_fpr = np.linspace(0,1,100)
	all_tpr = []
	classifier = svm.SVC(kernel='linear', probability=True)
	folds = 5
	skf = cross_validation.StratifiedKFold(annot, folds)
	for train,test in skf:
		probs = classifier.fit(data[train],annot[train]).predict_proba(data[test])
		fpr, tpr, thresholds = roc_curve(annot[test], probs[:,1])
		mean_tpr += interp(mean_fpr, fpr, tpr)
		mean_tpr[0] = 0.0
	mean_tpr /= folds 
	mean_tpr[-1] = 1.0
	print 'auc:',auc(mean_fpr, mean_tpr)
	return

def get_resolutions(infile='resolu.idx'):
	ihandle = open(infile, 'r')
	res = {}
	for line in ihandle:
		temp = line.split(';')
		if len(temp) == 2:
			res[temp[0].strip().lower()] = float(temp[1].strip())
	return res

def core_extended_features(head, data, annots, feats, title, boundary):
	legend_loc = 'upper right'
	div = float(max(feats.values()))
	tfeats = { x:y/div for x,y in feats.items()}
	core = [x for x,y in tfeats.items() if y == 1.0]

	for x in core:
		del tfeats[x]
		
	ftold = None
	tannot = annots[:,0]
	tann = np.array([int(x) for x in tannot == 'active'])
	pdbid = annots[:,-1]
	all_data,all_ids = {},[]
	#get data
	ft = {x:y for x,y in tfeats.items() if y >= boundary}
	ext = []
	if len(ft) > 0:
		ft = ft.items()
		ft.sort(key=lambda kv: float(kv[1]), reverse=True)
		a,b = zip(*ft)
		ext = list(a)
	n = len(ext)+1
	
	def testpick(event):
		ind = event.ind
		print np.take(pdbid,ind), np.take(tannot,ind)

	#plot figure
	if n > 1:
		f, ax = plt.subplots(1, n, sharey=True, figsize=(4*n,4))
		for t,axs in enumerate(ax):
			temp = core + ext[:t]
			tempi = [list(head).index(x) for x in temp]
			clrs = np.array([ 'r' if x == 'active' else 'b' for x in tannot ])
			pca = decomposition.PCA(n_components=2)
			pdat = pca.fit(data[:,tempi]).transform(data[:,tempi])
			#for c,i in zip('rb', ['active','inactive']):
			axs.scatter(pdat[:,0],pdat[:,1],c=clrs,picker=0.0001) 
			f.canvas.mpl_connect('pick_event', testpick)
			evaluate_features(data[:,tempi],tann,head)
			#af.append(AnnoteFinder(pdat[tannot==i,0],pdat[tannot==i,1],pdbid[tannot==i],axis=axs))
			#plt.connect('button_press_event', af[-1])
		print 'core:',core
		print 'ext:',ext
		ax[0].legend(loc=legend_loc)
		plt.suptitle('PCA of %s' % title)
		plt.tight_layout()
		plt.show()
		#plt.savefig('./pics/'+title+'.svg', dpi=300, format='svg')
	elif len(core) > 1:
		f = plt.figure()	
		tempi = [list(head).index(x) for x in core]
		pca = decomposition.PCA(n_components=2)
		pdat = pca.fit(data[:,tempi]).transform(data[:,tempi])
		clrs = np.array([ 'r' if x == 'active' else 'b' for x in tannot ])
		plt.scatter(pdat[:,0],pdat[:,1],c=clrs,picker=0.0001)
		f.canvas.mpl_connect('pick_event', testpick)
		evaluate_features(data[:,tempi],tann,head)
		plt.legend(loc=legend_loc)
		plt.title('PCA of %s' % title)
		plt.tight_layout()
		plt.show()
	else:
		ind = list(head).index(core[0])	
		act = data[tannot == 'active',ind]
		inact = data[tannot== 'inactive',ind]
		plt.figure()
		n,b,p = plt.hist([act,inact],50,color=['red','blue'],label=['active','inactive'],stacked=True)
		plt.legend()
		plt.title(core[0])
		plt.tight_layout()
		plt.show()
	return 

if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("input_file", help="input data file", type=str)
	parser.add_argument("feature_file", help="input feature file", type=str)
	parser.add_argument("-r", "--resolution", help="structure resolution limit", type=float)
	parser.add_argument("-b", "--boundary", help="relative frequency for features to be selected", type=float, default=0.75)
	args = parser.parse_args()

	#head,pdata,pannots = pickle.load(open('pca_data.p','r'))
	head,pdata,pannots = pickle.load(open(args.input_file,'r'))

	if args.resolution != None:
		res = get_resolutions()
		pres = np.array([res[x[-1][:-2]] for x in pannots])
		rdata = pdata[pres <= args.resolution,:]
		rannots = pannots[pres <= args.resolution,:]
	else:
		rdata = pdata
		rannots = pannots
	#features = find_all_features()

	group = rannots[:,2]
	#features = pickle.load(open('all_features_2.p','r'))[0]
	features = pickle.load(open(args.feature_file,'r'))
	#core_extended_features(head, pdata[group=='TK'], pannots[group=='TK'], dict(features['TK']), 'TK')
	print rdata.shape
	core_extended_features(head,rdata,rannots,dict(features['Kinome']),'Kinome', boundary=args.boundary)
	print '\n\n'
	
	tk = group == 'TK'
	cmgc = group == 'CMGC'
	both = np.array([a or b for a,b in zip(tk,cmgc)])
	print rdata[both,:].shape
	core_extended_features(head,rdata[both,:],rannots[both,:],dict(features['TK/CMGC']),'TK/CMGC', boundary=args.boundary)
	print '\n\n'
	for grp in ['TK','CMGC']:
		print rdata[group==grp].shape
		if len(dict(features[grp])) > 0:
			core_extended_features(head, rdata[group==grp], rannots[group==grp], dict(features[grp]), grp, boundary=args.boundary)
		else:
			print 'No features found for %s.' % (grp,)
		print '\n\n'

