#from sklearn.preprocessing import Imputer
from sklearn import preprocessing, cluster, decomposition, linear_model, svm, cross_validation, feature_selection
from sklearn.metrics import roc_curve, auc
from numpy.random import random_integers
from scipy import interp
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import csv
import pickle
from copy import copy
from itertools import combinations
from datetime import date
from progress.bar import Bar
import networkx as nx
from annotate import *

def get_data(infile='measures_20150128.csv'):
	temp,annots = [],[]
	reader = csv.reader(open(infile,'r'),delimiter=',',quotechar='"')
	head = reader.next()
	for row in reader:
		if len(row) > 0:
			trow = copy(row[1:])
			while '?' in trow:
				trow[trow.index('?')] = np.nan
			temp.append(trow)
			annots.append(row[0])
	pickle.dump([np.array(head),np.array(temp,dtype=float),np.array(annots)],open('data.p','w'))
	return np.array(head),np.array(temp,dtype=float),np.array(annots)

def get_annotations(infile='annotations.csv'):
        annots = []
        reader = csv.reader(open(infile,'r'),delimiter=',',quotechar='"')
        for row in reader:
                if len(row) > 0:
                        annots.append(row)
        return np.array(annots)

def get_pseudo_data(infile='pseudo_measures_20150403.csv'):
	temp,annots = [],[]
	reader = csv.reader(open(infile,'r'),delimiter=',',quotechar='"')
	#head = reader.next()
	for row in reader:
		if len(row) > 0:
                        temp.append(row[1:])
                        annots.append(row[0])
        #for item in temp:
            #while '?' in item:
                #item[item.index('?')] = np.nan
	nphead, nptemp, npannots = np.array(head),np.array(temp),np.array(annots)
        nptemp[nptemp == '?'] = np.nan
	return nphead,nptemp.astype(float),npannots

def preprocess(data, annots):
	#impute missing values
        tdata = copy(data)
        act = annots[:,0]
        act_data = data[act == 'active',:]
        inact_data = data[act == 'inactive',:]

	imp = preprocessing.Imputer(missing_values=np.nan, strategy='mean', axis=0)
        imp.fit_transform(data[act == 'active',:])
        imp.fit_transform(data[act == 'inactive',:])
	scaled_data = preprocessing.scale(data)
	return scaled_data, annots, tdata

def test(data,atype,head,feature_index=None):
	if feature_index is not None:
		for i,d in enumerate(feature_index):
			act = data[atype == 1,i]
			inact = data[atype == 0,i]
			plt.figure()
			n,b,p = plt.hist([act,inact],50,color=['red','blue'],label=['active','inactive'],stacked=True)
			plt.legend()
			plt.title(head[d])
			plt.tight_layout()
			plt.savefig('/home/dim/Dropbox/henrik/mltemp/%s.png' % (head[d],), dpi=300, format='png')
			plt.clf()
			#plt.show()

def evaluate_features(data, annot, head):
	#10 fold cross validation
	mean_tpr = 0.0
	mean_fpr = np.linspace(0,1,100)
	all_tpr = []
	classifier = svm.SVC(kernel='linear', probability=True)
	folds = 5
	skf = cross_validation.StratifiedKFold(annot, folds)
	#plt.figure()
	for train,test in skf:
		#print train,test,type(test)
		probs = classifier.fit(data[train],annot[train]).predict_proba(data[test])
		fpr, tpr, thresholds = roc_curve(annot[test], probs[:,1])
		#print '\t',auc(fpr,tpr)
		mean_tpr += interp(mean_fpr, fpr, tpr)
		mean_tpr[0] = 0.0
		#plt.plot(fpr, tpr, lw=1)
	mean_tpr /= folds 
	mean_tpr[-1] = 1.0
	print 'auc:',auc(mean_fpr, mean_tpr)
	#plt.plot([0,1],[0,1], '--', color=(0.5,0.5,0.5), label='Guess')
	#plt.plot(mean_fpr, mean_tpr, lw=4)
	#plt.xlim([-0.05,1.05])
	#plt.ylim([-0.05,1.05])
	#plt.xlabel('False Positive Rate')
	#plt.ylabel('True Positive Rate')
	#plt.show()
	return

def get_features(infile):
	res = []
	ihand = open(infile,'r')
	for line in ihand:
		if len(line) > 0:
			res.append(line.strip())
	return res

def get_resolutions(infile='resolu.idx'):
	ihandle = open(infile, 'r')
	res = {}
	for line in ihandle:
		temp = line.split(';')
		if len(temp) == 2:
			res[temp[0].strip().lower()] = float(temp[1].strip())
	return res

def get_pseudokinases(infile='pseudokinases.p'):
	return pickle.load(open(infile,'rU'))[0]

def get_pseudokinase_annotations(infile='pseudokinase_annotations.csv'):
        reader = csv.reader(open(infile,'r'),delimiter=',',quotechar='"')
        res = []
        #skip header
        reader.next()
        for row in reader:
                while len(row) < 5:
                    row.append(' ')
                res += [x for x in row]
        tres = np.array(res)
        res = tres.reshape((len(tres)/5,5))
        return res

def get_compounds(infile='compound.idx'):
	ihandle = open(infile, 'r')
	res = {}
	for line in ihandle:
		temp = line.split('\t')
		if len(temp) == 2:
			res[temp[0].lower()] = temp[1].strip()
	return res

def pk_join(pks, desc):
	ohand = open('pseudokinase_annotations.csv','w')
	writer = csv.writer(ohand, delimiter=',', quotechar='"')
	writer.writerow(['pdb_chn','VAIK','HRD','DFG','description'])
	for k,v in pks.items():
		tv = [1 if x != 0 else 0 for x in v]
		row = [k]+tv
		temp = k[:4]
		if temp in desc:
			row.append(desc[temp])
		writer.writerow(row)
	ohand.close()
	return

def make_pk_file():
	pks = get_pseudokinases()
	cmps = get_compounds()
	pk_join(pks, cmps)
	return
			

def eval_all(data, annot, head):
	clf = svm.SVC(kernel='linear')
	thead = list(head)
	'''
	res = []
	for a,b in combinations(thead,2):
		need = np.array([thead.index(a),thead.index(b)])
		scores = cross_validation.cross_val_score(clf, data[:,need], annot, cv=5)
		res.append((a,b,np.mean(scores),np.std(scores)))
	'''
	res = []
	have = ['dist_119_CA_140_CA','dist_141_CA_237_CA','dist_47_CA_117_CA']
	
	for i in thead:
		need = [thead.index(x) for x in have]
		need.append(thead.index(i))
		need.sort()
		temp = np.array(need)
		scores = cross_validation.cross_val_score(clf, data[:,temp], annot, cv=5)
		res.append((i,np.mean(scores),np.std(scores)))
	return res

def perform_pca(data,annots,title,legend_loc='upper right'):
	pca = decomposition.PCA(n_components=2)
	pdat = pca.fit(data).transform(data)
	f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, figsize=(14,6))
	tannot = annots[:,0]
	#print data.shape, tannot.shape
	for c,i in zip('rb', ['active','inactive']):
		ax1.scatter(pdat[tannot==i,0],pdat[tannot==i,1],c=c,label=i)

	tannot = annots[:,3]
	for c,i in zip('rbgy', ['helixC IN','helixC OUT', 'helixC DILATED', 'helixC DISORDER']):
		ax2.scatter(pdat[tannot==i,0],pdat[tannot==i,1],c=c,label=i)

	tannot = annots[:,4]
	for c,i in zip('rbg', ['dfg in','dfg out']):
		ax3.scatter(pdat[tannot==i,0],pdat[tannot==i,1],c=c,label=i)
	ax1.legend(loc=legend_loc)
	ax2.legend(loc=legend_loc)
	ax3.legend(loc=legend_loc)
	#plt.suptitle('PCA of %s' % title)
	plt.tight_layout()
	#plt.show()
	plt.savefig('./pics/'+title+'.svg', dpi=300, format='svg')
	return pdat 

def all_active(head,data,annots,annot_ind=0):
	#pdata = data
	#print 'all'	
	svc = svm.LinearSVC(C=0.001, penalty='l1', dual=False).fit(data,annots[:,annot_ind])
	svc_scores = zip(head[:-6],np.abs(svc.coef_[0]))
	svc_scores.sort(key=lambda kv: kv[1], reverse=True)
	logit = linear_model.LogisticRegression(penalty='l1',C=0.005).fit(data,annots[:,annot_ind])
	logit_scores = zip(head[:-6],np.abs(logit.coef_[0]))
	logit_scores.sort(key=lambda kv: kv[1], reverse=True)

	#print 'pca on all'
	#feats = get_features('all.feats')
	feats = [a for (a,b) in svc_scores if b > 0]
	feats += [a for (a,b) in logit_scores if b > 0]
	#featsi = [list(head).index(x) for x in feats]
	#featsi.sort()
	#allpca = perform_pca(data[:,featsi],annots,'Kinome')
	return list(feats)

def not_a_function():
	'''
	#identified features/classifiers
	clf = svm.SVC(kernel='linear')
	print len(featsi)
	print np.mean(cross_validation.cross_val_score(clf,pdata[:,featsi], annots[:,0], cv=5))
	print np.mean(cross_validation.cross_val_score(clf,tkdata[:,featsi], tkannot, cv=5))
	print np.mean(cross_validation.cross_val_score(clf,cmgcdata[:,featsi], cmgcannot, cv=5))
	print np.mean(cross_validation.cross_val_score(clf,agcdata[:,featsi], agcannot, cv=5))
	print np.mean(cross_validation.cross_val_score(clf,camkdata[:,featsi], camkannot, cv=5))
	print '\n'

	print len(tkfi)
	print np.mean(cross_validation.cross_val_score(clf,pdata[:,tkfi], annots[:,0], cv=5))
	print np.mean(cross_validation.cross_val_score(clf,tkdata[:,tkfi], tkannot, cv=5))
	print np.mean(cross_validation.cross_val_score(clf,cmgcdata[:,tkfi], cmgcannot, cv=5))
	print np.mean(cross_validation.cross_val_score(clf,agcdata[:,tkfi], agcannot, cv=5))
	print np.mean(cross_validation.cross_val_score(clf,camkdata[:,tkfi], camkannot, cv=5))
	print '\n'

	print len(cmgcfi)
	print np.mean(cross_validation.cross_val_score(clf,pdata[:,cmgcfi], annots[:,0], cv=5))
	print np.mean(cross_validation.cross_val_score(clf,tkdata[:,cmgcfi], tkannot, cv=5))
	print np.mean(cross_validation.cross_val_score(clf,cmgcdata[:,cmgcfi], cmgcannot, cv=5))
	print np.mean(cross_validation.cross_val_score(clf,agcdata[:,cmgcfi], agcannot, cv=5))
	print np.mean(cross_validation.cross_val_score(clf,camkdata[:,cmgcfi], camkannot, cv=5))
	print '\n'

	print len(agcfi)
	print np.mean(cross_validation.cross_val_score(clf,pdata[:,agcfi], annots[:,0], cv=5))
	print np.mean(cross_validation.cross_val_score(clf,tkdata[:,agcfi], tkannot, cv=5))
	print np.mean(cross_validation.cross_val_score(clf,cmgcdata[:,agcfi], cmgcannot, cv=5))
	print np.mean(cross_validation.cross_val_score(clf,agcdata[:,agcfi], agcannot, cv=5))
	print np.mean(cross_validation.cross_val_score(clf,camkdata[:,agcfi], camkannot, cv=5))
	print '\n'

	print len(camkfi)
	print np.mean(cross_validation.cross_val_score(clf,pdata[:,camkfi], annots[:,0], cv=5))
	print np.mean(cross_validation.cross_val_score(clf,tkdata[:,camkfi], tkannot, cv=5))
	print np.mean(cross_validation.cross_val_score(clf,cmgcdata[:,camkfi], cmgcannot, cv=5))
	print np.mean(cross_validation.cross_val_score(clf,agcdata[:,camkfi], agcannot, cv=5))
	print np.mean(cross_validation.cross_val_score(clf,camkdata[:,camkfi], camkannot, cv=5))
	print '\n'
	'''

def make_graph(data):
	from matplotlib.cm import jet, ScalarMappable
	from matplotlib.colors import Normalize
	g = nx.Graph()
	cnorm = Normalize(vmin=1, vmax=241)
	smap = ScalarMappable(norm=cnorm, cmap=jet)
	edge_list = []
	for k in data:
		tk = k.split('_')
		if len(tk) != 5:
			g.add_node(k)
		else:
			a,b = tk[1],tk[3]
			g.add_node(a)
			g.add_node(b)
			g.add_edge(a,b)
	pos = nx.spring_layout(g)
	nxcols,glabels = [],{}
	for i,node in enumerate(g.nodes()):
		if '_' not in node:
			nxcols.append(smap.to_rgba(int(node)))
		else:
			nxcols.append((0,1,0,0))
	nx.draw_networkx_nodes(g,pos,node_color=nxcols)
	nx.draw_networkx_labels(g,pos)
	nx.draw_networkx_edges(g,pos)
	plt.show()
	return 

def core_extended_features(head, data, annots, feats, title):
	legend_loc = 'upper right'
	tfeats = { x:y/2000.0 for x,y in feats.items()}
	core = [x for x,y in tfeats.items() if y == 1.0]
		
	ftold = None
	tannot = annots[:,0]
	tann = np.array([int(x) for x in tannot == 'active'])
	pdbid = annots[:,-1]
	all_data,all_ids = {},[]
	tt = { a:b for a,b in zip(pdbid,tannot)}
	#get data
	boundary = 0.75
	ft = {x:y for x,y in tfeats.items() if y >= boundary}
	ft = ft.items()
	ft.sort(key=lambda kv: float(kv[1]), reverse=True)
	a,b = zip(*ft)
	ext = list(a[len(core):])
	n = len(ext)+1
	
	def testpick(event):
		ind = event.ind
		print ind,np.take(pdbid,ind), np.take(tannot,ind)

	#plot figure
	if n > 1:
		f, ax = plt.subplots(1, n, sharey=True, figsize=(4*n,4))
		for t,axs in enumerate(ax):
			temp = core + ext[:t]
			tempi = [list(head).index(x) for x in temp]
			pca = decomposition.PCA(n_components=2)
			pdat = pca.fit_transform(data[:,tempi])
			print pdat[201],pdbid[201],tannot[201],tann[201]
			print pdat[119],pdbid[119],tannot[119],tann[119]
		
			clrs = np.array([ 'r' if x == 'active' else 'b' for x in tannot ])
			#for c,i in zip('rb', ['active','inactive']):
				#pdbs = pdbid[tannot==i]
			axs.scatter(pdat[:,0],pdat[:,1],c=clrs,picker=0.0001) 
			f.canvas.mpl_connect('pick_event', testpick)
			evaluate_features(data[:,tempi],tann,head)
			#af.append(AnnoteFinder(pdat[tannot==i,0],pdat[tannot==i,1],pdbid[tannot==i],axis=axs))
			#plt.connect('button_press_event', af[-1])
		print 'core:',core
		print 'ext:',ext
		#ax[0].legend(circs,loc=legend_loc)
		plt.suptitle('PCA of %s' % title)
		plt.tight_layout()
		plt.show()
		#plt.savefig('./pics/'+title+'.svg', dpi=300, format='svg')
	elif len(core) > 1:
		fig = plt.figure()	
		tempi = [list(head).index(x) for x in core]
		pca = decomposition.PCA(n_components=2)
		pdat = pca.fit(data[:,tempi]).transform(data[:,tempi])
		for c,i in zip('rb', ['active','inactive']):
			plt.scatter(pdat[tannot==i,0],pdat[tannot==i,1],c=c,label=i) 
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
	return tt 

def process_data():
	head,data,annots = pickle.load(open('data.p','r'))
	pdata,pannots = preprocess(data,annots)	
	pickle.dump([head,pdata,pannots],open('processed_data_final.p','w'))
        return

def find_all_features(head,pdata,pannots,res_bound=None,num=1000,outfile=None):
	#head,pdata,pannots = pickle.load(open('processed_data.p','r'))
	if res_bound is not None:
		res = get_resolutions()
		pres = np.array([res[x[-1][:-2]] for x in pannots])
		rdata = pdata[pres <= res_bound,:]
		rannots = pannots[pres <= res_bound,:]
	else:
		rdata = pdata
		rannots = pannots
	#num = 10
	#wanted_groups = ['TK','CMGC','AGC','CAMK']
	wanted_groups = ['TK','CMGC']
	group = rannots[:,2]
	result = {}
	#each group individually
	for grp in wanted_groups:
		res = []
		#print pdata.shape, group.shape
		tdata = rdata[group == grp,:]
		tannots = rannots[group == grp,:]
		act = tannots[:,0]
		tact = np.count_nonzero(act == 'active')
		print tdata.shape, tannots.shape, '%d active, %d inactive' % (tact,len(act)-tact)
		if tact < len(act):
			bar = Bar(grp, max=num)
			for i in range(num):
				feats = all_active(head,tdata,tannots,0)
				res += list(feats)
				bar.next()
			bar.finish()
			temp = [ (a,res.count(a)) for a in set(res) ]
			temp.sort(key = lambda kv: kv[1], reverse=True)
			print len(temp)
			result[grp] = temp
	#entire kinome
	res = []
	act = rannots[:,0]
	tact = np.count_nonzero(act == 'active')
	print rdata.shape, rannots.shape, '%d active, %d inactive' % (tact,len(act)-tact)
	grp = 'Kinome'
	bar = Bar(grp, max=num)
	for i in range(num):
		feats = all_active(head,rdata,rannots,0)
		res += list(feats)
		bar.next()
	bar.finish()
	temp = [ (a,res.count(a)) for a in set(res) ]
	temp.sort(key = lambda kv: kv[1], reverse=True)
	result[grp] = temp
	#just TK & CMGC
	grp = 'TK/CMGC'
	tk = group == 'TK'
	cmgc = group == 'CMGC'
	both = np.array([a or b for a,b in zip(tk,cmgc)])
	tdata = rdata[both,:]
	tannots = rannots[both,:]
	res = []
	act = tannots[:,0]
	tact = np.count_nonzero(act == 'active')
	print tdata.shape, tannots.shape, '%d active, %d inactive' % (tact,len(act)-tact)
	bar = Bar(grp, max=num)
	for i in range(num):
		feats = all_active(head,tdata,tannots,0)
		res += list(feats)
		bar.next()
	bar.finish()
	temp = [ (a,res.count(a)) for a in set(res) ]
	temp.sort(key = lambda kv: kv[1], reverse=True)
	result[grp] = temp

	#save features for later
	today = date.today()
	savedate = ''.join(today.isoformat().split('-'))
	if res_bound is not None and outfile is None:
		pickle.dump(result,open('all_features_%s_%s.p' % (savedate, res_bound),'w'))
	elif outfile is None:
		pickle.dump(result,open('all_features_%s.p' % (savedate, ), 'w'))
	else:
		pickle.dump(result,open(outfile,'w'))
	return result

if __name__ == '__main__':
	#features = find_all_features()
        pseudo = get_pseudokinase_annotations()
	head,pdata,pannots = pickle.load(open('processed_data.p','r'))

	head = head[:-6]
	group = pannots[:,2]
        pdb = pannots[:,-1]
        #remove pseudokinases
        pseudo_i = [list(pdb).index(x) for x in pseudo[:,0] if x in pdb]
        a,b = pdata.shape
        temp = range(a)
        not_pseudo = [x for x in temp if x not in pseudo_i]
        tdata = pdata[not_pseudo,:]
        tannots = pannots[not_pseudo,:]
        test_data = pdata[pseudo_i,:]
        test_annots = pannots[pseudo_i,:]

        test_head,test_data,test_annots = get_pseudo_data()

        #res = find_all_features(head,tdata,tannots,num=50,outfile='no_pseudo_training_features.p')
        #res = pickle.load(open('no_pseudo_training_features.p','r'))
        res = pickle.load(open('all_features_v2.p','r'))
        for k,v in res.items():
            print k
	    #svc = svm.LinearSVC(dual=False).fit(tdata,tannots[:,0])
            feats = { x:y for x,y in res[k] }
	    tfeats = { x:y/100.0 for x,y in feats.items()}
	    core = [x for x,y in tfeats.items() if y == 1.0]
            core_i = [list(head).index(x) for x in core]
            ttdata = tdata[:,core_i]
	    svc = svm.LinearSVC(dual=False).fit(ttdata,tannots[:,0])
            ttest_data = test_data[:,core_i]
            rest = svc.predict(ttest_data)
            print zip(test_annots,rest)

	#features = pickle.load(open('all_features_2.p','r'))[0]
	#features = pickle.load(open('all_features_20150331_2.0.p','r'))[0]
	#core_extended_features(head, pdata[group=='TK'], pannots[group=='TK'], dict(features['TK']), 'TK')
	#res = core_extended_features(head,pdata,pannots,dict(features['Kinome']),'Kinome')
	'''	
	tk = group == 'TK'
	cmgc = group == 'CMGC'
	both = np.array([a or b for a,b in zip(tk,cmgc)])
	core_extended_features(head,pdata[both,:],pannots[both,:],dict(features['TK/CMGC']),'TK/CMGC')
	for grp in ['TK','CMGC']:
		core_extended_features(head, pdata[group==grp], pannots[group==grp], dict(features[grp]), grp)

	#for title in features.keys():
		#core_extended_features(head, data, annots, features, title)
	'''
