from CGE_GmsPython import *
import GamsProduction

class Production(GmsPython):
	def __init__(self, f=None, tree=None, ns=None, s=None, glob=None, s_kwargs = None, g_kwargs = None, dur_kwargs = None):
		""" Initialize from a pickle file 'f' or nesting tree 'tree'. """
		super().__init__(name=tree.name if tree else None, f=f, s=s, glob=glob, ns=ns, s_kwargs = s_kwargs, g_kwargs=g_kwargs)
		if f is None:
			self.readTree(tree)
			self.addDurables(**noneInit(dur_kwargs, {}))
			self.addPriceWedge()

	def addDurables(self, dur = None, f = None, dur2inv = None):
		self.ns.update({k: f"{k}_{self.name}" for k in ('dur','inv')})
		self.s.db[self.n('dur')] = noneInit(dur, pd.MultiIndex.from_product([self.get('output').levels[0], []], names = ['s','n']))
		self.s.db[self.n('dur2inv')] = pd.MultiIndex.from_frame(self.get('dur').to_frame(index=False).assign(nn=lambda x: 'I_'+x.n)) if dur2inv is None else dur2inv
		self.s.db[self.n('inv')] = self.get('dur2inv').droplevel(self.ns['n']).unique().rename({self.ns['nn']:self.ns['n']})
		self.s.db['n'] = self.get('n').union(self.get('inv').levels[-1])
		self.s.db[self.n('input')] = self.get('input').difference(self.get('dur')).union(self.get('inv'))
		self.m[f"{self.name}_IC"] = Submodule(**{'f': noneInit(f, 'sqrAdjCosts')})

	def addPriceWedge(self,f=None):
		self.ns.update({k: f"{k}_{self.name}" for k in ['s','output_n']})
		self.s.db[self.n('output_n')] = self.get('output').levels[-1]
		self.s.db[self.n('s')] = self.get('output').levels[0]
		self.m[f"{self.name}_pWedge"] = Submodule(**{'f': noneInit(f,'priceWedge')})

	def readTree(self,tree):
		gpyDB_wheels.robust.robust_merge_dbs(self.s.db,tree.db,priority='second')
		self.ns.update(tree.ns)
		[self.readTree_i(t) for t in tree.trees.values()];
		self.addCalibrationSubsets(tree)

	def readTree_i(self,t):
		self.addModule(t,**{k:v for k,v in t.__dict__.items() if k in ('name','ns','f','io','sp')})

	def addCalibrationSubsets(self,tree):
		self.ns.update({k: k+'_'+self.name for k in ('exomu','endo_qS','endo_qD','endo_pS')})
		exomu_inp, exomu_out = self.getExomuFromTree(tree)
		self.s.db[self.ns['exomu']] = exomu_inp.union(exomu_out)
		self.s.db[self.ns['endo_qD']] = exomu_inp.droplevel('n').rename(['s','n']).union(gpyDB_wheels.adj.rc_pd(exomu_out.droplevel('nn'),('not',self.get('output'))))
		self.s.db[self.ns['endo_qS']] = gpyDB_wheels.adj.rc_pd(exomu_out.droplevel('nn'),self.get('output'))
		self.s.db[self.ns['endo_pS']] = self.uniqueFromMap(self.get('output'),gb=['s'])

	def TroubleNodes(self,map_spinp, map_spout, output=None):
		""" Identify nodes that are branches in input tree and output tree"""
		return gpyDB_wheels.adj.rc_pd(map_spinp.droplevel('n').rename(['s','n']), c = gpyDB_wheels.adj.rc_pd(map_spout.droplevel('nn'), c = None if output is None else ('not', output)))

	def uniqueFromMap(self,map_,gb=('s','n')):
		""" MultiIndex-like groupby statement with function 'first' """
		return pd.MultiIndex.from_frame(map_.to_frame(index=False).groupby([s for s in gb]).first().reset_index()).reorder_levels(map_.names)

	def getCleanExoMu(self,map_spinp, map_spout, trouble):
		""" return elements to be added to exomu_inp, exomu_out"""
		return self.uniqueFromMap(gpyDB_wheels.adj.rc_pd(map_spinp, c= ('not', trouble.rename(['s','nn'])))), self.uniqueFromMap(gpyDB_wheels.adj.rc_pd(map_spout, c= ('not', trouble)),gb=('s','nn'))

	def getExomuFromTree(self,tree,maxiter=10):
		map_spinp = self.get('map_spinp').copy()
		map_spout = self.get('map_spout').copy()
		trouble = self.TroubleNodes(map_spinp, map_spout,output=self.get('output'))
		exomu_inp, exomu_out = self.getCleanExoMu(map_spinp, map_spout, trouble)
		i = 0
		while not trouble.empty:
			map_spinp = gpyDB_wheels.adj.rc_pd(map_spinp, ('not', exomu_inp.droplevel('nn')))
			map_spout = gpyDB_wheels.adj.rc_pd(map_spout, ('not', exomu_out.droplevel('n')))
			trouble = self.TroubleNodes(map_spinp, map_spout)
			exomu_inp_i, exomu_out_i = self.getCleanExoMu(map_spinp,map_spout,trouble)
			exomu_inp, exomu_out = exomu_inp.union(exomu_inp_i), exomu_out.union(exomu_out_i)
			i += 1
			if i == maxiter:
				print("Algorithm for exomu sets failed; increase maxiter or consider nesting for loops.")
				break
		return exomu_inp, exomu_out

	def initDB(self,m=None):
		return gpyDB_wheels.robust.robust_merge_dbs(self.s.db,self.initSymbols(m=m),priority='first')

	def initSymbols(self,m=None):
		return {self.n('pS'): gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('output',m=m)],self.s.db), name=self.n('pS'))),
				self.n('qS'): gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('output',m=m)],self.s.db), name=self.n('qS'))),
				self.n('pD'): gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('int',m=m).union(self.get('input',m=m)).union(self.get('dur',m=m))],self.s.db), name = self.n('pD'))),
				self.n('qD'): gpy(pd.Series(0.5, index = adjMultiIndexDB.mergeDomains([self.get('t'),self.get('int',m=m).union(self.get('input',m=m)).union(self.get('dur',m=m))],self.s.db), name = self.n('qD'))),
				self.n('mu'): gpy(pd.Series(1, index = self.get('map',m=m), name = self.n('mu'))),
				self.n('eta'): gpy(pd.Series(0.5, index = self.get('knout',m=m), name = self.n('eta'))),
				self.n('sigma'): gpy(pd.Series(0.5, index = self.get('kninp',m=m), name = self.n('sigma'))),
				self.n('qnorm_out'): gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('knout',m=m)],self.s.db),name=self.n('qnorm_out')),**{'type':'parameter'}),
				self.n('qnorm_inp'): gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('kninp',m=m)],self.s.db),name=self.n('qnorm_inp')),**{'type':'parameter'}),
				self.n('qiv_out'): gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('spout',m=m)],self.s.db), name = self.n('qiv_out'))),
				self.n('qiv_inp'): gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('spinp',m=m)],self.s.db), name = self.n('qiv_inp'))),
				self.n('Rrate'): gpy(pd.Series(self.get('R_LR'), index = self.get('t'), name = self.n('Rrate'))),
				self.n('rDepr'): gpy(pd.Series(0.075, index = adjMultiIndexDB.mergeDomains([self.get('t'),self.get('dur',m=m)],self.s.db), name = self.n('rDepr',m=m))),
				self.n('icpar1'): gpy(pd.Series(0.1, index = self.get('dur',m=m), name= self.n('icpar1',m=m))),
				self.n('icpar2'): gpy(pd.Series(0.095, index = self.get('dur',m=m), name= self.n('icpar2',m=m))),
				self.n('K_tvc'): gpy(pd.Series(0, index = self.get('dur',m=m), name=self.n('K_tvc',m=m))),
				self.n('ic'): gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('output',m=m)], self.s.db), name= self.n('ic',m=m))),
				self.n('p'): gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('output_n',m=m)],self.s.db), name=self.n('p'))),
				self.n('markup'): gpy(pd.Series(1, index = self.get('s',m=m)), name = self.n('markup')),
				self.n('tauS'): gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('output',m=m)],self.s.db), name=self.n('tauS'))),
				self.n('tauD'): gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('input',m=m)], self.s.db), name=self.n('tauD'))),
				self.n('outShare'): gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('output',m=m)], self.s.db), name = self.n('outShare'))),
				self.n('TotalTax'): gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('s')], self.s.db), name = self.n('TotalTax'))),
				self.n('tauLump'): gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('s')], self.s.db), name = self.n('tauLump')))}

	def initDurables(self):
		gpyDB_wheels.robust.robust_merge_dbs(self.s.db,
		{self.n('rDepr'): gpy(pd.Series(0.075, index = adjMultiIndexDB.mergeDomains([self.get('t'),self.get('dur')],self.s.db), name = self.n('rDepr'))),
		 self.n('icpar1'): gpy(pd.Series(0.1, index = self.get('dur'), name= self.n('icpar1'))),
		 self.n('K_tvc') : gpy(pd.Series(0, index = self.get('dur'), name=self.n('K_tvc')))},
		priority='first')
		self.adjustToSteadyState()

	def adjustToSteadyState(self,m=None):
		gpyDB_wheels.robust.robust_merge_dbs(self.s.db,
		{self.n('qD'): gpy(adjMultiIndex.applyMult((self.get('g_LR')+self.get('rDepr',m=m))*gpyDB_wheels.adj.rc_pd(self.get('qD'), self.get('dur',m=m)), self.get('dur2inv',m=m)).droplevel('n').rename_axis(index={'nn':'n'}), **{'name': self.n('qD')}),
		 self.n('icpar2'): gpy(self.get('rDepr',m=m).xs(self.get('t0')[0],level=self.n('t')), **{'name': self.n('icpar2')}),
		 self.n('pD'): gpy(adjMultiIndex.applyMult(gpyDB_wheels.adj.rc_pd(self.get('pD'), self.get('dur',m=m))/(self.get('Rrate')/(1+self.get('infl_LR'))-(1-self.get('rDepr',m=m))), self.get('dur2inv',m=m)).droplevel('n').rename_axis(index={'nn':'n'}), **{'name':self.n('pD')})}
		,priority='second')

	def groups(self,m=None):
		return {g.name: g for g in self.groups_(m=m)}
	def states(self,m=None):
		return {k: self.s.standardInstance(state=k) | {attr: getattr(self,attr)()[k] for attr in ('g_endo','g_exo','blocks','args')} for k in ('B','C')}
	def args(self):
		return {k: {self.name+'_Blocks': '\n'.join([getattr(GamsProduction, module.f)(name, self.name) for name,module in self.m.items()])} for k in ('B','C')}
	def blocks(self):
		return {k: OrdSet([f"B_{name}" for name in self.m]) for k in ('B','C')}
	def g_endo(self):
		return {'B': OrdSet([f"G_{self.name}_endo_always", f"G_{self.name}_exo_in_calib", f"G_{self.name}_endo_dur"]),
				'C': OrdSet([f"G_{self.name}_endo_always", f"G_{self.name}_endo_in_calib",f"G_{self.name}_endo_dur"])}
	def g_exo(self):
		return {'B': OrdSet([f"G_{self.name}_exo_always", f"G_{self.name}_endo_in_calib",f"G_{self.name}_exo_dur"]),
				'C': OrdSet([f"G_{self.name}_exo_always", f"G_{self.name}_exo_in_calib", f"G_{self.name}_exo_dur"])}
	def groups_(self,m=None):
		return [GmsPy.Group(f"G_{self.name}_exo_always", 
		v = 	[('qS', self.g('output')), 
				 ('pD', self.g('input')), 
				 ('sigma', self.g('kninp')), 
				 ('eta', self.g('knout')), 
				 ('mu', self.g('exomu')),
				 ('tauS', self.g('output')),
				 ('tauD', self.g('input')),
				 ('tauLump', ('and', [self.g('s'), self.g('tx0E')])),
				 ('p', ('not', self.g('output_n')))],
		neg_v = [('qS', ('and', [self.g('endo_qS'),self.g('t0')]))]
				),
				GmsPy.Group(f"G_{self.name}_endo_always",
		v = 	[('pD', self.g('int')),
				 ('pS', self.g('output')),
				 ('p' , ('and', [self.g('output_n'), self.g('tx0')])),
				 ('qD', ('and', [('or', [self.g('int'), self.g('input')]), self.g('tx0')])),
				 ('qD', ('and', [self.g('endo_qD'), self.g('t0')])),
				 ('qiv_inp', self.g('spinp')),
				 ('qiv_out', self.g('spout')),
				 ('outShare', self.g('output')),
				 ('TotalTax', ('and', [self.g('s'), self.g('tx0E')]))]),
				GmsPy.Group(f"G_{self.name}_exo_in_calib",
		v = 	[('qD', ('and', [('or', [self.g('int'), self.g('input')]), self.g('t0')])),
				 ('p' , ('and', [self.g('output_n'), self.g('t0')])),
				 ('TotalTax', ('and', [self.g('s'), self.g('t0')]))],
		neg_v = [('qD', ('and', [self.g('endo_qD'), self.g('t0')]))]),
				GmsPy.Group(f"G_{self.name}_endo_in_calib",
		v = 	[('mu', ('and', [self.g('map'), ('not', self.g('exomu'))])),
				 ('qS', ('and', [self.g('endo_qS'), self.g('t0')])),
				 ('tauLump', ('and', [self.g('s'), self.g('t0')])),
				 ('markup', self.g('s'))]),
				GmsPy.Group(f"G_{self.name}_exo_dur", 
		v = [('Rrate', None),
			 ('rDepr', self.g('dur')),
			 ('icpar1', self.g('dur')),
			 ('icpar2', self.g('dur')),
			 ('K_tvc' , self.g('dur')),
			 ('qD', ('and', [self.g('dur'), self.g('t0')]))
			 ]),
				GmsPy.Group(f"G_{self.name}_endo_dur",
		v = [('qD', ('and', [self.g('dur'), self.g('tx0')])),
			 ('pD', ('and', [self.g('dur'), self.g('txE')])),
			 ('ic', ('and', [self.g('output'), self.g('txE')]))
		]
		)]


class Production_ExoMu(Production):
	def __init__(self, f=None, tree=None, ns=None, s=None, glob=None, s_kwargs = None, g_kwargs = None, dur_kwargs = None):
		super().__init__(f=f, tree=tree, ns = ns, s = s, glob=glob, s_kwargs = s_kwargs, g_kwargs = g_kwargs, dur_kwargs=dur_kwargs)

	def addDurables(self, dur = None, f = None, dur2inv = None):
		self.ns.update({k: f"{k}_{self.name}" for k in ('dur','inv')})
		self.s.db[self.n('dur')] = noneInit(dur, pd.MultiIndex.from_product([self.get('output').levels[0], []], names = ['s','n']))
		self.s.db[self.n('dur2inv')] = pd.MultiIndex.from_frame(self.get('dur').to_frame(index=False).assign(nn=lambda x: 'I_'+x.n)) if dur2inv is None else dur2inv
		self.s.db[self.n('inv')] = self.get('dur2inv').droplevel(self.ns['n']).unique().rename({self.ns['nn']:self.ns['n']})
		self.s.db['n'] = self.get('n').union(self.get('inv').levels[-1])
		self.s.db[self.n('input')] = self.get('input').difference(self.get('dur')).union(self.get('inv'))
		self.m[f"{self.name}_IC"] = Submodule(**{'f': noneInit(f, 'sqrAdjCosts')})

	def addPriceWedge(self,f=None):
		self.ns.update({k: f"{k}_{self.name}" for k in ['s','output_n']})
		self.s.db[self.n('output_n')] = self.get('output').levels[-1]
		self.s.db[self.n('s')] = self.get('output').levels[0]
		self.m[f"{self.name}_pWedge"] = Submodule(**{'f': noneInit(f,'priceWedge')})

	def readTree(self,tree):
		gpyDB_wheels.robust.robust_merge_dbs(self.s.db,tree.db,priority='second')
		self.ns.update(tree.ns)
		[self.readTree_i(t) for t in tree.trees.values()];
		self.addCalibrationSubsets(tree)

	def readTree_i(self,t):
		self.addModule(t,**{k:v for k,v in t.__dict__.items() if k in ('name','ns','f','io','sp')})

	def addCalibrationSubsets(self,tree):
		""" endo_mu is a subset of share parameters to endogenize when calibrating to IO data. 
			endo_qD is a subset of demand variables to keep endogeonus when calibrating, even though they might be a part the 'inputs' subset """
		self.ns.update({k: k+'_'+self.name for k in ('endo_mu','endo_qD')})
		endo_mu_all = gpyDB_wheels.adj.rc_pd(self.get('map'), self.get('input').rename(['s','nn'])) # identify all elements that can be endogenized
		spinp_Trouble =gpyDB_wheels.adj.rc_pd(self.get('map_spinp'), ('not', self.get('map_spinp').difference(endo_mu_all).droplevel('nn').unique())) # all nests that might be an issue
		mu_out = gpyDB_wheels.adj.rc_pd(self.get('map'), self.get('output')) # All potentially endogenizeable share parameters on the output side - the ones with multiple outputs:
		s_numberOut = mu_out.to_frame(index=False).groupby('s').nunique()['n']
		mout_s =s_numberOut[s_numberOut>1].index
		mu_out = gpyDB_wheels.adj.rc_pd(mu_out, mout_s)
		endo_mu = self.uniqueFromMap(endo_mu_all.difference(mu_out), gb = ['s','nn']) # suggest endogenous combination that does not include elements  from mu_out
		x = self.uniqueFromMap(endo_mu_all, gb = ['s','nn']) # suggest endogenous combination from all endo_mu_all
		endo_mu = endo_mu.union(gpyDB_wheels.adj.rc_pd(x, ('not', endo_mu.droplevel('n')))) # if there is a combination of (s,nn) in x that is missing, add it here. 
		exo_mu  = self.uniqueFromMap(spinp_Trouble, gb = ['s','n']) # suggest making this element exogenous
		exo_mu_out = self.uniqueFromMap(mu_out.intersection(endo_mu), gb ='s') # suggest NOT making this endogenous. 
		y = self.uniqueFromMap(mu_out, gb = 's') # suggest not making this endogenous either - selecting one for each
		exo_mu_out = exo_mu_out.union(gpyDB_wheels.adj.rc_pd(y, ('not', exo_mu_out.droplevel('nn')))) # if there is a combination of (s,n) in y that is missing, add it here.
		endo_mu_out = self.uniqueFromMap(gpyDB_wheels.adj.rc_pd(mu_out, ('not', exo_mu_out.droplevel('nn'))), gb = ['s','n'])
		if not len(endo_mu_out) == sum(mu_out.to_frame(index=False).groupby('s').nunique()['n']-1):
			print(f"Check that grouping of endogenous/exogenous mu parameters in {self.name}; unexpected number of parameters were endogenized")
		self.s.db[self.ns['endo_mu']] = endo_mu.union(endo_mu_out).difference(exo_mu)
		self.s.db[self.ns['endo_qD']] = exo_mu.droplevel('n').rename(['s','n'])
		
	def uniqueFromMap(self,map_,gb=('s','n')):
		""" MultiIndex-like groupby statement with function 'first' """
		return pd.MultiIndex.from_frame(map_.to_frame(index=False).groupby([s for s in gb]).first().reset_index()).reorder_levels(map_.names)

	def groups_(self,m=None):
		return [GmsPy.Group(f"G_{self.name}_exo_always", 
		v = 	[('qS', self.g('output')), 
				 ('pD', self.g('input')), 
				 ('sigma', self.g('kninp')), 
				 ('eta', self.g('knout')),
				 ('mu', self.g('map')),
				 ('tauS', self.g('output')),
				 ('tauD', self.g('input')),
				 ('tauLump', ('and', [self.g('s'), self.g('tx0E')])),
				 ('p', ('not', self.g('output_n')))],
		neg_v = [('mu', self.g('endo_mu'))]
				),
				GmsPy.Group(f"G_{self.name}_endo_always",
		v = 	[('pD', self.g('int')),
				 ('pS', self.g('output')),
				 ('p' , ('and', [self.g('output_n'), self.g('tx0')])),
				 ('qD', self.g('int')),
				 ('qD', ('and', [self.g('input'), self.g('tx0')])),
				 ('qD', ('and', [self.g('endo_qD'), self.g('t0')])),
				 ('qiv_inp', self.g('spinp')),
				 ('qiv_out', self.g('spout')),
				 ('outShare', self.g('output')),
				 ('TotalTax', ('and', [self.g('s'), self.g('tx0E')]))]),
				GmsPy.Group(f"G_{self.name}_exo_in_calib",
		v = 	[('qD', ('and', [self.g('input'), self.g('t0')])),
				 ('p' , ('and', [self.g('output_n'), self.g('t0')])),
				 ('TotalTax', ('and', [self.g('s'), self.g('t0')]))],
		neg_v = [('qD', ('and', [self.g('endo_qD'), self.g('t0')]))]),
				GmsPy.Group(f"G_{self.name}_endo_in_calib",
		v = 	[('mu', self.g('endo_mu')),
				 ('tauLump', ('and', [self.g('s'), self.g('t0')])),
				 ('markup', self.g('s'))]),
				GmsPy.Group(f"G_{self.name}_exo_dur", 
		v = [('Rrate', None),
			 ('rDepr', self.g('dur')),
			 ('icpar1', self.g('dur')),
			 ('icpar2', self.g('dur')),
			 ('K_tvc' , self.g('dur')),
			 ('qD', ('and', [self.g('dur'), self.g('t0')]))
			 ]),
				GmsPy.Group(f"G_{self.name}_endo_dur",
		v = [('qD', ('and', [self.g('dur'), self.g('tx0')])),
			 ('pD', ('and', [self.g('dur'), self.g('txE')])),
			 ('ic', ('and', [self.g('output'), self.g('txE')]))
		]
		)]

class Inventory(GmsPython):
	def __init__(self, f = None, name = None, db_IO = None, itory = None, s=None, glob=None, s_kwargs = None, g_kwargs = None):
		""" Initialize from name, io data, and subset of inventory"""
		super().__init__(name=f"{name}_itory", f=f, s=s, glob=glob, g_kwargs=g_kwargs, s_kwargs = noneInit(s_kwargs, {}) | {'db': db_IO})
		if f is None:
			self.initItory(itory = itory)
	def initItory(self, itory = None):
		self.s.db['s_itory'] = noneInit(itory, pd.Index(['itory'], name = 's'))
		self.s.db['d_itory'] = gpyDB_wheels.adj.rc_pd(self.s.db.get('qD'), self.get('s_itory')).index.droplevel('t').unique()
		self.s.db['inventoryAR'] = gpy(pd.Series(0.5, index = self.get('d_itory'), name = 'inventoryAR'), **{'type': 'parameter'})
	def states(self, m= None):
		return {'B': self.s.standardInstance(state='B') | {attr: getattr(self,attr)() for attr in ('g_endo','g_exo','blocks','args')}}
	def args(self):
		return {f"{self.name}_Blocks": self.equationText}
	def blocks(self):
		return OrdSet([f"B_{self.name}"])
	def g_endo(self):
		return OrdSet([f"G_{self.name}_endo"])
	def g_exo(self):
		return OrdSet([f"G_{self.name}_exo"])
	def groups(self,m=None):
		return {g.name: g for g in [GmsPy.Group(f"G_{self.name}_exo",  v = [('qD', ('and', [self.g('t0'), self.g('d_itory')]))]), 
									GmsPy.Group(f"G_{self.name}_endo", v = [('qD', ('and', [self.g('tx0E'), self.g('d_itory')]))])]}
	@property
	def equationText(self):
		return f"""
$BLOCK B_{self.name}
	E_{self.name}[t,s,n]$(d_itory[s,n] and tx0E[t])..	qD[t,s,n] =E= inventoryAR[s,n] * qD[t-1,s,n];
$ENDBLOCK
"""
