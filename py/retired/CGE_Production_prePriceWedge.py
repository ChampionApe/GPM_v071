from CGE_GmsPython import *
import GamsProduction

class Production(GmsPython):
	def __init__(self, f=None, tree=None, ns=None, s=None, glob=None, s_kwargs = None, g_kwargs = None, dur_kwargs = None):
		""" Initialize from a pickle file 'f' or nesting tree 'tree'. """
		super().__init__(name=tree.name if tree else None, f=f, s=s, glob=glob, ns=ns, s_kwargs = s_kwargs, g_kwargs=g_kwargs)
		if f is None:
			self.readTree(tree)
			self.addDurables(**noneInit(dur_kwargs, {}))

	def addDurables(self, dur = None, f = None, dur2inv = None):
		self.ns.update({k: f"{k}_{self.name}" for k in ('dur','inv')})
		self.s.db[self.n('dur')] = noneInit(dur, pd.MultiIndex.from_product([self.get('output').levels[0], []], names = ['s','n']))
		self.s.db[self.n('dur2inv')] = pd.MultiIndex.from_frame(self.get('dur').to_frame(index=False).assign(nn=lambda x: 'I_'+x.n)) if dur2inv is None else dur2inv
		self.s.db[self.n('inv')] = self.get('dur2inv').droplevel(self.ns['n']).unique().rename({self.ns['nn']:self.ns['n']})
		self.s.db['n'] = self.get('n').union(self.get('inv').levels[-1])
		self.s.db[self.n('input')] = self.get('input').difference(self.get('dur')).union(self.get('inv'))
		self.m['IC'] = Submodule(**{'f': noneInit(f, 'sqrAdjCosts')})

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
		return pd.MultiIndex.from_frame(map_.to_frame(index=False).groupby([self.n(s) for s in gb]).first().reset_index()).reorder_levels(map_.names)

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
		return {self.n('pS'): gpyDB.gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('output',m=m)],self.s.db), name=self.n('pS'))),
				self.n('qS'): gpyDB.gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('output',m=m)],self.s.db), name=self.n('qS'))),
				self.n('pD'): gpyDB.gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('int',m=m).union(self.get('input',m=m)).union(self.get('dur',m=m))],self.s.db), name = self.n('pD'))),
				self.n('qD'): gpyDB.gpy(pd.Series(0.5, index = adjMultiIndexDB.mergeDomains([self.get('t'),self.get('int',m=m).union(self.get('input',m=m)).union(self.get('dur',m=m))],self.s.db), name = self.n('qD'))),
				self.n('mu'): gpyDB.gpy(pd.Series(1, index = self.get('map',m=m), name = self.n('mu'))),
				self.n('eta'): gpyDB.gpy(pd.Series(0.5, index = self.get('knout',m=m), name = self.n('eta'))),
				self.n('sigma'): gpyDB.gpy(pd.Series(0.5, index = self.get('kninp',m=m), name = self.n('sigma'))),
				self.n('qnorm_out'): gpyDB.gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('knout',m=m)],self.s.db),name=self.n('qnorm_out')),**{'type':'parameter'}),
				self.n('qnorm_inp'): gpyDB.gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('kninp',m=m)],self.s.db),name=self.n('qnorm_inp')),**{'type':'parameter'}),
				self.n('qiv_out'): gpyDB.gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('spout',m=m)],self.s.db), name = self.n('qiv_out'))),
				self.n('qiv_inp'): gpyDB.gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('spinp',m=m)],self.s.db), name = self.n('qiv_inp'))),
				self.n('Rrate'): gpyDB.gpy(pd.Series(self.get('R_LR'), index = self.get('t'), name = self.n('Rrate'))),
				self.n('rDepr'): gpyDB.gpy(pd.Series(0.05, index = adjMultiIndexDB.mergeDomains([self.get('t'),self.get('dur',m=m)],self.s.db), name = self.n('rDepr',m=m))),
				self.n('icpar1'): gpyDB.gpy(pd.Series(0.1, index = self.get('dur',m=m), name= self.n('icpar1',m=m))),
				self.n('icpar2'): gpyDB.gpy(pd.Series(0.05, index = self.get('dur',m=m), name= self.n('icpar2',m=m))),
				self.n('K_tvc'): gpyDB.gpy(pd.Series(0, index = self.get('dur',m=m), name=self.n('K_tvc',m=m))),
				self.n('ic'): gpyDB.gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('output',m=m)], self.s.db), name= self.n('ic',m=m)))}

	def initDurables(self):
		gpyDB_wheels.robust.robust_merge_dbs(self.s.db,
		{self.n('rDepr'): gpyDB.gpy(pd.Series(0.05, index = adjMultiIndexDB.mergeDomains([self.get('t'),self.get('dur')],self.s.db), name = self.n('rDepr'))),
		 self.n('icpar1'): gpyDB.gpy(pd.Series(0.1, index = self.get('dur'), name= self.n('icpar1'))),
		 self.n('K_tvc') : gpyDB.gpy(pd.Series(0, index = self.get('dur'), name=self.n('K_tvc')))},
		priority='first')
		self.adjustToSteadyState()

	def adjustToSteadyState(self,m=None):
		gpyDB_wheels.robust.robust_merge_dbs(self.s.db,
		{self.n('qD'): gpyDB.gpy(adjMultiIndex.applyMult((self.get('g_LR')+self.get('rDepr',m=m))*gpyDB_wheels.adj.rc_pd(self.get('qD'), self.get('dur',m=m)), self.get('dur2inv',m=m)).droplevel('n').rename_axis(index={'nn':'n'}), **{'name': self.n('qD')}),
		 self.n('icpar2'): gpyDB.gpy(self.get('rDepr',m=m).xs(self.get('t0')[0],level=self.n('t')), **{'name': self.n('icpar2')}),
		 self.n('pD'): gpyDB.gpy(adjMultiIndex.applyMult(gpyDB_wheels.adj.rc_pd(self.get('pD'), self.get('dur',m=m))/(self.get('Rrate')/(1+self.get('infl_LR'))-(1-self.get('rDepr',m=m))), self.get('dur2inv',m=m)).droplevel('n').rename_axis(index={'nn':'n'}), **{'name':self.n('pD')})}
		,priority='second')

	def groups(self,m=None):
		return {g.name: g for g in self.groups_(m=m)}
	def states(self,m=None):
		return {k: self.s.standardInstance(state=k) | {attr: getattr(self,attr)()[k] for attr in ('g_endo','g_exo','blocks','args')} for k in ('B','C')}
	def args(self):
		return {k: {self.name+'_Blocks': '\n'.join([getattr(GamsProduction, module.f)(f"{self.name}_{name}", name, self.name) for name,module in self.m.items()])} for k in ('B','C')}
	def blocks(self):
		return {k: OrdSet([f"B_{self.name}_{name}" for name in self.m]) for k in ('B','C')}
	def g_endo(self):
		return {'B': OrdSet([f"G_{self.name}_endo_always", f"G_{self.name}_exo_in_calib", f"G_{self.name}_endo_dur"]),
				'C': OrdSet([f"G_{self.name}_endo_always", f"G_{self.name}_endo_in_calib",f"G_{self.name}_endo_dur"])}
	def g_exo(self):
		return {'B': OrdSet([f"G_{self.name}_exo_always", f"G_{self.name}_endo_in_calib",f"G_{self.name}_exo_dur"]),
				'C': OrdSet([f"G_{self.name}_exo_always", f"G_{self.name}_exo_in_calib", f"G_{self.name}_exo_dur"])}
	def groups_(self,m=None):
		return [GmsPy.Group(f"G_{self.name}_exo_always", 
		v = 	[('qS', self.g('output',m=m)), 
				 ('pD', self.g('input',m=m)), 
				 ('sigma', self.g('kninp',m=m)), 
				 ('eta', self.g('knout',m=m)), 
				 ('mu', self.g('exomu',m=m))],
		neg_v = [('qS', ('and', [self.g('endo_qS',m=m),self.g('t0')]))]
				),
				GmsPy.Group(f"G_{self.name}_endo_always",
		v = 	[('pD', self.g('int',m=m)),
				 ('pS', ('and', [self.g('output',m=m), self.g('tx0')])),
				 ('pS', ('and', [self.g('endo_pS',m=m), self.g('t0')])),
				 ('qD', ('and', [('or', [self.g('int',m=m), self.g('input',m=m)]), self.g('tx0')])),
				 ('qD', ('and', [self.g('endo_qD',m=m), self.g('t0')])),
				 ('qiv_inp', self.g('spinp')),
				 ('qiv_out', self.g('spout'))]),
				GmsPy.Group(f"G_{self.name}_exo_in_calib",
		v = 	[('qD', ('and', [('or', [self.g('int',m=m), self.g('input',m=m)]), self.g('t0')])),
				 ('pS', ('and', [self.g('output',m=m), self.g('t0')]))],
		neg_v = [('qD', ('and', [self.g('endo_qD',m=m), self.g('t0')])),
				 ('pS', ('and', [self.g('endo_pS',m=m), self.g('t0')]))]),
				GmsPy.Group(f"G_{self.name}_endo_in_calib",
		v = 	[('mu', ('and', [self.g('map',m=m), ('not', self.g('exomu',m=m))])),
				 ('qS', ('and', [self.g('endo_qS',m=m), self.g('t0')]))]),
				GmsPy.Group(f"G_{self.name}_exo_dur", 
		v = [('Rrate', None),
			 ('rDepr', self.g('dur',m=m)),
			 ('icpar1', self.g('dur',m=m)),
			 ('icpar2', self.g('dur',m=m)),
			 ('K_tvc' , self.g('dur',m=m)),
			 ('qD', ('and', [self.g('dur',m=m), self.g('t0')]))
			 ]),
				GmsPy.Group(f"G_{self.name}_endo_dur",
		v = [('qD', ('and', [self.g('dur',m=m), self.g('tx0')])),
			 ('pD', ('and', [self.g('dur',m=m), self.g('txE')])),
			 ('ic', ('and', [self.g('output',m=m), self.g('txE')]))
		]
		)]