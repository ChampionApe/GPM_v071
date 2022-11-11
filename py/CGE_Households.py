from CGE_GmsPython import *
import GamsHouseholds

class SimpleRamsey(GmsPython):
	def __init__(self, f=None, tree=None, ns=None, s=None, glob=None, s_kwargs = None, g_kwargs = None, kwargs = None):
		""" Initialize from a pickle file 'f' or nesting tree 'tree'. """
		super().__init__(name=tree.name if tree else None, f=f, s=s, glob=glob, ns=ns, s_kwargs = s_kwargs, g_kwargs=g_kwargs)
		if f is None:
			self.readTree(tree)
			self.adjust(**kwargs)

	def adjust(self, L2C = None, svngs = None, f = None):
		""" add mapping from labor component L to top of consumption nest C. """

		self.ns.update({k: f"{k}_{self.name}" for k in ('labor','L2C','svngs','s','output_n','input_n')})
		self.s.db[self.n('L2C')] = noneInit(L2C, pd.MultiIndex.from_tuples([], names = ['s','n','nn']))
		self.s.db[self.n('labor')] = L2C.droplevel('nn').unique()
		self.s.db[self.n('output_n')] = self.get('labor').levels[-1]
		self.s.db[self.n('input_n')]  = self.get('input').levels[-1]
		self.s.db[self.n('s')] = self.get('output').levels[0]
		self.s.db[self.n('svngs')] = noneInit(svngs, pd.MultiIndex.from_tuples([], names = ['s','n']))
		self.s.db['n'] = self.get('n').union(self.get('output_n')).union(self.get('svngs').levels[-1])
		self.m[f"{self.name}_labor"] = Submodule(**{'f': noneInit(f, 'IsoFrisch')})
		self.m[f"{self.name}_dynamic"] = Submodule(**{'f': 'simpleDynamic'})
		self.m[f"{self.name}_pw"] = Submodule(**{'f': 'priceWedge'})

	def readTree(self,tree):
		robust.robust_merge_dbs(self.s.db,tree.db,priority='second')
		self.ns.update(tree.ns)
		[self.readTree_i(t) for t in tree.trees.values()];
		self.addCalibrationSubsets(tree)

	def readTree_i(self,t):
		self.addModule(t,**{k:v for k,v in t.__dict__.items() if k in ('name','ns','f','io','sp')})

	def addCalibrationSubsets(self,tree):
		self.ns.update({k: k+'_'+self.name for k in ['endo_mu']})
		self.s.db[self.n('endo_mu')] = adj.rc_pd(self.get('map'), self.get('input').rename({'n':'nn'}))

	def initDB(self,m=None):
		return robust.robust_merge_dbs(self.s.db,self.initSymbols(),priority='first')

	def initSymbols(self):
		return {self.n('pS'): gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('labor')],self.s.db), name=self.n('pS'))),
				self.n('qS'): gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('labor')],self.s.db), name=self.n('qS'))),
				self.n('pD'): gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('input').union(self.get('int')).union(self.get('output'))],self.s.db), name = self.n('pD'))),
				self.n('qD'): gpy(pd.Series(0.5, index = adjMultiIndexDB.mergeDomains([self.get('t'),self.get('input').union(self.get('int')).union(self.get('output'))],self.s.db), name = self.n('qD'))),
				self.n('mu'): gpy(pd.Series(1, index = self.get('map'), name = self.n('mu'))),
				self.n('eta'): gpy(pd.Series(0.5, index = self.get('knout'), name = self.n('eta'))),
				self.n('sigma'): gpy(pd.Series(0.5, index = self.get('kninp'), name = self.n('sigma'))),
				self.n('qnorm_out'): gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('knout')],self.s.db),name=self.n('qnorm_out')),**{'type':'parameter'}),
				self.n('qnorm_inp'): gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('kninp')],self.s.db),name=self.n('qnorm_inp')),**{'type':'parameter'}),
				self.n('qiv_out'): gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('spout')],self.s.db), name = self.n('qiv_out'))),
				self.n('qiv_inp'): gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('spinp')],self.s.db), name = self.n('qiv_inp'))),
				self.n('Rrate'): gpy(pd.Series(self.get('R_LR'), index = self.get('t'), name = self.n('Rrate'))),
				self.n('iRate'): gpy(pd.Series(self.get('R_LR')*(1+self.get('infl_LR')), index = self.get('t'), name = self.n('iRate'))),
				self.n('crra'): gpy(pd.Series(2, index = self.get('output'), name = self.n('crra'))),
				self.n('frisch'): gpy(pd.Series(0.25, index = self.get('labor'), name = self.n('frisch'))),
				self.n('Lscale'): gpy(pd.Series(1, index = self.get('labor'), name = self.n('Lscale'))),
				self.n('disc'): gpy(pd.Series(1/self.get('R_LR'), index = self.get('s'), name = self.n('disc'))),
				self.n('vD'): gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('t'),self.get('svngs')],self.s.db), name = self.n('vD'))),
				self.n('h_tvc'): gpy(pd.Series(0, index = self.get('svngs'), name = self.n('h_tvc'))),
				self.n('tauLump'): gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('txE'),self.get('s')], self.s.db), name = self.n('tauLump'))),
				self.n('tauS'): gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('labor')], self.s.db), name = self.n('tauS'))),
				self.n('tauD'): gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('input')], self.s.db), name = self.n('tauD'))),
				self.n('TotalTax'): gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('s')], self.s.db), name = self.n('TotalTax'))),
				self.n('sp'): gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('s')],self.s.db), name = self.n('sp'))),
				self.n('p'): gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('output_n').union(self.get('input_n'))], self.s.db), name = self.n('p')))}

	def groups(self,m=None):
		return {g.name: g for g in self.groups_(m=m)}
	def states(self,m=None):
		return {k: self.s.standardInstance(state=k) | {attr: getattr(self,attr)()[k] for attr in ('g_endo','g_exo','blocks','args')} for k in ('B','C')}
	def args(self):
		return {k: {self.name+'_Blocks': '\n'.join([getattr(GamsHouseholds, module.f)(name, self.name) for name,module in self.m.items()])} for k in ('B','C')}
	def blocks(self):
		return {k: OrdSet([f"B_{name}" for name in self.m]) for k in ('B','C')}
	def g_endo(self):
		return {'B': OrdSet([f"G_{self.name}_endo_always", f"G_{self.name}_exo_in_calib"]),
				'C': OrdSet([f"G_{self.name}_endo_always", f"G_{self.name}_endo_in_calib"])}
	def g_exo(self):
		return {'B': OrdSet([f"G_{self.name}_exo_always", f"G_{self.name}_endo_in_calib"]),
				'C': OrdSet([f"G_{self.name}_exo_always", f"G_{self.name}_exo_in_calib"])}
	def groups_(self,m=None):
		return [GmsPy.Group(f"G_{self.name}_exo_always", 
			v =[('sigma', self.g('kninp')),
				('eta', self.g('knout')),
				# ('disc',self.g('s')),
				('crra', self.g('output')),
				('h_tvc', self.g('svngs')),
				('mu', ('and', [self.g('map'), ('not', self.g('endo_mu'))])),
				('Rrate', None), ('iRate', None), ('frisch', self.g('labor')), 
				('tauD', self.g('input')), ('tauS', self.g('labor')), ('tauLump', ('and', [self.g('s'), self.g('tx0E')])), 
				('vD', ('and', [self.g('svngs'), self.g('t0')])),
				('p', ('or', [self.g('output_n'), self.g('input_n')])),
				]),
				GmsPy.Group(f"G_{self.name}_endo_always",
			v =[('pD', ('or', [self.g('int'), self.g('input')])),
				('pD', ('and', [self.g('output'), self.g('tx0E')])),
				('qS', ('and',[self.g('labor'), self.g('tx0E')])),
				('qD', ('and',[self.g('input'), self.g('tx0E')])),
				('qD', ('or', [self.g('int'), self.g('output')])),
				('qiv_inp', self.g('spinp')),
				('qiv_out', self.g('spout')),
				('sp', self.g('s')),
				('pS', self.g('labor')),
				('TotalTax', ('and', [self.g('s'), self.g('tx0E')])),
				('vD', ('and', [self.g('svngs'), self.g('tx0')]))]
			),
				GmsPy.Group(f"G_{self.name}_exo_in_calib", 
			v =[('qD', ('and', [self.g('input'), self.g('t0')])),
				('qS', ('and', [self.g('labor'), self.g('t0')])),
				('pD', ('and', [self.g('output'),self.g('t0')])),
				('TotalTax', ('and', [self.g('s'), self.g('t0')]))]),
				GmsPy.Group(f"G_{self.name}_endo_in_calib",
			v =[('mu', self.g('endo_mu')),
				('Lscale', self.g('labor')),
				('tauLump', ('and', [self.g('s'), self.g('t0')])),
				# ('crra', self.g('output'))
				# ('h_tvc', self.g('svngs'))
				('disc', self.g('s'))
				])
			]