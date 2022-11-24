from CGE_GmsPython import *
import GamsGovernment

class balancedBudget(GmsPython):
	def __init__(self, f=None, tree=None, ns=None, s=None, glob=None, s_kwargs=None, g_kwargs=None, kwargs=None):
		""" Initialize from a pickle file 'f' or nesting tree 'tree'. """
		super().__init__(name=tree.name if tree else None, f=f, s=s, glob=glob, ns=ns, s_kwargs=s_kwargs, g_kwargs=g_kwargs)
		if f is None:
			self.readTree(tree)
			self.adjust(**noneInit(kwargs,{}))

	def adjust(self, standAlone=True, balanceInstrument='laborTax'):
		self.ns.update({k: f"{k}_{self.name}" for k in ('s', 'input_n')})
		self.s.db[self.n('input_n')] = self.get('input').levels[-1]
		self.s.db[self.n('s')] = self.get('input').levels[0]
		self.m[f"{self.name}_dynamic"] = Submodule(**{'f': 'simpleDynamics'})
		self.m[f"{self.name}_bb"] = Submodule(**{'f': 'balancedBudget'})
		self.m[f"{self.name}_calib"] = Submodule(**{'f': 'calibrationFlat'})
		self.balanceInstrument = balanceInstrument
	
	@property	
	def getBalanceInstrument(self):
		if self.balanceInstrument == 'laborTax':
			return ('tauS', ('and', [self.g('labor'), self.g('tx0E')]))

	def readTree(self, tree):
		robust.robust_merge_dbs(self.s.db, tree.db, priority='second')
		self.ns.update(tree.ns)
		[self.readTree_i(t) for t in tree.trees.values()]
		self.addCalibrationSubsets(tree)

	def readTree_i(self, t):
		self.addModule(t, **{k: v for k, v in t.__dict__.items() if k in ('name', 'ns', 'f', 'io', 'sp')})

	def addCalibrationSubsets(self, tree):
		self.ns.update({k: k + '_' + self.name for k in ['endo_mu']})
		self.s.db[self.n('endo_mu')] = adj.rc_pd(self.get('map'), self.get('input').rename({'n': 'nn'}))

	def initDB(self, m=None):
		return robust.robust_merge_dbs(self.s.db, self.initSymbols(), priority='first')

	def initSymbols(self):
		return {self.n('pD'): gpy(pd.Series(1, index=adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('input').union(self.get('int')).union(self.get('output'))], self.s.db), name=self.n('pD'))),
				self.n('qD'): gpy(pd.Series(0.5, index=adjMultiIndexDB.mergeDomains([self.get('t'), self.get('input').union(self.get('int')).union(self.get('output'))], self.s.db), name=self.n('qD'))),
				self.n('mu'): gpy(pd.Series(1, index=self.get('map'), name=self.n('mu'))),
				self.n('sigma'): gpy(pd.Series(0.5, index=self.get('kninp'), name=self.n('sigma'))),
				self.n('qnorm_inp'): gpy(pd.Series(0, index=adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('kninp')], self.s.db), name=self.n('qnorm_inp')), **{'type': 'parameter'}),
				self.n('qiv_inp'): gpy(pd.Series(1, index=adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('spinp')], self.s.db), name=self.n('qiv_inp'))),
				self.n('Rrate'): gpy(pd.Series(self.get('R_LR'), index=self.get('t'), name=self.n('Rrate'))),
				self.n('iRate'): gpy(pd.Series(self.get('R_LR') * (1 + self.get('infl_LR')), index=self.get('t'), name=self.n('iRate'))),
				self.n('vAssets'): gpy(pd.Series(1, index=pd.MultiIndex.from_product([self.get('t'), self.get('s'), pd.Index(['total'], name='a')]), name=self.n('vAssets'))),
				self.n('tauD'): gpy(pd.Series(0.01, index=adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('input')], self.s.db), name=self.n('tauD'))),
				self.n('TotalTax'): gpy(pd.Series(0.01, index=adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('d_TotalTax')], self.s.db), name=self.n('TotalTax'))),
				self.n('sp'): gpy(pd.Series(0, index=adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('s')], self.s.db), name=self.n('sp'))),
				self.n('p'): gpy(pd.Series(1, index=adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('input_n')], self.s.db), name=self.n('p'))),
				self.n('tauD0'): gpy(pd.Series(0.01, index=adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('input')], self.s.db), name=self.n('tauD0'))),
				self.n('tauG_calib'): gpy(1, name=self.n('tauG_calib')),
				self.n('tauS'): gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('labor')], self.s.db), name = self.n('tauS'))),
				self.n('jG_budget'): gpy(0, name = 'jG_budget', **{'type': 'parameter'})}

	def groups(self, m=None):
		return {g.name: g for g in self.groups_(m=m)}

	def states(self, m=None):
		return {k: self.s.standardInstance(state=k) | {attr: getattr(self, attr)()[k] for attr in ('g_endo', 'g_exo', 'blocks', 'args')} for k in ('B', 'C', 'B_standAlone', 'C_standAlone')}
	def args(self):
		return {k: {self.name + '_Blocks': '\n'.join([getattr(GamsGovernment, module.f)(name, self.name) for name, module in self.m.items()])} for k in ('B', 'C', 'B_standAlone', 'C_standAlone')}
	def blocks(self):
		full = {'C': OrdSet([f"B_{name}" for name in self.m])}
		full['B'] = full['C']-OrdSet([f"B_{self.name}_calib"])
		standAlone = {k+'_standAlone': full[k]-OrdSet([f"B_{self.name}_bb"]) for k in ('B','C')}
		return full | standAlone

	def g_endo(self):
		full = {'B': OrdSet([f"G_{self.name}_endo_always", f"G_{self.name}_exo_in_calib", f"G_{self.name}_balanceBudget"]),
				'C': OrdSet([f"G_{self.name}_endo_always", f"G_{self.name}_endo_in_calib", f"G_{self.name}_balanceBudget"])}
		standAlone = {k+'_standAlone': full[k]-OrdSet([f"G_{self.name}_balanceBudget"]) for k in ('B','C')}
		return full | standAlone

	def g_exo(self):
		full = {'B': OrdSet([f"G_{self.name}_exo_always", f"G_{self.name}_endo_in_calib"]),
				'C': OrdSet([f"G_{self.name}_exo_always", f"G_{self.name}_exo_in_calib"])}
		standAlone = {k+'_standAlone': full[k] for k in ('B','C')}
		return full | standAlone

	def groups_(self, m=None):
		return [GmsPy.Group(f"G_{self.name}_exo_always",
							v=[('sigma', self.g('kninp')),
								('mu', ('and', [self.g('map'),('not', self.g('endo_mu'))])),
								('Rrate', None), ('iRate', None),
								('vAssets', ('and', [self.g('s'), self.g('t0')])),
								('p', self.g('input_n')),
								('TotalTax', ('and', [self.g('d_TotalTax'), ('not', self.g('s'))])),
								('qD', self.g('output')),
								('tauD0', self.g('input'))
							   ]),
				GmsPy.Group(f"G_{self.name}_endo_always",
							v=[('pD', ('or', [self.g('int'), self.g('input'), self.g('output')])),
								('qD', ('and', [self.g('input'), self.g('tx0E')])),
								('qD', self.g('int')),
								('qiv_inp', self.g('spinp')),
								('sp', self.g('s')),
								('TotalTax',('and', [self.g('s'), self.g('tx0E')])),
								('vAssets', ('and', [self.g('s'), self.g('tx0')]))]
							),
				GmsPy.Group(f"G_{self.name}_exo_in_calib",
							v=[('qD', ('and', [self.g('input'), self.g('t0')])),
							   ('TotalTax', ('and', [self.g('s'), self.g('t0')]))]),
				GmsPy.Group(f"G_{self.name}_endo_in_calib",
							v=[('mu', self.g('endo_mu')),
								('tauG_calib', None),
								('tauD', self.g('input'))
							   ]),
				GmsPy.Group(f"G_{self.name}_balanceBudget", 
							v=[self.getBalanceInstrument])
				]