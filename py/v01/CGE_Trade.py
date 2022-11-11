from CGE_GmsPython import *
import GamsTrade

class SimpleArmington(GmsPython):
	def __init__(self, f = None, name = None, ns = None, s = None, glob = None, s_kwargs = None, g_kwargs = None, kwargs=None):
		super().__init__(name=name, f=f, s=s, glob=glob, ns=ns, s_kwargs = s_kwargs, g_kwargs=g_kwargs)
		if f is None:
			self.adjust(**kwargs)

	def adjust(self, sfor_ndom=None, dom2for = None):
		""" The mapping sfor_ndom is a local mapping for this module """
		self.ns.update({k: f"{k}_{self.name}" for k in ['s','sfor_ndom', 'nOut']})
		self.s.db[self.n('sfor_ndom')] = sfor_ndom
		if 'dom2for' not in self.s.db.symbols:
			self.s.db['dom2for'] = dom2for
		self.s.db[self.n('nOut')] = gpyDB_wheels.adj.rc_pd(self.get('dom2for'), self.get('sfor_ndom')).droplevel('n').rename('n').unique()
		self.s.db[self.n('s')] = self.get('sfor_ndom').levels[0]

	def initDB(self,m=None):
		return gpyDB_wheels.robust.robust_merge_dbs(self.s.db,self.initSymbols(),priority='first')

	def initSymbols(self):
		return {'p': gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('nOut')], self.s.db), name = self.n('p'))), 
				'qD': gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('sfor_ndom')], self.s.db), name = self.n('qD'))),
				'pD': gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('sfor_ndom')], self.s.db), name = self.n('pD'))),
				'Fscale': gpy(pd.Series(1, index = self.get('sfor_ndom'), name = self.n('Fscale'))),
				'tauD': gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('sfor_ndom')], self.s.db), name = self.n('tauD'))),
				'tauLump': gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('s')], self.s.db), name = self.n('tauLump'))),
				'TotalTax': gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('s')], self.s.db), name = self.n('TotalTax'))),
				'sigma': gpy(pd.Series(5, index = self.get('sfor_ndom'), name = self.n('sigma')))}

	def states(self,m=None):
		return {k: self.s.standardInstance(state=k) | {attr: getattr(self,attr)()[k] for attr in ('g_endo','g_exo','blocks','args')} for k in ('B','C')}
	def args(self):
		return {k: {f"{self.name}_Blocks": GamsTrade.Armington(self.name)} for k in ('B','C')}
	def blocks(self):
		return {k: OrdSet([f"B_{self.name}"]) for k in ('B','C')}
	def groups(self,m=None):
		return {g.name: g for g in self.groups_(m=m)}
	def g_endo(self):
		return {'B': OrdSet([f"G_{self.name}_endo_always", f"G_{self.name}_exo_in_calib"]),
				'C': OrdSet([f"G_{self.name}_endo_always", f"G_{self.name}_endo_in_calib"])}
	def g_exo(self):
		return {'B': OrdSet([f"G_{self.name}_exo_always", f"G_{self.name}_endo_in_calib"]),
				'C': OrdSet([f"G_{self.name}_exo_always", f"G_{self.name}_exo_in_calib"])}
	def groups_(self,m=None):
		return [GmsPy.Group(f"G_{self.name}_exo_always", 
			v =[('p', self.g('nOut')), ('sigma', self.g('sfor_ndom')), ('pD', self.g('sfor_ndom')), ('tauD', self.g('sfor_ndom')), ('tauLump', ('and', [self.g('s'), self.g('tx0E')]))]),
				GmsPy.Group(f"G_{self.name}_endo_always", 
			v =[('qD', ('and', [self.g('sfor_ndom'), self.g('tx0E')])), ('TotalTax', ('and', [self.g('s'), self.g('tx0E')]))]),
				GmsPy.Group(f"G_{self.name}_exo_in_calib", 
			v =[('qD', ('and', [self.g('sfor_ndom'), self.g('t0')])), ('TotalTax', ('and', [self.g('s'), self.g('t0')]))]),
				GmsPy.Group(f"G_{self.name}_endo_in_calib",
			v =[('Fscale', self.g('sfor_ndom')), ('tauLump', ('and', [self.g('s'), self.g('t0')]))])
			]