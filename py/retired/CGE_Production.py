from CGE_GmsPython import *
import GamsProduction

class Production(GmsPython):
	def __init__(self, f=None, tree=None, ns=None, s=None, glob=None, s_kwargs = None, g_kwargs = None):
		""" Initialize from a pickle file 'f' or nesting tree 'tree'. """
		super().__init__(name=tree.name if tree else None, f=f, s=s, glob=glob, ns=ns, s_kwargs = s_kwargs, g_kwargs=g_kwargs)
		if f is None:
			self.readTree(tree)

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
		self.s.db[self.ns['endo_qD']] = exomu_inp.droplevel(self.n('n')).rename([self.n('s'),self.n('n')]).union(gpyDB_wheels.adj.rc_pd(exomu_out.droplevel(self.n('nn')),('not',self.get('output'))))
		self.s.db[self.ns['endo_qS']] = gpyDB_wheels.adj.rc_pd(exomu_out.droplevel(self.n('nn')),self.get('output'))
		self.s.db[self.ns['endo_pS']] = self.uniqueFromMap(self.get('output'),gb=['s'])

	def TroubleNodes(self,map_spinp, map_spout, output=None):
		""" Identify nodes that are branches in input tree and output tree"""
		return gpyDB_wheels.adj.rc_pd(map_spinp.droplevel(self.n('n')).rename([self.n('s'),self.n('n')]), c = gpyDB_wheels.adj.rc_pd(map_spout.droplevel(self.n('nn')), c = None if output is None else ('not', output)))

	def uniqueFromMap(self,map_,gb=('s','n')):
		""" MultiIndex-like groupby statement with function 'first' """
		return pd.MultiIndex.from_frame(map_.to_frame(index=False).groupby([self.n(s) for s in gb]).first().reset_index()).reorder_levels(map_.names)

	def getCleanExoMu(self,map_spinp, map_spout, trouble):
		""" return elements to be added to exomu_inp, exomu_out"""
		return self.uniqueFromMap(gpyDB_wheels.adj.rc_pd(map_spinp, c= ('not', trouble.rename([self.n('s'),self.n('nn')])))), self.uniqueFromMap(gpyDB_wheels.adj.rc_pd(map_spout, c= ('not', trouble)),gb=('s','nn'))

	def getExomuFromTree(self,tree,maxiter=10):
		map_spinp = self.get('map_spinp').copy()
		map_spout = self.get('map_spout').copy()
		trouble = self.TroubleNodes(map_spinp, map_spout,output=self.get('output'))
		exomu_inp, exomu_out = self.getCleanExoMu(map_spinp, map_spout, trouble)
		i = 0
		while not trouble.empty:
			map_spinp = gpyDB_wheels.adj.rc_pd(map_spinp, ('not', exomu_inp.droplevel(self.n('nn'))))
			map_spout = gpyDB_wheels.adj.rc_pd(map_spout, ('not', exomu_out.droplevel(self.n('n'))))
			trouble = self.TroubleNodes(map_spinp, map_spout)
			exomu_inp_i, exomu_out_i = self.getCleanExoMu(map_spinp,map_spout,trouble)
			exomu_inp, exomu_out = exomu_inp.union(exomu_inp_i), exomu_out.union(exomu_out_i)
			i += 1
			if i == maxiter:
				print("Algorithm for exomu sets failed; increase maxiter or consider nesting for loops.")
				break
		return exomu_inp, exomu_out

	@property
	def default_variables(self):
		return ('pS','pD','qS','qD','mu','sigma','eta','qnorm_out','qnorm_inp','qiv_out','qiv_inp')

	def initDB(self,m=None):
		return gpyDB_wheels.robust.robust_merge_dbs(self.s.db,{self.n(s): self.initSymbol(s,m=m) for s in self.default_variables},priority='first')

	def initSymbol(self,s,m=None):
		if s == 'pS':
			return gpyDB.gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('t'),self.get('output',m=m)],self.s.db), name=self.n(s)))
		elif s == 'qS':
			return gpyDB.gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('t'),self.get('output',m=m)],self.s.db), name=self.n(s)))
		elif s == 'pD':
			return gpyDB.gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('t'),self.get('int',m=m).union(self.get('input',m=m))],self.s.db), name = self.n(s)))
		elif s == 'qD':
			return gpyDB.gpy(pd.Series(0.5, index = adjMultiIndexDB.mergeDomains([self.get('t'),self.get('int',m=m).union(self.get('input',m=m))],self.s.db), name = self.n(s)))
		elif s == 'mu':
			return gpyDB.gpy(pd.Series(1, index = self.get('map',m=m), name = self.n(s)))
		elif s == 'eta':
			return gpyDB.gpy(pd.Series(0.5, index = self.get('knout',m=m), name = self.n(s)))
		elif s == 'sigma':
			return gpyDB.gpy(pd.Series(0.5, index = self.get('kninp',m=m), name = self.n(s)))
		elif s == 'qnorm_out':
			return gpyDB.gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('t'),self.get('knout',m=m)],self.s.db),name=self.n(s)),**{'type':'parameter'})
		elif s == 'qnorm_inp':
			return gpyDB.gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('t'),self.get('kninp',m=m)],self.s.db),name=self.n(s)),**{'type':'parameter'})
		elif s =='qiv_out':
			return gpyDB.gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('t'),self.get('spout',m=m)],self.s.db), name = self.n(s)))
		elif s =='qiv_inp':
			return gpyDB.gpy(pd.Series(1, index = adjMultiIndexDB.mergeDomains([self.get('t'),self.get('spinp',m=m)],self.s.db), name = self.n(s)))

	def groups(self,m=None):
		return {g.name: g for g in self.groups_(m=m)}
	def states(self,m=None):
		return {k: self.s.standardInstance(state=k) | {attr: getattr(self,attr)()[k] for attr in ('g_endo','g_exo','blocks','args')} for k in ('B','C')}
	def args(self):
		return {k: {self.name+'_Blocks': '\n'.join([getattr(GamsProduction, module.f)(self.name+'_'+name,name,inclusiveVal=False) for name,module in self.m.items()])} for k in ('B','C')}
	def blocks(self):
		return {k: OrdSet([f"B_{self.name}_{name}" for name in self.m]) for k in ('B','C')}
	def g_endo(self):
		return {'B': OrdSet([f"G_{self.name}_endo_always", f"G_{self.name}_exo_in_calib"]),
				'C': OrdSet([f"G_{self.name}_endo_always", f"G_{self.name}_endo_in_calib"])}
	def g_exo(self):
		return {'B': OrdSet([f"G_{self.name}_exo_always", f"G_{self.name}_endo_in_calib"]),
				'C': OrdSet([f"G_{self.name}_exo_always", f"G_{self.name}_exo_in_calib"])}
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
		v = 	[('mu', ('not', self.g('exomu',m=m))),
				 ('qS', ('and', [self.g('endo_qS',m=m), self.g('t0')]))])]
