from CGE_GmsPython import *
import GamsGovernment

class BalancedBudget(GmsPython):
	def __init__(self, f = None, name = None, ns = None, s = None, glob = None, s_kwargs = None, g_kwargs = None, kwargs=None):
		super().__init__(name=name, f=f, s=s, glob=glob, ns=ns, s_kwargs = s_kwargs, g_kwargs=g_kwargs)

	def initDB(self,m=None):
		return gpyDB_wheels.robust.robust_merge_dbs(self.s.db,self.initSymbols(),priority='first')

	def initSymbols(self):
		return {'TotalTax': gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('txE'), self.get('s_Tax')], self.s.db), name = self.n('TotalTax'))),
				'vD': gpy(pd.Series(0, index = adjMultiIndexDB.mergeDomains([self.get('t'), self.get('gsvngs')], self.s.db), name = self.n('vD'))),
				'iRate': gpy(pd.Series(self.get('R_LR')*(1+self.get('infl_LR')), index = self.get('t'), name = self.n('iRate')))}

	def states(self, m= None):
		return {'B': self.s.standardInstance(state='B') | {attr: getattr(self,attr)() for attr in ('g_endo','g_exo','blocks','args')}}
	def args(self):
		return {f"{self.name}_Blocks": GamsGovernment.BalancedBudget(self.name)}
	def blocks(self):
		return OrdSet([f"B_{self.name}"])
	def g_endo(self):
		return OrdSet([f"G_{self.name}_endo"])
	def g_exo(self):
		return OrdSet([f"G_{self.name}_exo"])

# NOTE: Still missing groups here. 