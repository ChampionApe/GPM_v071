from CGE_GmsPython import *

class Equi(GmsPython):
	def __init__(self, f = None, name = None, db_IO = None, itory = None, s=None, glob=None, s_kwargs = None, g_kwargs = None):
		""" Initialize from name, io data, and subset of inventory"""
		super().__init__(name=f"{name}_equi", f=f, s=s, glob=glob, g_kwargs=g_kwargs, s_kwargs = noneInit(s_kwargs, {}) | {'db': db_IO})

	def states(self,m=None):
		return {k: self.s.standardInstance(state=k) | {attr: getattr(self,attr)()[k] for attr in ('g_endo','g_exo','blocks','args')} for k in ('B','C')}
	def args(self):
		return {k: {f"{self.name}_Blocks": self.equationText} for k in ('B','C')}
	def blocks(self):
		return {'B': OrdSet([f"B_{self.name}"]), 'C': OrdSet()}
	def g_endo(self):
		return {'B': OrdSet([f"G_{self.name}_endo"]), 'C': OrdSet()}
	def g_exo(self):
		return {'B': OrdSet(), 'C': OrdSet()}
	def groups(self,m=None):
		return {g.name: g for g in [GmsPy.Group(f"G_{self.name}_endo",  v = [('qS', self.g('d_qSEqui')), ('p', self.g('d_pEqui'))])]}
	@property
	def equationText(self):
		return f"""
$BLOCK B_{self.name}
	E_equi_{self.name}[t,n]$(nEqui[n] and txE[t])..	 sum(s$(d_qS[s,n]), qS[t,s,n]) =E= sum(s$(d_qD[s,n]), qD[t,s,n]);
$ENDBLOCK
"""
