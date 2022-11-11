# 0.1: Auxiliary functions used for scale-preserving nests:
def _CES(px,py,mu,sigma,norm=None):
	return f"{mu} * ({py}/{px})**({sigma})" if norm is None else f"{mu} * ({py}/({px}*(1+{norm})))**({sigma})"
def _exp(px,py,mu,sigma,norm=None):
	return f"{mu} * exp(({py}-{px})*{sigma})" if norm is None else f"{mu} * exp(({py}-{px})*{sigma}-{norm})"

# 0.2: Normalized demand
def _Fnorm_input_demand(ftype, name):
	f = globals()['_'+ftype]
	return f"""E_q_{name}[t,s,n]$(branch_{name}[s,n] and txE[t])..		qD[t,s,n] =E= sum(nn$(map_{name}[s,nn,n]), qS[t,s,nn]*{f("pD[t,s,n]","pD[t,s,nn]","mu[s,nn,n]","sigma[s,nn]")} / sum(nnn$(map_{name}[s,nn,nnn]), {f("pD[t,s,nnn]","pD[t,s,nn]","mu[s,nn,nnn]","sigma[s,nn]")}));"""

def _Fnorm_output_demand(ftype, name):
	f = globals()['_'+ftype]
	return f"""E_q_{name}[t,s,n]$(branch_{name}[s,n] and txE[t])..		qD[t,s,n] =E= sum(nn$(map_{name}[s,n,nn]), qD[t,s,nn]*{f("pD[t,s,nn]","pD[t,s,n]","mu[s,n,nn]","eta[s,nn]")} / sum(nnn$(map_{name}[s,nnn,nn]), {f("pD[t,s,nnn]","pD[t,s,nn]","mu[s,nnn,nn]","eta[s,nn]")}));"""

# 0.2: Zero profit equations:
def zp_input(name):
	return f"""E_zp_{name}[t,s,n]$(knot_{name}[s,n] and txE[t])..	pD[t,s,n]*qD[t,s,n] =E= sum(nn$(map_{name}[s,n,nn]), qD[t,s,nn]*pD[t,s,nn]);"""
def zp_output(name):
	return f"""E_zp_{name}[t,s,n]$(knot_{name}[s,n] and txE[t])..	pD[t,s,n]*qD[t,s,n] =E= sum(nn$(map_{name}[s,nn,n]), qD[t,s,nn]*pD[t,s,nn]);"""

# 1: Input type nests:
# 1.1: CES nest:
def CES(name,m,**kwargs):
	return f"""
$BLOCK B_{name}
	{zp_input(name)}
	E_q_{name}[t,s,n]$(branch_{name}[s,n] and txE[t])..	qD[t,s,n] =E= sum(nn$(map_{name}[s,nn,n]), mu[s,nn,n] * (pD[t,s,nn]/pD[t,s,n])**(sigma[s,nn]) * qD[t,s,nn]);
$ENDBLOCK
"""

# 1.2: Scale-preserving nests:
def Fnorm_input(ftype,name,m):
	return f"""
$BLOCK B_{name}
	{zp_input(name)}
	{_Fnorm_input_demand(ftype,name)}
$ENDBLOCK
"""
def CES_norm(name,m):
	return Fnorm_input('CES',name,m)
def MNL(name,m):
	return Fnorm_input('exp',name,m)

# 2: Output type nests:
# 2.1: CET function:
def CET(name,m,**kwargs):
	return f"""
$BLOCK B_{name}
	{zp_output(name)}
	E_demand_{name}[t,s,n]$(branch_{name}[s,n] and txE[t])..	qD[t,s,n] =E= sum(nn$(map_{name}[s,n,nn]), mu[s,n,nn] * (pD[t,s,n]/pD[t,s,nn])**(eta[s,nn]) * qD[t,s,nn]);
$ENDBLOCK
"""

# 2.2: scale-preserving nests: 
def Fnorm_output(ftype,name,m):
	return f"""
$BLOCK B_{name}
	{zp_output(name)}
	{_Fnorm_output_demand(ftype,name)}
$ENDBLOCK
"""

def CET_norm(name,m):
	return Fnorm_output('CES',name,m)
def MNL_out(name,m):
	return Fnorm_output('exp',name,m)

# 3: Labor supply function
def IsoFrisch(name, m):
	return f"""
$BLOCK B_{name}
	E_labor_{name}[t,s,n]$(labor_{m}[s,n] and txE[t])..	qS[t,s,n]	=E=	Lscale[s,n] * ( sum(nn$(L2C_{m}[s,n,nn]), pS[t,s,n]/(pD[t,s,nn]*(qD[t,s,nn]**(crra[s,nn]))) )**(frisch[s,n]) );
$ENDBLOCK
"""

# 4: Ramsey:
def simpleDynamic(name, m):
	return f"""
$BLOCK B_{name}
	E_lom_{name}[t,s,n]$(svngs_{m}[s,n] and txE[t])..		vD[t+1,s,n] =E= (vD[t,s,n]*iRate[t]+sp[t,s])/((1+g_LR)*(1+infl_LR));
	E_euler_{name}[t,s,n]$(output_{m}[s,n] and tx0E[t])..	qD[t,s,n]	=E= qD[t-1,s,n]*( (disc[s]*iRate[t]*pD[t-1,s,n]/(pD[t,s,n]*(1+infl_LR)))**(1/crra[s,n]) )/ (1+g_LR);
	E_tvc_{name}[t,s,n]$(svngs_{m}[s,n] and tE[t])..		vD[t,s,n]	=E= (1+h_tvc[s,n])*vD[t-1,s,n];
$ENDBLOCK
"""

def priceWedge(name,m):
	return f"""
$BLOCK B_{name}
	E_pw_{name}[t,s,n]$(labor_{m}[s,n] and txE[t])..	pS[t,s,n] 		=E= p[t,n]-tauS[t,s,n];
	E_TaxRev_{name}[t,s]$(s_{m}[s] and txE[t])..		TotalTax[t,s]	=E= tauLump[t,s]+sum(n$(input_{m}[s,n]), tauD[t,s,n] * qD[t,s,n])+sum(n$(labor_{m}[s,n]), tauS[t,s,n]*qS[t,s,n]);
	E_sp_{name}[t,s]$(s_{m}[s] and txE[t])..			sp[t,s]			=E= sum(n$(labor_{m}[s,n]), pS[t,s,n]*qS[t,s,n]) - sum(n$(input_{m}[s,n]), pD[t,s,n]*qD[t,s,n])-tauLump[t,s];
$ENDBLOCK
"""


# 6: System of value shares
def valueShares():
	return f"""
$BLOCK B_ValueShares
	E_Out_knot[t,s,n]$(knotOutTree[s,n])..								vD[t,s,n]	=E= sum(nn$(map[s,nn,n] and branchOut[s,nn]), vD[t,s,n])+sum(nn$(map[s,nn,n] and branchNOut[s,nn]), vD[t,s,n]);
	E_Out_shares_o[t,s,n,nn]$(mapOut[s,n,nn] and branchOut[s,n])..		mu[t,s,n,nn]=E= vD[t,s,n]/vD[t,s,nn];
	E_Out_shares_no[t,s,n,nn]$(mapOut[s,n,nn] and branchNOut[s,n])..	mu[t,s,n,nn]=E= vD[t,s,n]/vD[t,s,nn];
	E_Inp_knot_o[t,s,n]$(knotOut[s,n])..								vD[t,s,n]	=E= sum(nn$(map[s,n,nn]), vD[t,s,nn]);
	E_Inp_knot_no[t,s,n]$(knotNOut[s,n])..								vD[t,s,n]	=E= sum(nn$(map[s,n,nn]), vD[t,s,nn]);
	E_Inp_shares2o[t,s,n,nn]$(mapInp[s,n,nn] and branch2Out[s,nn])..	mu[t,s,n,nn]=E= vD[t,s,nn]/vD[t,s,n];
	E_Inp_shares2no[t,s,n,nn]$(mapInp[s,n,nn] and branch2NOut[s,nn])..	mu[t,s,n,nn]=E= vD[t,s,nn]/vD[t,s,n];
$ENDBLOCK
"""
