# 0.1: Auxiliary functions used for scale-preserving nests:
def _CES(px,py,mu,sigma,norm=None):
	return f"{mu} * ({py}/{px})**({sigma})" if norm is None else f"{mu} * ({py}/({px}*(1+{norm})))**({sigma})"
def _exp(px,py,mu,sigma,norm=None):
	return f"{mu} * exp(({py}-{px})*{sigma})" if norm is None else f"{mu} * exp(({py}-{px})*{sigma}-{norm})"

def _Fnorm_input_demand(ftype,name):
	f = globals()['_'+ftype]
	return f"""E_q_out_{name}[t,s,n]$(branch2o_{name}[s,n] and txE[t])..	qD[t,s,n] =E= sum(nn$(map_{name}[s,nn,n]), qS[t,s,nn]*{f("pD[t,s,n]","pS[t,s,nn]","mu[s,nn,n]","sigma[s,nn]")} / sum(nnn$(map_{name}[s,nn,nnn]), {f("pD[t,s,nnn]","pS[t,s,nn]","mu[s,nn,nnn]","sigma[s,nn]")}));
	E_q_nout_{name}[t,s,n]$(branch2no_{name}[s,n] and txE[t])..		qD[t,s,n] =E= sum(nn$(map_{name}[s,nn,n]), qD[t,s,nn]*{f("pD[t,s,n]","pD[t,s,nn]","mu[s,nn,n]","sigma[s,nn]")} / sum(nnn$(map_{name}[s,nn,nnn]), {f("pD[t,s,nnn]","pD[t,s,nn]","mu[s,nn,nnn]","sigma[s,nn]")}));"""

def _Fnorm_output_demand(ftype,name):
	f = globals()['_'+ftype]
	return f"""E_q_out_{name}[t,s,n]$(branch_o_{name}[s,n] and txE[t])..	qS[t,s,n] =E= sum(nn$(map_{name}[s,n,nn]), qD[t,s,nn]*{f("pD[t,s,nn]","pS[t,s,n]","mu[s,n,nn]","eta[s,nn]")} / (sum(nnn$(map_{name}[s,nnn,nn] and branch_o_{name}[s,nnn]), {f("pS[t,s,nnn]","pD[t,s,nn]","mu[s,nnn,nn]","eta[s,nn]")})+sum(nnn$(map_{name}[s,nnn,nn] and branch_no_{name}[s,nnn]), {f("pD[t,s,nnn]","pD[t,s,nn]","mu[s,nnn,nn]","eta[s,nn]")})));
	E_q_nout_{name}[t,s,n]$(branch_no_{name}[s,n] and txE[t])..		qD[t,s,n] =E= sum(nn$(map_{name}[s,n,nn]), qD[t,s,nn]*{f("pD[t,s,nn]","pD[t,s,n]","mu[s,n,nn]","eta[s,nn]")} / (sum(nnn$(map_{name}[s,nnn,nn] and branch_o_{name}[s,nnn]), {f("pS[t,s,nnn]","pD[t,s,nn]","mu[s,nnn,nn]","eta[s,nn]")})+sum(nnn$(map_{name}[s,nnn,nn] and branch_no_{name}[s,nnn]), {f("pD[t,s,nnn]","pD[t,s,nn]","mu[s,nnn,nn]","eta[s,nn]")})));"""

def _Fnorm_input_with_InclusiveValue(ftype,name):
	f = globals()['_'+ftype]
	return f"""E_q_out_{name}[t,s,n]$(branch2o_{name}[s,n] and txE[t])..	qD[t,s,n] =E= sum(nn$(map_{name}[s,nn,n]), qS[t,s,nn]*{f("pD[t,s,n]","pS[t,s,nn]","mu[s,nn,n]","sigma[s,nn]",norm='qnorm[t,s,nn]')} / qiv_inp[t,s,nn]);
	E_q_nout_{name}[t,s,n]$(branch2no_{name}[s,n] and txE[t])..		qD[t,s,n] =E= sum(nn$(map_{name}[s,nn,n]), qD[t,s,nn]*{f("pD[t,s,n]","pD[t,s,nn]","mu[s,nn,n]","sigma[s,nn]",norm='qnorm[t,s,nn]')} / qiv_inp[t,s,nn]);
	E_inclVal_out_{name}[t,s,n]$(knot_o_{name}[s,n] and txE[t])..	qiv_inp[t,s,n] =E= sum(nn$(map_{name}[s,n,nn]), {f("pD[t,s,nn]","pS[t,s,n]","mu[s,n,nn]","sigma[s,n]",norm='qnorm[t,s,n]')});
	E_inclVal_nout_{name}[t,s,n]$(knot_no_{name}[s,n] and txE[t])..	qiv_inp[t,s,n] =E= sum(nn$(map_{name}[s,n,nn]), {f("pD[t,s,nn]","pD[t,s,n]","mu[s,n,nn]","sigma[s,n]",norm='qnorm[t,s,n]')});"""

def _Fnorm_output_with_InclusiveValue(ftype,name):
	f = globals()['_'+ftype]
	return f"""E_q_out_{name}[t,s,n]$(branch_o_{name}[s,n] and txE[t])..	qS[t,s,n] =E= sum(nn$(map_{name}[s,n,nn]), qD[t,s,nn]*{f("pD[t,s,nn]","pS[t,s,n]","mu[s,n,nn]","eta[s,nn]")} / qiv_out[t,s,nn]);
	E_q_nout_{name}[t,s,n]$(branch_no_{name}[s,n] and txE[t])..		qD[t,s,n] =E= sum(nn$(map_{name}[s,n,nn]), qD[t,s,nn]*{f("pD[t,s,nn]","pD[t,s,n]","mu[s,n,nn]","eta[s,nn]")} / qiv_out[t,s,nn]);
	E_inclVal_out_{name}[t,s,n]$(knot_{name}[s,n] and txE[t])..		qiv_out[t,s,n]=E= sum(nn$(map_{name}[s,nn,n] and branch_o_{name}[s,nn]), {f("pS[t,s,nn]","pD[t,s,n]","mu[s,nn,n]","eta[s,n]",norm='qnorm[t,s,n]')})+sum(nn$(map_{name}[s,nn,n] and branch_no_{name}[s,nn]), {f("pD[t,s,nn]","pD[t,s,n]","mu[s,nn,n]","eta[s,n]",norm='qnorm[t,s,n]')});"""

# 0.2: Zero profit equations:
def zp_input(name):
	return f"""E_zp_out_{name}[t,s,n]$(knot_o_{name}[s,n] and txE[t])..	pS[t,s,n]*qS[t,s,n] =E= sum(nn$(map_{name}[s,n,nn]), qD[t,s,nn]*pD[t,s,nn]);
	E_zp_nout_{name}[t,s,n]$(knot_no_{name}[s,n] and txE[t])..	pD[t,s,n]*qD[t,s,n] =E= sum(nn$(map_{name}[s,n,nn]), qD[t,s,nn]*pD[t,s,nn]);"""
def zp_output(name):
	return f"""E_zp_{name}[t,s,n]$(knot_{name}[s,n] and txE[t])..	pD[t,s,n]*qD[t,s,n] =E= sum(nn$(map_{name}[s,nn,n] and branch_o_{name}[s,nn]), qS[t,s,nn]*pS[t,s,nn])+sum(nn$(map_{name}[s,nn,n] and branch_no_{name}[s,nn]), qD[t,s,nn]*pD[t,s,nn]);"""


# 1: Input type nests:
# 1.1: CES nest:
def CES(name,m,**kwargs):
	return f"""
$BLOCK B_{name}
	{zp_input(name)}
	E_q_out_{name}[t,s,n]$(branch2o_{name}[s,n] and txE[t])..	qD[t,s,n] =E= sum(nn$(map_{name}[s,nn,n]), mu[s,nn,n] * (pS[t,s,nn]/pD[t,s,n])**(sigma[s,nn]) * qS[t,s,nn]);
	E_q_nout_{name}[t,s,n]$(branch2no_{name}[s,n] and txE[t])..	qD[t,s,n] =E= sum(nn$(map_{name}[s,nn,n]), mu[s,nn,n] * (pD[t,s,nn]/pD[t,s,n])**(sigma[s,nn]) * qD[t,s,nn]);
$ENDBLOCK
"""

# 1.2: Scale-preserving nests:
def Fnorm_input(ftype,name,m,inclusiveVal = False):
	return f"""
$BLOCK B_{name}
	{zp_input(name)}
	{_Fnorm_input_with_InclusiveValue(ftype,name) if inclusiveVal else _Fnorm_input_demand(ftype,name)}
$ENDBLOCK
"""
def CES_norm(name,m,inclusiveVal = False):
	return Fnorm_input('CES',name,m,inclusiveVal=inclusiveVal)
def MNL(name,m,inclusiveVal = False):
	return Fnorm_input('exp',name,m,inclusiveVal=inclusiveVal)

# 2: Output type nests:
# 2.1: CET function:
def CET(name,m,**kwargs):
	return f"""
$BLOCK B_{name}
	{zp_output(name)}
	E_demand_out_{name}[t,s,n]$(branch_o_{name}[s,n] and txE[t])..		qS[t,s,n] =E= sum(nn$(map_{name}[s,n,nn]), mu[s,n,nn] * (pS[t,s,n]/pD[t,s,nn])**(eta[s,nn]) * qD[t,s,nn]);
	E_demand_nout_{name}[t,s,n]$(branch_no_{name}[s,n] and txE[t])..	qD[t,s,n] =E= sum(nn$(map_{name}[s,n,nn]), mu[s,n,nn] * (pD[t,s,n]/pD[t,s,nn])**(eta[s,nn]) * qD[t,s,nn]);
$ENDBLOCK
"""
# 2.2: scale-preserving nests: 
def Fnorm_output(ftype,name,m,inclusiveVal = False):
	return f"""
$BLOCK B_{name}
	{zp_output(name)}
	{_Fnorm_output_with_InclusiveValue(ftype,name) if inclusiveVal else _Fnorm_output_demand(ftype,name)}
$ENDBLOCK
"""

def CET_norm(name,m,inclusiveVal = False):
	return Fnorm_output('CES',name,m,inclusiveVal=inclusiveVal)
def MNL_out(name,m,inclusiveVal = False):
	return Fnorm_output('exp',name,m,inclusiveVal=inclusiveVal)

# 3: Adjustment costs / installation cost equations:
def sqrAdjCosts(name, m):
	return f"""
$BLOCK B_{name}
	E_lom_{name}[t,s,n]$(dur_{m}[s,n] and txE[t])..		qD[t+1,s,n]	=E= (qD[t,s,n]*(1-rDepr[t,s,n])+sum(nn$(dur2inv[s,n,nn]), qD[t,s,nn]))/(1+g_LR);
	E_pk_{name}[t,s,n]$(dur_{m}[s,n] and tx02E[t])..	pD[t,s,n]	=E= sqrt(sqr(sum(nn$(dur2inv[s,n,nn]), Rrate[t]*pD[t-1,s,nn]*(1+icpar[s,n]*(qD[t-1,s,nn]/qD[t-1,s,n]-(rDepr[t-1,s,n]+g_LR)))/(1+infl_LR)+pD[t,s,nn]*(icpar[s,n]*0.5*(sqr(rDepr[t,s,n]+g_LR)-sqr(qD[t,s,nn]/qD[t,s,n]))-(1-rDepr[t,s,n])*(1+icpar[s,n]*(qD[t,s,nn]/qD[t,s,n]-(rDepr[t,s,n]+g_LR)))))));
	E_pkT_{name}[t,s,n]$(dur_{m}[s,n] and t2E[t])..		pD[t,s,n]	=E= sum(nn$(dur2inv[s,n,nn]), Rrate[t]*pD[t-1,s,nn] * (1+icpar[s,n]*(qD[t-1,s,nn]/qD[t-1,s,n]-(rDepr[t-1,s,n]+g_LR)))/(1+infl_LR) + (rDepr[t,s,n]-1)*pD[t,s,nn]);
	E_Ktvc_{name}[t,s,n]$(dur_{m}[s,n] and tE[t])..		qD[t,s,n]	=E= (1+K_tvc[s,n])*qD[t-1,s,n];
	E_instcost_{name}[t,s]$(s_{m}[s] and txE[t])..		ic[t,s] 	=E= sum([n,nn]$(dur2inv[s,n,nn]), pD[t,s,nn] * icpar[s,n]*0.5*qD[t,s,n]*sqr(qD[t,s,nn]/qD[t,s,n]-(rDepr[t,s,n]+g_LR)));
$ENDBLOCK
"""


# 4: Introduce price wedge with mark-up, unit-tax, and installation costs
def priceWedge(name,m):
	return f"""
$BLOCK B_{name}
	E_pwInp_{name}[t,s,n]$(input_{m}[s,n] and txE[t])..			pD[t,s,n]		=E= p[t,n]+tauD[t,s,n];	
	E_pwOut_{name}[t,s,n]$(output_{m}[s,n] and txE[t])..		p[t,n] 			=E= (1+markup[s])*(pS[t,s,n]+tauS[t,s,n]+(outShare[t,s,n]/qS[t,s,n])*(ic[t,s]+tauLump[t,s]));
	E_outShare_{name}[t,s,n]$(output_{m}[s,n] and txE[t])..		outShare[t,s,n] =E= qS[t,s,n]*pS[t,s,n]/(sum(nn$(output_{m}[s,nn]), qS[t,s,nn]*pS[t,s,nn]));
	E_TaxRev_{name}[t,s]$(s_{m}[s] and txE[t])..				TotalTax[t,s]	=E= tauLump[t,s]+sum(n$(input_{m}[s,n]), tauD[t,s,n] * qD[t,s,n])+sum(n$(output_{m}[s,n]), tauS[t,s,n]*qS[t,s,n]);
$ENDBLOCK
"""


# 5: System of value shares
def valueShares():
	return f"""
$BLOCK B_ValueShares
	E_Out_knot[t,s,n]$(knotOutTree[s,n])..								vD[t,s,n]	=E= sum(nn$(map[s,nn,n] and branchOut[s,nn]), vS[t,s,n])+sum(nn$(map[s,nn,n] and branchNOut[s,nn]), vD[t,s,n]);
	E_Out_shares_o[t,s,n,nn]$(mapOut[s,n,nn] and branchOut[s,n])..		mu[t,s,n,nn]=E= vS[t,s,n]/vD[t,s,nn];
	E_Out_shares_no[t,s,n,nn]$(mapOut[s,n,nn] and branchNOut[s,n])..	mu[t,s,n,nn]=E= vD[t,s,n]/vD[t,s,nn];
	E_Inp_knot_o[t,s,n]$(knotOut[s,n])..								vS[t,s,n]	=E= sum(nn$(map[s,n,nn]), vD[t,s,nn]);
	E_Inp_knot_no[t,s,n]$(knotNOut[s,n])..								vD[t,s,n]	=E= sum(nn$(map[s,n,nn]), vD[t,s,nn]);
	E_Inp_shares2o[t,s,n,nn]$(mapInp[s,n,nn] and branch2Out[s,nn])..	mu[t,s,n,nn]=E= vD[t,s,nn]/vS[t,s,n];
	E_Inp_shares2no[t,s,n,nn]$(mapInp[s,n,nn] and branch2NOut[s,nn])..	mu[t,s,n,nn]=E= vD[t,s,nn]/vD[t,s,n];
$ENDBLOCK
"""
