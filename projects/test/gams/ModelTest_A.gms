$ONEOLCOM
$EOLCOM #


OPTION SYSOUT=OFF, SOLPRINT=OFF, LIMROW=0, LIMCOL=0, DECIMALS=6;


# ----------------------------------------------------------------------------------------------------
#  Define function: SolveEmptyNLP
# ----------------------------------------------------------------------------------------------------

sets
	alias_set
	alias_map2
	n
	s
;

alias(n,nn,nnn);

sets
	alias_[alias_set,alias_map2]
	t[t]
	t0[t]
	tE[t]
	tx0[t]
	txE[t]
	map_ModelTest_A[s,n,nn]
	map_spinp_ModelTest_A[s,n,nn]
	map_spout_ModelTest_A[s,n,nn]
	knout_ModelTest_A[s,n]
	kninp_ModelTest_A[s,n]
	spout_ModelTest_A[s,n]
	spinp_ModelTest_A[s,n]
	input_ModelTest_A[s,n]
	output_ModelTest_A[s,n]
	int_ModelTest_A[s,n]
	map_Tree1[s,n,nn]
	knot_Tree1[s,n]
	branch_Tree1[s,n]
	knot_o_Tree1[s,n]
	knot_no_Tree1[s,n]
	branch2o_Tree1[s,n]
	branch2no_Tree1[s,n]
	map_Tree2[s,n,nn]
	knot_Tree2[s,n]
	branch_Tree2[s,n]
	branch_o_Tree2[s,n]
	branch_no_Tree2[s,n]
	exomu_ModelTest_A[s,n,nn]
	endo_qD_ModelTest_A[s,n]
	endo_qS_ModelTest_A[s,n]
	endo_pS_ModelTest_A[s,n]
;
$GDXIN %gdxInput%
$onMulti
$load alias_set
$load alias_map2
$load n
$load s
$load alias_
$load t
$load t0
$load tE
$load tx0
$load txE
$load map_ModelTest_A
$load map_spinp_ModelTest_A
$load map_spout_ModelTest_A
$load knout_ModelTest_A
$load kninp_ModelTest_A
$load spout_ModelTest_A
$load spinp_ModelTest_A
$load input_ModelTest_A
$load output_ModelTest_A
$load int_ModelTest_A
$load map_Tree1
$load knot_Tree1
$load branch_Tree1
$load knot_o_Tree1
$load knot_no_Tree1
$load branch2o_Tree1
$load branch2no_Tree1
$load map_Tree2
$load knot_Tree2
$load branch_Tree2
$load branch_o_Tree2
$load branch_no_Tree2
$load exomu_ModelTest_A
$load endo_qD_ModelTest_A
$load endo_qS_ModelTest_A
$load endo_pS_ModelTest_A
$GDXIN
$offMulti;

variables
	R_LR
	g_LR
	infl_LR
	sigma[s,n]
	mu[s,n,nn]
	eta[s,n]
	pS[t,s,n]
	pD[t,s,n]
	qS[t,s,n]
	qD[t,s,n]
	qnorm_out[t,s,n]
	qnorm_inp[t,s,n]
	qiv_out[t,s,n]
	qiv_inp[t,s,n]
;
$GDXIN %gdxInput%
$onMulti
$load R_LR
$load g_LR
$load infl_LR
$load sigma
$load mu
$load eta
$load pS
$load pD
$load qS
$load qD
$load qnorm_out
$load qnorm_inp
$load qiv_out
$load qiv_inp
$GDXIN
$offMulti;




# -----------------------------------------B_ModelTest_A_Tree1----------------------------------------
#  Initialize B_ModelTest_A_Tree1 equation block
# ----------------------------------------------------------------------------------------------------
EQUATION E_zp_out_Tree1[t,s,n];
E_zp_out_Tree1[t,s,n]$(knot_o_tree1[s,n] and txe[t]).. 	pS[t,s,n]*qS[t,s,n]  =E=  sum(nn$(map_Tree1[s,n,nn]), qD[t,s,nn]*pD[t,s,nn]);
EQUATION E_zp_nout_Tree1[t,s,n];
E_zp_nout_Tree1[t,s,n]$(knot_no_tree1[s,n] and txe[t]).. 	pD[t,s,n]*qD[t,s,n]  =E=  sum(nn$(map_Tree1[s,n,nn]), qD[t,s,nn]*pD[t,s,nn]);
EQUATION E_q_out_Tree1[t,s,n];
E_q_out_Tree1[t,s,n]$(branch2o_tree1[s,n] and txe[t]).. 	qD[t,s,n]  =E=  sum(nn$(map_Tree1[s,nn,n]), mu[s,nn,n] * (pS[t,s,nn]/pD[t,s,n])**(sigma[s,nn]) * qS[t,s,nn]);
EQUATION E_q_nout_Tree1[t,s,n];
E_q_nout_Tree1[t,s,n]$(branch2no_tree1[s,n] and txe[t]).. 	qD[t,s,n]  =E=  sum(nn$(map_Tree1[s,nn,n]), mu[s,nn,n] * (pD[t,s,nn]/pD[t,s,n])**(sigma[s,nn]) * qD[t,s,nn]);

# ----------------------------------------------------------------------------------------------------
#  Define B_ModelTest_A_Tree1 model
# ----------------------------------------------------------------------------------------------------
Model B_ModelTest_A_Tree1 /
E_zp_out_Tree1, E_zp_nout_Tree1, E_q_out_Tree1, E_q_nout_Tree1
/;




# -----------------------------------------B_ModelTest_A_Tree2----------------------------------------
#  Initialize B_ModelTest_A_Tree2 equation block
# ----------------------------------------------------------------------------------------------------
EQUATION E_zp_Tree2[t,s,n];
E_zp_Tree2[t,s,n]$(knot_tree2[s,n] and txe[t]).. 	pD[t,s,n]*qD[t,s,n]  =E=  sum(nn$(map_Tree2[s,nn,n] and branch_o_Tree2[s,nn]), qS[t,s,nn]*pS[t,s,nn])+sum(nn$(map_Tree2[s,nn,n] and branch_no_Tree2[s,nn]), qD[t,s,nn]*pD[t,s,nn]);
EQUATION E_q_out_Tree2[t,s,n];
E_q_out_Tree2[t,s,n]$(branch_o_tree2[s,n] and txe[t]).. 	qS[t,s,n]  =E=  sum(nn$(map_Tree2[s,n,nn]), qD[t,s,nn]*mu[s,n,nn] * exp((pS[t,s,n]-pD[t,s,nn])*eta[s,nn]) / (sum(nnn$(map_Tree2[s,nnn,nn] and branch_o_Tree2[s,nnn]), mu[s,nnn,nn] * exp((pD[t,s,nn]-pS[t,s,nnn])*eta[s,nn]))+sum(nnn$(map_Tree2[s,nnn,nn] and branch_no_Tree2[s,nnn]), mu[s,nnn,nn] * exp((pD[t,s,nn]-pD[t,s,nnn])*eta[s,nn]))));
EQUATION E_q_nout_Tree2[t,s,n];
E_q_nout_Tree2[t,s,n]$(branch_no_tree2[s,n] and txe[t]).. 		qD[t,s,n]  =E=  sum(nn$(map_Tree2[s,n,nn]), qD[t,s,nn]*mu[s,n,nn] * exp((pD[t,s,n]-pD[t,s,nn])*eta[s,nn]) / (sum(nnn$(map_Tree2[s,nnn,nn] and branch_o_Tree2[s,nnn]), mu[s,nnn,nn] * exp((pD[t,s,nn]-pS[t,s,nnn])*eta[s,nn]))+sum(nnn$(map_Tree2[s,nnn,nn] and branch_no_Tree2[s,nnn]), mu[s,nnn,nn] * exp((pD[t,s,nn]-pD[t,s,nnn])*eta[s,nn]))));

# ----------------------------------------------------------------------------------------------------
#  Define B_ModelTest_A_Tree2 model
# ----------------------------------------------------------------------------------------------------
Model B_ModelTest_A_Tree2 /
E_zp_Tree2, E_q_out_Tree2, E_q_nout_Tree2
/;


qS.fx[t,s,n]$(((output_ModelTest_A[s,n] and ( not ((endo_qS_ModelTest_A[s,n] and t0[t])))) or (endo_qS_ModelTest_A[s,n] and t0[t]))) = qS.l[t,s,n];
pD.fx[t,s,n]$(input_ModelTest_A[s,n]) = pD.l[t,s,n];
sigma.fx[s,n]$(kninp_ModelTest_A[s,n]) = sigma.l[s,n];
eta.fx[s,n]$(knout_ModelTest_A[s,n]) = eta.l[s,n];
mu.fx[s,n,nn]$((exomu_ModelTest_A[s,n,nn] or ( not (exomu_ModelTest_A[s,n,nn])))) = mu.l[s,n,nn];
pD.lo[t,s,n]$(int_ModelTest_A[s,n]) = -inf;
pD.up[t,s,n]$(int_ModelTest_A[s,n]) = inf;
pS.lo[t,s,n]$((((output_ModelTest_A[s,n] and tx0[t]) or (endo_pS_ModelTest_A[s,n] and t0[t])) or ((output_ModelTest_A[s,n] and t0[t]) and ( not ((endo_pS_ModelTest_A[s,n] and t0[t])))))) = -inf;
pS.up[t,s,n]$((((output_ModelTest_A[s,n] and tx0[t]) or (endo_pS_ModelTest_A[s,n] and t0[t])) or ((output_ModelTest_A[s,n] and t0[t]) and ( not ((endo_pS_ModelTest_A[s,n] and t0[t])))))) = inf;
qD.lo[t,s,n]$(((((int_ModelTest_A[s,n] or input_ModelTest_A[s,n]) and tx0[t]) or (endo_qD_ModelTest_A[s,n] and t0[t])) or (((int_ModelTest_A[s,n] or input_ModelTest_A[s,n]) and t0[t]) and ( not ((endo_qD_ModelTest_A[s,n] and t0[t])))))) = -inf;
qD.up[t,s,n]$(((((int_ModelTest_A[s,n] or input_ModelTest_A[s,n]) and tx0[t]) or (endo_qD_ModelTest_A[s,n] and t0[t])) or (((int_ModelTest_A[s,n] or input_ModelTest_A[s,n]) and t0[t]) and ( not ((endo_qD_ModelTest_A[s,n] and t0[t])))))) = inf;
qiv_inp.lo[t,s,n]$(spinp_ModelTest_A[s,n]) = -inf;
qiv_inp.up[t,s,n]$(spinp_ModelTest_A[s,n]) = inf;
qiv_out.lo[t,s,n]$(spout_ModelTest_A[s,n]) = -inf;
qiv_out.up[t,s,n]$(spout_ModelTest_A[s,n]) = inf;

# ----------------------------------------------------------------------------------------------------
#  Define ModelTest_A_B model
# ----------------------------------------------------------------------------------------------------
Model ModelTest_A_B /
E_zp_out_Tree1, E_zp_nout_Tree1, E_q_out_Tree1, E_q_nout_Tree1, E_zp_Tree2, E_q_out_Tree2, E_q_nout_Tree2
/;


solve ModelTest_A_B using CNS;