$BLOCK B_ModelTest_B_Tree1
	E_zp_out_Tree1[t,s,n]$(knot_o_Tree1[s,n] and txE[t])..	pS[t,s,n]*qS[t,s,n] =E= sum(nn$(map_Tree1[s,n,nn]), qD[t,s,nn]*pD[t,s,nn]);
	E_zp_nout_Tree1[t,s,n]$(knot_no_Tree1[s,n] and txE[t])..	pD[t,s,n]*qD[t,s,n] =E= sum(nn$(map_Tree1[s,n,nn]), qD[t,s,nn]*pD[t,s,nn]);
	E_q_out_Tree1[t,s,n]$(branch2o_Tree1[s,n] and txE[t])..	qD[t,s,n] =E= sum(nn$(map_Tree1[s,nn,n]), mu[s,nn,n] * (pS[t,s,nn]/pD[t,s,n])**(sigma[s,nn]) * qS[t,s,nn]);
	E_q_nout_Tree1[t,s,n]$(branch2no_Tree1[s,n] and txE[t])..	qD[t,s,n] =E= sum(nn$(map_Tree1[s,nn,n]), mu[s,nn,n] * (pD[t,s,nn]/pD[t,s,n])**(sigma[s,nn]) * qD[t,s,nn]);
$ENDBLOCK


$BLOCK B_ModelTest_B_Tree2
	E_zp_Tree2[t,s,n]$(knot_Tree2[s,n] and txE[t])..	pD[t,s,n]*qD[t,s,n] =E= sum(nn$(map_Tree2[s,nn,n] and branch_o_Tree2[s,nn]), qS[t,s,nn]*pS[t,s,nn])+sum(nn$(map_Tree2[s,nn,n] and branch_no_Tree2[s,nn]), qD[t,s,nn]*pD[t,s,nn]);
	E_q_out_Tree2[t,s,n]$(branch_o_Tree2[s,n] and txE[t])..	qS[t,s,n] =E= sum(nn$(map_Tree2[s,n,nn]), qD[t,s,nn]*mu[s,n,nn] * exp((pS[t,s,n]-pD[t,s,nn])*eta[s,nn]) / (sum(nnn$(map_Tree2[s,nnn,nn] and branch_o_Tree2[s,nnn]), mu[s,nnn,nn] * exp((pD[t,s,nn]-pS[t,s,nnn])*eta[s,nn]))+sum(nnn$(map_Tree2[s,nnn,nn] and branch_no_Tree2[s,nnn]), mu[s,nnn,nn] * exp((pD[t,s,nn]-pD[t,s,nnn])*eta[s,nn]))));
	E_q_nout_Tree2[t,s,n]$(branch_no_Tree2[s,n] and txE[t])..		qD[t,s,n] =E= sum(nn$(map_Tree2[s,n,nn]), qD[t,s,nn]*mu[s,n,nn] * exp((pD[t,s,n]-pD[t,s,nn])*eta[s,nn]) / (sum(nnn$(map_Tree2[s,nnn,nn] and branch_o_Tree2[s,nnn]), mu[s,nnn,nn] * exp((pD[t,s,nn]-pS[t,s,nnn])*eta[s,nn]))+sum(nnn$(map_Tree2[s,nnn,nn] and branch_no_Tree2[s,nnn]), mu[s,nnn,nn] * exp((pD[t,s,nn]-pD[t,s,nnn])*eta[s,nn]))));
$ENDBLOCK
