def BalancedBudget(name):
	return f"""
$BLOCK B_{name}
	E_bb_{name}[t,s]$(s_G[s] and tx0E[t])..		sum(ss$(s_Tax[ss]), TotalTax[t,ss])+sum(n$(gsvngs[s,n]), vD[t,s,n]*(iRate[t]-1)) = 0;
$ENDBLOCK
"""
