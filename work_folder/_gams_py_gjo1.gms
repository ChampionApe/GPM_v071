$ONEOLCOM
$EOLCOM #
;
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


sets
	alias_[alias_set,alias_map2]
	active[n,s]
	activeRow[n]
	activeCol[s]
;
$GDXIN %rname%
$onMulti
$load alias_set
$load alias_map2
$load n
$load s
$load alias_
$load active
$load activeRow
$load activeCol
$GDXIN
$offMulti;

parameters
	vD0[n,s]
	deltaRow[n]
	deltaCol[s]
	rowSum[n]
	colSum[s]
;
$GDXIN %rname%
$onMulti
$load vD0
$load deltaRow
$load deltaCol
$load rowSum
$load colSum
$GDXIN
$offMulti;

variables
	vD[n,s]
	etaRow[n,s]
	etaCol[n,s]
	object
;
$GDXIN %rname%
$onMulti
$load vD
$load etaRow
$load etaCol
$load object
$GDXIN
$offMulti;




# ------------------------------------------------B_RAS-----------------------------------------------
#  Initialize B_RAS equation block
# ----------------------------------------------------------------------------------------------------
EQUATION E_object;
E_object.. 							object	  =E=  sum([n,s]$(active[n,s]), Sqr(etaRow[n,s]-1)+Sqr(etaCol[n,s]-1));
EQUATION E_vD[n,s];
E_vD[n,s]$(active[n,s]).. 	vD[n,s]	  =E=  vD0[n,s]*(1-etaRow[n,s]*deltaRow[n]-etaCol[n,s]*deltaCol[s]);
EQUATION E_colSum[s];
E_colSum[s]$(activecol[s]).. 	colSum[s] =E=  sum(n$(active[n,s]), vD[n,s]);
EQUATION E_rowSum[n];
E_rowSum[n]$(activerow[n]).. 	rowSum[n] =E=  sum(s$(active[n,s]), vD[n,s]);

# ----------------------------------------------------------------------------------------------------
#  Define B_RAS model
# ----------------------------------------------------------------------------------------------------
Model B_RAS /
E_object, E_vD, E_colSum, E_rowSum
/;



vD.lo[n,s]$(active[n,s]) = -inf;
vD.up[n,s]$(active[n,s]) = inf;
etaRow.lo[n,s]$(active[n,s]) = -inf;
etaRow.up[n,s]$(active[n,s]) = inf;
etaCol.lo[n,s]$(active[n,s]) = -inf;
etaCol.up[n,s]$(active[n,s]) = inf;

# ----------------------------------------------------------------------------------------------------
#  Define someName_B model
# ----------------------------------------------------------------------------------------------------
Model someName_B /
E_object, E_vD, E_colSum, E_rowSum
/;


vD.lo[n,s]$(active[n,s]) = 0; solve someName_B minimizing object using QCP;