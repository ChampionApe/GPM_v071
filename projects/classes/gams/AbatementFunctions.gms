# Root File for model
OPTION SYSOUT=OFF, SOLPRINT=OFF, LIMROW=0, LIMCOL=0, DECIMALS=6;
$SETLOCAL qmark ";


$FUNCTION load_level({group}, {gdx}):
  $offlisting
  $GROUP __load_group {group};
  $LOOP __load_group:
    parameter load_{name}{sets} "";
    load_{name}{sets}$({conditions}) = 0;
  $ENDLOOP
  execute_load {gdx} $LOOP __load_group: load_{name}={name}.l $ENDLOOP;
  $LOOP __load_group:
    {name}.l{sets}$({conditions}) = load_{name}{sets};
  $ENDLOOP
  $onlisting
$ENDFUNCTION

$FUNCTION load_fixed({group}, {gdx}):
  $offlisting
  $GROUP __load_group {group};
  $LOOP __load_group:
    parameter load_{name}{sets} "";
    load_{name}{sets}$({conditions}) = 0;
  $ENDLOOP
  execute_load {gdx} $LOOP __load_group: load_{name}={name}.l $ENDLOOP;
  $LOOP __load_group:
    {name}.fx{sets}$({conditions}) = load_{name}{sets};
  $ENDLOOP
  $onlisting
$ENDFUNCTION

$MACRO TechLogit(p,lambda_,sigma_,mu,map_) sum(nn$(map_), (1/(1+exp((p-mu-lambda_)/sigma_))))
$MACRO TechNormal(p,lambda_,sigma_,mu,map_) sum(nn$(map_),errorf((lambda_-p+mu)/sigma_))
$MACRO TechLognormal(p,lambda_,sigma_,mu,map_) sum(nn$(map_), errorf((log(lambda_)-log(p)+Sqr(sigma_)/2+mu)/sigma_))

$FUNCTION TechGeneral({p},{lambda_},{sigma_},{mu},{map_}): 
  $IF %techtype% == 'logit': TechLogit( ({p}), ({lambda_}), ({sigma_}), ({mu}), ({map_})) $ENDIF
  $IF %techtype% == 'normal': TechNormal( ({p}),({lambda_}),({sigma_}),({mu}),({map_})) $ENDIF
  $IF %techtype% == 'lognormal': TechLogNormal(({p}),({lambda_}),({sigma_}),({mu}),({map_})) $ENDIF
$ENDFUNCTION

$FUNCTION SolveEmptyNLP({name})
variable randomnameobj;  
randomnameobj.L = 0;

EQUATION E_EmptyNLPObj;
E_EmptyNLPObj..    randomnameobj  =E=  0;

Model M_SolveEmptyNLP /
E_EmptyNLPObj, {name}
/;
solve M_SolveEmptyNLP using NLP min randomnameobj;
$ENDFUNCTION