$PROBLEM    Jumbe Resistance Model

$INPUT      TIME DV ID DOSE MODEL
; ------------------------Load in first efficacy study only
$DATA       Olaparib_monotherapy_data.csv
IGNORE=@
;--------------------Tumour Parameters to be estimated-----
$SUBROUTINES ADVAN13 TOL=8

$MODEL      NCOMPARTMENTS=2
$OMEGA  0.0   FIX;     	omega for the extra effect of 100mg/kg
$OMEGA  0.0   FIX;     	omega for the 50mg/kg dose effect
$OMEGA  0.64  FIX;     	omega for the growth rate
$OMEGA  0.8  ;   	omega for the resistance fraction
$OMEGA  0.4    ;   	omega for the initial tumour volume
$OMEGA  0.5  ;    	omega on the PDX level for the resistance fraction
$THETA  1  ;      	THETA for the extra effect of 100mg/kg
$THETA  -1  ;       	THETA for 50mg/kg
$THETA  -3.2    FIX;    THETA for growth rate
$THETA  -2       ;    	THETA for olaparib resistance fraction
$THETA -1.4      ;   	THETA for baseline tumour volume
$SIGMA  0.3       ;     Additive Error

$LEVEL
MODEL=(6[4])

$PK
;main parameter list
Ks=EXP(THETA(3)+ETA(3))
Dtr= Ks + EXP(THETA(2)+(THETA(1)+ETA(1))*((DOSE-50)/50)+ETA(2))
Fr_0=EXP(THETA(4)+ETA(4)+ETA(6))/(1+EXP(THETA(4)+ETA(4)+ETA(6)))
INIT=EXP(THETA(5)+ETA(5))

;set up initial condition for tumour model

IF(Fr_0.GT.1) THEN
  Fr_0=1.0
ENDIF


A_0(1)=INIT*(1-Fr_0)
A_0(2)=INIT*Fr_0


;differential equations
$DES

Vs=A(1)
Vr=A(2)
V=A(1) + A(2)


;Jumbe model


;labelled differential equations

IF (V.LE.0) THEN
  V= 0.0
  DVsDT=0
  DVrDT=0
ELSE
  DVsDT=(Ks-Dtr)*Vs*(V)**(-1/3)
  DVrDT=Ks*Vr*(V)**(-1/3)
ENDIF

;differential equations for NONMEM
DADT(1)=DVsDT
DADT(2)=DVrDT

$ERROR (ONLY OBSERVATIONS)
Vol=A(1)+A(2)
IPRE=log(Vol)
Y=IPRE+ERR(1)
W=IPRE
IF (W.LT.1E-3) W=1E-3
IRES=DV-IPRE
IWRE=IRES


;output table


$ESTIMATION METHOD=CONDITIONAL MAXEVAL=9999 PRINT=5 NOHABORT
            FNLETA=0 SLOW NONINFETA = 1
$COVARIANCE SLOW MATRIX=S

$TABLE ID TIME MDV IPRE IWRE IRES CWRES NOPRINT ONEHEADER FILE=sdtab
$TABLE ID NOPRINT NOAPPEND ONEHEADER FILE=cotab
$TABLE ID DOSE Dtr Fr_0 Ks INIT ETA1 ETA2 ETA3 ETA4 ETA5 ETA6
     NOPRINT NOAPPEND ONEHEADER FILE=patab