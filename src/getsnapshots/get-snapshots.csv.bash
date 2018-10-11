#!/bin/bash

simulate=0   # don't write csv file

# define variables
variables='Min Cycle Frame Omega M66cAcB M66cBcG M66cGsD alpha alphap beta gamma phi tau E16cAcB E16rot E16orie L199rot K70cDcE K70cGcD K70cBcG I65cAcB I65rot Q109rot S111rot S62rot L119rot hbE16acy hbI65acy bE16oh I197rot Q213rot hbE215imi hbW239o64 hbW239o66 hbW239E16 hbE16W239 hbW240imi hbW240E215 hbK70E148 hbK70imi hbW304imi hbW239acy hbQ213acy hbK70W301 dW304Trp93 dW240Q42 M18cBcG M18cGsD hbW230Q109 dW301Q109 hbW304K70 o66orie bb66dihe hbW260phe hbS146phe hbR95imi hbW240n2 hbW230S111 hbW230o66 hbW230E16 hbW230o64 hbW230W239 hbW301W230 hbW301W304 hbW301S111 hbW301o66 hbW304E148 W301rmsd W304rmsd'

nvar=$(echo $variables | wc -w)
declare -A var
for i in $variables ; do
  var[$i]="NA"
done
echo "$nvar variables: ${!var[*]}"
# check completeness of variables
nassigned=$(echo ${var[*]} | wc -w)
if [ $nassigned -ne $nvar ]; then echo "number of assigned variables ($nassigned) is not $nvar" ; head -n 1 $csv; tail -n 1 $csv ; exit ; fi

# write first line of csv file with all necesary labels
if [ $simulate -eq 1 ]; then
  echo "just a simulation..."
  csv=/dev/null
else
#  csv=snapshots-s0.csv
#  csv=snapshots-old.csv
#  csv=snapshots-a.csv
  csv=snapshots-ref.csv
  echo "${!var[*]}" | sed 's/  */,/g' > $csv
fi

#Mins="a_qm12 a_s0-charges-fix-met105-lys-ile-glu-ser a_s0-charges-fix-met105-lys-ile-glu-ser-hb v38s0"
#Mins="a a_fix-lys-ile a_fix-lys-ile-glu a_fix-met105 a_fix-met105-glu a_fix-met105-lys a_fix-met105-lys-glu a_fix-met105-lys-ile a_fix-met105-lys-ile-glu a_fix-met105-lys-ile-glu-ser a_fix-lys-ile-glu-ser a_fix-met105-lys-ile-glu-ser-hb a_fix-lys-ile-glu-ser-hb"
#Mins="v36b v36c v38a v38b v38c v38d v38f v38g v38h v39d v38i v38j v38k v38l v38m v38n v38o v38p v38q v38r v38s v38t v38u v38v v38w v38x v38y"
Mins="v38s"

for Min in $Mins ; do
echo "$Min ..."

for ((Cycle=2;Cycle<=92;Cycle++)) ; do
for ((Frame=10;Frame<=200;Frame+=10)) ; do

# initialize variables (atom indices from ~/bin/get-other-dihedrals)
for subscript in $variables ; do
  var[$subscript]="NA"
done

# assign all variables from final excited-state optimization (cro-opt)
opt=cro-opt
i=${Cycle}-$Frame
log=snapshots-$Min/${opt}-${i}.log
sge=snapshots-$Min/${opt}-${i}.sge
mndoout=snapshots-$Min/${opt}-${i}-mndo.out
optc=snapshots-$Min/${opt}-${i}-opt.c

if ! [ -f $log ] && ! [ -f $sge ]; then
  continue
fi
if ! [ -f $log ]; then
  echo "    Warning: found $sge but not $log." ; continue
fi
if [ `grep -c 'Converged' $log` -eq 0 ]; then
  echo "    Warning: $log not converged." ; continue
fi
if ! [ -f $mndoout ]; then
  echo "$mndoout missing although $log is converged!" ; exit
elif [ $(grep -ce 'State  2' $mndoout ) -ne 1 ]; then
  echo "no excitation energy found in $mndoout" ; exit
fi
if ! [ -f $optc ]; then
  echo "$optc missing although $log signals convergence!" ; exit
fi

if [ $simulate -ne 0 ]; then
  echo "    $optc: Min=$Min Cycle=$Cycle Frame=$Frame"
  continue
else
  echo -n "."
fi

if [ $HOSTNAME = ganbo141 ]; then
  c2xyz $optc > /scratch/wanko/283950/tmp.xyz || exit
  xyz=/scratch/wanko/283950/tmp.xyz
else
  c2xyz $optc > /gscratch/wanko/tmp.xyz ; xyz=/gscratch/wanko/tmp.xyz
fi

var[Min]=$Min
var[Cycle]=$Cycle
var[Frame]=$Frame
var[Omega]=$(grep 'State  2' $mndoout | awk '{printf("%.3f",$6)}')
var[M66cAcB]=`getdihedral_fast 959 947 948 951 $xyz | awk '$1<0{printf("%.2f",360+$1)} $1>=0{printf("%.2f",$1)}'`
var[M66cBcG]=`getdihedral_fast 947 948 951 954 $xyz | awk '$1<0{printf("%.2f",360+$1)} $1>=0{printf("%.2f",$1)}'`
var[M66cGsD]=`getdihedral_fast 948 951 954 955 $xyz | awk '$1<0{printf("%.2f",360+$1)} $1>=0{printf("%.2f",$1)}'`
var[alpha]=`getdihedral_fast 946 947 959 961 $xyz | awk '{printf("%.2f",$1)}'`
var[alphap]=`getdihedral_fast 948 947 959 961 $xyz | awk '$1<0{printf("%.2f",360+$1)} $1>=0{printf("%.2f",$1)}'`
var[beta]=`getdihedral_fast 944 946 947 959 $xyz | awk '$1<0{printf("%.2f",360+$1)} $1>=0{printf("%.2f",$1)}'`
var[gamma]=`getdihedral_fast 945 944 946 947 $xyz | awk '{printf("%.2f",$1)}'`
var[phi]=`getdihedral_fast 973 972 970 964 $xyz | awk '$1<0{printf("%.2f",0+$1)} $1>=0{printf("%.2f",$1)}'`
var[tau]=`getdihedral_fast 972 970 964 960 $xyz | awk '$1<0{printf("%.2f",0+$1)} $1>=0{printf("%.2f",$1)}'`
var[E16cAcB]=`getdihedral_fast 193 195 197 200 $xyz | awk '$1<0{printf("%.2f",360+$1)} $1>=0{printf("%.2f",$1)}'` # N-CA-CB-CG
var[E16rot]=`getdihedral_fast 195 197 200 203 $xyz | awk '$1<0{printf("%.2f",360+$1)} $1>=0{printf("%.2f",$1)}'`  # CA-CB-CG-CD
var[E16orie]=`getdihedral_fast 197 200 203 205 $xyz | awk '$1<0{printf("%.2f",0+$1)} $1>=0{printf("%.2f",$1)}'`   # CB-CG-CD-OE2
var[L199rot]=`getdihedral_fast 3049 3051 3054 3055 $xyz | awk '$1<0{printf("%.2f",0+$1)} $1>=0{printf("%.2f",$1)}'` #CA-CB-CG-HG
var[K70cDcE]=`getdihedral_fast 1001 1004 1007 1010 $xyz | awk '$1<0{printf("%.2f",360+$1)} $1>=0{printf("%.2f",$1)}'`
var[K70cGcD]=`getdihedral_fast 998 1001 1004 1007 $xyz | awk '$1<0{printf("%.2f",360+$1)} $1>=0{printf("%.2f",$1)}'`
var[K70cBcG]=`getdihedral_fast 996 998 1001 1004 $xyz | awk '$1<0{printf("%.2f",360+$1)} $1>=0{printf("%.2f",$1)}'`
var[I65cAcB]=`getdihedral_fast 927 929 931 937 $xyz | awk '$1<0{printf("%.2f",0+$1)} $1>=0{printf("%.2f",$1)}'` #  N CA CB CG1
var[I65rot]=`getdihedral_fast 929 931 937 940 $xyz | awk '$1<0{printf("%.2f",360+$1)} $1>=0{printf("%.2f",$1)}'`  # CA CB CG1 CD
var[Q109rot]=`getdihedral_fast 1623 1625 1628 1631 $xyz | awk '$1<0{printf("%.2f",360+$1)} $1>=0{printf("%.2f",$1)}'` # CA CB CG CD
var[S111rot]=`getdihedral_fast 1650 1652 1654 1657 $xyz | awk '$1<0{printf("%.2f",360+$1)} $1>=0{printf("%.2f",$1)}'` #N CA CB OG
var[S62rot]=`getdihedral_fast 885 887 889 892 $xyz | awk '$1<0{printf("%.2f",360+$1)} $1>=0{printf("%.2f",$1)}'` #N CA CB OG
var[L119rot]=`getdihedral_fast 1763 1765 1771 1774 $xyz | awk '$1<0{printf("%.2f",0+$1)} $1>=0{printf("%.2f",$1)}'`
var[hbE16acy]=`getdistance_fast 206 945 $xyz | awk '{printf("%.3f",$1)}'`
var[hbI65acy]=`getdistance_fast 928 946 $xyz | awk '{printf("%.3f",$1)}'` # artifact of qm3 in W12: stong HB from I65_HN (MM) to ACY_N
var[hbE215imi]=`getdistance_fast 3317 960 $xyz | awk '{printf("%.3f",$1)}'` # PROT 215 HE2  PROT 66 N2 (2.08 in qm12 MD)
d1=`getdistance_fast 3491 926 $xyz` ; d2=`getdistance_fast 3492 926 $xyz`
var[hbW239o64]=`echo $d1 $d2 | awk '{if($1<$2) d=$1; else d=$2; printf("%.3f",d)}'` # W239_H-64_O
d1=`getdistance_fast 3491 969 $xyz` ; d2=`getdistance_fast 3492 969 $xyz`
var[hbW239o66]=`echo $d1 $d2 | awk '{if($1<$2) d=$1; else d=$2; printf("%.3f",d)}'` # W239_H-CRO_O
d1=`getdistance_fast 3491 204 $xyz` ; d2=`getdistance_fast 3492 204 $xyz`
var[hbW239E16]=`echo $d1 $d2 | awk '{if($1<$2) d=$1; else d=$2; printf("%.3f",d)}'` # W239_H - E16_OE1  (W239->E16->acy)
var[hbE16W239]=`getdistance_fast 206 3490 $xyz | awk '{printf("%.3f",$1)}'`              # E16_HE2 - W239_O  (E16->W239->acy)
d1=`getdistance_fast 3494 963 $xyz` ; d2=`getdistance_fast 3495 963 $xyz`
var[hbW240imi]=`echo $d1 $d2 | awk '{if($1<$2) d=$1; else d=$2; printf("%.3f",d)}'` # W240_H-CRO_O2 (imi)
d1=`getdistance_fast 3494 3315 $xyz` ; d2=`getdistance_fast 3495 3315 $xyz`
var[hbW240E215]=`echo $d1 $d2 | awk '{if($1<$2) d=$1; else d=$2; printf("%.3f",d)}'` # W240_H-E215_OE2
var[bE16oh]=`getdistance_fast 205 206 $xyz | awk '{printf("%.3f",$1)}'`   # the E16_O-H bond weakens in partial PT to ACY65
var[I197rot]=`getdihedral_fast 3008 3010 3016 3019 $xyz | awk '$1<0{printf("%.2f",360+$1)} $1>=0{printf("%.2f",$1)}'` # superpos in all MDs
var[Q213rot]=`getdihedral_fast 3266 3268 3270 3273 $xyz | awk '$1<0{printf("%.2f",0+$1)} $1>=0{printf("%.2f",$1)}'` # superpos in all MDs
d1=`getdistance_fast 1010 2229 $xyz` ; d2=`getdistance_fast 1010 2230 $xyz` # PROT 70 NZ  PROT 148 OE1`
var[hbK70E148]=`echo $d1 $d2 | awk '{if($1<$2) d=$1; else d=$2; printf("%.3f",d)}'` # K70_NZ-E148_OE*
var[hbK70imi]=`getdistance_fast 1010 963 $xyz | awk '{printf("%.3f",$1)}'` # PROT 70 NZ  PROT 66_O2
d1=`getdistance_fast 3686 963 $xyz` ; d2=`getdistance_fast 3687 963 $xyz`
var[hbW304imi]=`echo $d1 $d2 | awk '{if($1<$2) d=$1; else d=$2; printf("%.3f",d)}'` # W304_H - CRO_O2 (imi)
d1=`getdistance_fast 3491 945 $xyz` ; d2=`getdistance_fast 3492 945 $xyz` 
var[hbW239acy]=`echo $d1 $d2 | awk '{if($1<$2) d=$1; else d=$2; printf("%.3f",d)}'` # W239_H - I65_O  (E16->W239->acy)
d1=`getdistance_fast 3279 945 $xyz` ; d2=`getdistance_fast 3280 945 $xyz`
var[hbQ213acy]=`echo $d1 $d2 | awk '{if($1<$2) d=$1; else d=$2; printf("%.3f",d)}'` # Q213_HE2 - I65_O (acy)
var[hbK70W301]=`getdistance_fast 995 3676 $xyz | awk '{printf("%.3f",$1)}'` # PROT 70 HN - W301_O
var[dW304Trp93]=`getdistance_fast 1384 3685 $xyz | awk '{printf("%.3f",$1)}'` # PROT 93 NE1 - W304_O
var[dW240Q42]=`getdistance_fast 576 3493 $xyz | awk '{printf("%.3f",$1)}'` # PROT 42 HE22 - W240_O
var[M18cBcG]=`getdihedral_fast 229 231 234 237 $xyz | awk '$1<0{printf("%.2f",0+$1)} $1>=0{printf("%.2f",$1)}'` # CA CB CG SD
var[M18cGsD]=`getdihedral_fast 231 234 237 238 $xyz | awk '$1<0{printf("%.2f",0+$1)} $1>=0{printf("%.2f",$1)}'` # CB CG SD CE
var[hbW230Q109]=`getdistance_fast 1635 3463 $xyz | awk '{printf("%.3f",$1)}'` # Q109_HE22 - W230_OH2
var[dW301Q109]=`getdistance_fast 1635 3676 $xyz | awk '{printf("%.3f",$1)}'` # Q109_HE22 - W301_OH2
var[hbW304K70]=`getdistance_fast 995 3685 $xyz | awk '{printf("%.3f",$1)}'` # K70_HN - W304_OH2
var[o66orie]=`getdihedral_fast 969 968 965 961 $xyz | awk '$1<0{printf("%.2f",360+$1)} $1>=0{printf("%.2f",$1)}'`  # O C CA3 N3
var[bb66dihe]=`getdihedral_fast 968 965 961 959 $xyz | awk '$1<0{printf("%.2f",360+$1)} $1>=0{printf("%.2f",$1)}'`  # C CA3 N3 C1
d1=`getdistance_fast 3555 982 $xyz` ; d2=`getdistance_fast 3554 982 $xyz`
var[hbW260phe]=`echo $d1 $d2 | awk '{if($1<$2) d=$1; else d=$2; printf("%.3f",d)}'` # W260 - phe (OH)
var[hbS146phe]=`getdistance_fast 2204 982 $xyz | awk '{printf("%.2f",$1)}'`              # Ser146_HG1 - phe (OH)
d1=`getdistance_fast 1433 963 $xyz` ; d2=`getdistance_fast 1431 963 $xyz`
var[hbR95imi]=`echo $d1 $d2 | awk '{if($1<$2) d=$1; else d=$2; printf("%.3f",d)}'` # Arg95_HH21/HH12 - imi (O2)
d1=`getdistance_fast 3495 960 $xyz` ; d2=`getdistance_fast 3494 960 $xyz`
var[hbW240n2]=`echo $d1 $d2 | awk '{if($1<$2) d=$1; else d=$2; printf("%.3f",d)}'` # W240 - imi (N2)
d1=`getdistance_fast 3464 1657 $xyz` ; d2=`getdistance_fast 3465 1657 $xyz`
var[hbW230S111]=`echo $d1 $d2 | awk '{if($1<$2) d=$1; else d=$2; printf("%.3f",d)}'` # SOLV 230 H1/H2  PROT 111 OG
d1=`getdistance_fast 3464 969 $xyz` ; d2=`getdistance_fast 3465 969 $xyz`
var[hbW230o66]=`echo $d1 $d2 | awk '{if($1<$2) d=$1; else d=$2; printf("%.3f",d)}'` # SOLV 230 H1/H2  PROT  66 O
d1=`getdistance_fast 3464 204 $xyz` ; d2=`getdistance_fast 3465 204 $xyz`
var[hbW230E16]=`echo $d1 $d2 | awk '{if($1<$2) d=$1; else d=$2; printf("%.3f",d)}'` # SOLV 230 H1/H2  PROT  16 OE1
d1=`getdistance_fast 3464 926 $xyz` ; d2=`getdistance_fast 3465 926 $xyz`
var[hbW230o64]=`echo $d1 $d2 | awk '{if($1<$2) d=$1; else d=$2; printf("%.3f",d)}'` # SOLV 230 H1/H2  PROT  64 O
d1=`getdistance_fast 3464 3490 $xyz` ; d2=`getdistance_fast 3465 3490 $xyz`
var[hbW230W239]=`echo $d1 $d2 | awk '{if($1<$2) d=$1; else d=$2; printf("%.3f",d)}'` # SOLV 230 H1/H2  SOLV 239 OH2
d1=`getdistance_fast 3677 3463 $xyz` ; d2=`getdistance_fast 3678 3463 $xyz`
var[hbW301W230]=`echo $d1 $d2 | awk '{if($1<$2) d=$1; else d=$2; printf("%.3f",d)}'` # SOLV 301 H1/H2  SOLV 230 OH2
d1=`getdistance_fast 3677 3685 $xyz` ; d2=`getdistance_fast 3678 3685 $xyz`
var[hbW301W304]=`echo $d1 $d2 | awk '{if($1<$2) d=$1; else d=$2; printf("%.3f",d)}'` # SOLV 301 H1/H2  SOLV 304 OH2
d1=`getdistance_fast 3677 1657 $xyz` ; d2=`getdistance_fast 3678 1657 $xyz`
var[hbW301S111]=`echo $d1 $d2 | awk '{if($1<$2) d=$1; else d=$2; printf("%.3f",d)}'` # SOLV 301 H1/H2  PROT 111 OG
d1=`getdistance_fast 3677 969 $xyz` ; d2=`getdistance_fast 3678 969 $xyz`
var[hbW301o66]=`echo $d1 $d2 | awk '{if($1<$2) d=$1; else d=$2; printf("%.3f",d)}'` # SOLV 301 H1/H2  PROT  66 O
d1=`getdistance_fast 3686 2229 $xyz` ; d2=`getdistance_fast 3687 2229 $xyz`
var[hbW304E148]=`echo $d1 $d2 | awk '{if($1<$2) d=$1; else d=$2; printf("%.3f",d)}'` # SOLV 304 H1/H2  PROT 148 OE1
var[W301rmsd]=`awk 'NR-2==3676{print sqrt(($2-x)^2+($3-y)^2+($4-z)^2)}' x=-4.203 y=5.095 z=-2.508 $xyz`
var[W304rmsd]=`awk 'NR-2==3685{print sqrt(($2-x)^2+($3-y)^2+($4-z)^2)}' x=-4.539 y=3.184 z=-0.377 $xyz`


# Write variables to csv.
echo ${var[*]} | sed 's/  */,/g' >> $csv

# check completeness of variables
nassigned=$(echo ${var[*]} | wc -w)
if [ $nassigned -ne $nvar ]; then echo "number of assigned variables ($nassigned) is not $nvar" ; head -n 1 $csv; tail -n 1 $csv ; exit ; fi

done
done
echo "$Min DONE"
done

echo "variables $!{var[*]} written to $csv"
