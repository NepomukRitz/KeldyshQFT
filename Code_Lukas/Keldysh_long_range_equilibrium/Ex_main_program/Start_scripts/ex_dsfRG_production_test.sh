L=5
Lu=0
N=30
Nff=1500
NfbP=1500
NfbX=1500
num_freq_pre=30000
NL_full=15

Vg=0.25
h=0.0
mu=-1.475
T=0.0
U0=0.65
U1=0.0
Xi=5.0

accp=0.0001
accx=0.0001
accs=0.0001
tol=0.000001
Lambda_ini=100000.0
Lambda_fin=0.000000002


PROGRAM=/p/home/jusers/weidinger1/juwels/bin/ex_dsfRG_production
OPTIONS="-L $L -Lu $Lu -N $N -Nff $Nff -NfbP $NfbP -NfbX $NfbX -num_freq_pre $num_freq_pre -NL_full $NL_full -Vg $Vg -h $h -mu $mu -T $T -U0 $U0 -U1 $U1 -Xi $Xi -accp $accp -accx $accx -accs $accs -tol $tol -Lambda_ini $Lambda_ini -Lambda_fin $Lambda_fin" 

DIR_CURRENT=$PWD
DIR_OUTPUT="/p/scratch/chmu26/hmu261/Ex_DATA/Production_test"

cd $DIR_OUTPUT
$PROGRAM $OPTIONS
cd $DIR_CURRENT
