mindep=10
maxdep=10
model=sp6
#model=ak135
model_ray=raysp6prem
gcarclst=`gawk 'BEGIN {for(i=30;i<=30;i=i+1) printf "%5.2f ",i}'`
gawk -v gcarc="$gcarclst"  -v maxdep="$maxdep" -v mindep="$mindep" -v model=$model -v model_ray=$model_ray 'BEGIN {for (dep=maxdep;dep>=mindep;dep-=5) {
	print "time ./multitel3.pl -M"model"/"model_ray"/"dep" -O. -D "gcarc ;} }' | sh

./multisyn -M4.5/355/80/-70 -D1 -A33.5 -OPAS.z -Gsp6direct_10/30.00.grn.0 -Csp6core_10/30.00.grn.0

#########if you want to have green's function with both direct and core phases################################
for ((i=$mindep;i<=$maxdep;i=i+5))
 do
     if [ -d ./${model}_$i ];then
        rm -r ${model}_$i
     fi
     mkdir ${model}_$i
     saclst t3 f ${model}direct_$i/*.[258] | gawk -v model=$model '{split($1,aa,"_");print"cuterr fillz;cut "$2-30,$2+600 ;print"r "$1,model"core_"aa[2];print"w over";print"r "$1;print"addf "model"core_"aa[2];print"w "model"_"aa[2]}END{print"quit"}' | sac
     saclst t3 f ${model}direct_$i/*.[036147] | gawk -v model=$model '{split($1,aa,"_");print"cuterr fillz;cut "$2-30,$2+600 ;print"r "$1,model"core_"aa[2];print"w over";print"r "$1;print"addf "model"core_"aa[2];print"w "model"_"aa[2]}END{print"quit"}' | sac
#     rm -r ${model}direct_$i ${model}core_$i
done
