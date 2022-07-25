user=ensro

host=mysql-ens-sta-5

port=4684
db=$1
file=/nfs/production/flicek/ensembl/genebuild/ftricomi/nextflow/busco_dataset.txt
#mysql-ens-sta-5,4684,ensro,canis_lupus_gca905319855v2_core_106_1,mysql-ens-sta-5,4684,ensro,canis_lupus_gca905319855v2_core_106_1
mapfile result < <(mysql -N -u $user  -h $host -P $port -D $db -se'select meta_value from meta where meta_key="species.classification" order by meta_id')
#echo ${result[2]}
tentative=0
for i in "${result[@]}";
do      

	species=$(echo "$i" | tr '[:upper:]' '[:lower:]')	
	if grep $species $file;then
	
	  #if  $tentative -ne 2;then 
      	    
   #a=$(grep $i $file)
   #echo grep "$i" $file
            echo $i
	    #$tentative=$tentative+1
            echo $db	$(grep $species $file) >> test.txt
            break;
         #fi
    fi
        

done


#string='My long string'
#if [[ $string == *"My long"* ]]; then
#  echo "It's there!"
#fi
