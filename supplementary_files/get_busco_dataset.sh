user=$1
host=$2
port=$3
dbname=$4
file=$ENSCODE/ensembl-genes-nf/supplementary_files/busco_dataset.txt
#mysql-ens-sta-5,4684,ensro,canis_lupus_gca905319855v2_core_106_1,mysql-ens-sta-5,4684,ensro,canis_lupus_gca905319855v2_core_106_1
#match the first correspondence between the list of metavalues and the list of busco datasets
mapfile result < <(mysql -N -u $user  -h $host -P $port -D $db -se'select meta_value from meta where meta_key="species.classification" order by meta_id')
tentative=0
for i in "${result[@]}";
do      
	species=$(echo "$i" | tr '[:upper:]' '[:lower:]')	
	if grep $species $file;then
            echo $i
            #echo $db	$(grep $species $file) >> test.txt
            break;
    fi
done
