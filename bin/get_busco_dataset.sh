user=$1
host=$2
port=$3
dbname=$4
file=$5
#match the first correspondence between the list of metavalues and the list of busco datasets
mapfile result < <(mysql -u $user -h $host -P $port -D $dbname -NB -e "SELECT meta_value FROM meta where meta_key='species.classification' ORDER BY meta_id")

for i in "${result[@]}";
do      
	species=$(echo "$i" | tr '[:upper:]' '[:lower:]')	
	FILE=($(grep -e "^${species}_" $file))
	if [[ $? -eq 0 ]];then
            echo ${FILE[0]} | cut -d '.' -f1
            break;
    fi
done
