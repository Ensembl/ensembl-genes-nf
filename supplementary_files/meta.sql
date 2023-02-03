select replace(group_concat(meta_value order by field(meta_key,'species.scientific_name','assembly.accession','species.annotation_source') SEPARATOR '/'), ' ','_') from meta where meta_key in ('species.scientific_name','assembly.accession','species.annotation_source');

