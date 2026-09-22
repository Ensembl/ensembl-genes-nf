process CHECK_DRY_RUN_OUTPUT {
    tag "$db_name"

    input:
        tuple val(db_name),
              path(dry_run_sql_out)

    script:
	"""
	awk '
		NF==3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ {
		  if ($2 != $3) {
		  	print "MISMATCH:", $1, $2, "!=", $3; bad=1
		  }
		}
		END { if (!bad) print "All staged_count/matched_count pairs match" }
		' ${dry_run_sql_out}
	"""
}
