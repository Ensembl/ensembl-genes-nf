#!/usr/bin/env nextflow
/*
See the NOTICE file distributed with this work for additional information
regarding copyright ownership.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

process CREATE_SLICES {
        input:
        file genome_file from params.genome_file

        output:
        file 'slice_ids.txt' into slice_ids_ch

        script:
        """
        python -c "
        from ensembl.tools.anno.utils._utils import get_seq_region_length, get_slice_id
        genome_file = '$genome_file'
        seq_region_to_length = get_seq_region_length(genome_file, 5000)
        slice_ids_per_region = get_slice_id(seq_region_to_length, slice_size=1000000, overlap=0, min_length=5000)
        with open('slice_ids.txt', 'w') as f:
            for slice_id in slice_ids_per_region:
                f.write(f'{slice_id[0]}\\t{slice_id[1]}\\t{slice_id[2]}\\n')
        "
        """
    }
