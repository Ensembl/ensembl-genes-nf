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

process REPEATMODELER {
    label 'repeatmodeler' 
    tag 'repeatmodeler'
    storeDir "${params.cacheDir}/repeatmodeler/"

    input:
    val path(rm_database)    

    output:
    val path(repeatmodeler_output)    

    script:
    """
    RepeatModeler -engine ncbi -numAddlRounds 1 -pa 10 -database ${rm_database}/repeatmodeler_db
    """
}