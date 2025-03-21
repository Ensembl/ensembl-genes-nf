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

import groovy.json.JsonSlurper

def make_publish_dir(publish_dir, project, name) {
    list = [publish_dir]
    if (project == "ensembl" && name != "") {
        list = list + [name]
    }
    return list.join("/")
}

def read_json(json_path) {
    slurp = new JsonSlurper()
    json_file = file(json_path)
    text = json_file.text

    not_a_lazy_val = slurp.parseText(text)
    return not_a_lazy_val
}

def index_metadata(json_path) {
    //Import meta JSON info from coredb
    metadata = read_json(json_path)

    // define core_metadata + export for pipeline meta propagation
    return [
        insdc_acc: metadata.get("assembly").get("accession"),
        taxonomy_id: metadata.get("species").get("taxonomy_id"),
        dbname: metadata.get("database").get("name"),
        production_name: metadata.get("species").get("production_name"),
        organism_name: metadata.get("species").get("scientific_name"),
        annotation_source: metadata.get("species").get("annotation_source"),
    ]
}