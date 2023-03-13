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


def build_ncbi_path(gca, assembly_name) {
  final gca_splitted = gca.replaceAll("_","").tokenize(".")[0].split("(?<=\\G.{3})").join('/')
  return  'https://ftp.ncbi.nlm.nih.gov/genomes/all'  + '/' + gca_splitted + '/' + "$gca" +'_' + assembly_name.replaceAll(" ","_") + '/' + "$gca" + '_' + assembly_name.replaceAll(" ","_") + '_genomic.fna.gz'
}


def concatString(string1, string2, string3){
 return string1 + '_'+string2 + '_'+string3
}

def getDataset(busco_lineage){
 String fileContents = new File(busco_lineage).text
 return fileContents
}
