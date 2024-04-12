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

// If the project is ensembl, append the name to the publish_dir
// Otherwise just use the publish_dir as is
def make_publish_dir(publish_dir, project, name) {
    list = [publish_dir]
    if (project == "ensembl" && name != "") {
        list = list + [name]
    }
    return list.join("/")
}

















import groovy.sql.Sql
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter

def checkTaxonomy(String jdbcUrl, String username, String password, String taxonId) {
    def sql = Sql.newInstance(jdbcUrl, username, password)
    
    try {
        def query = "SELECT * FROM meta  WHERE taxon_id = '${taxonId}'"
        def result = sql.rows(query)
        return result.size() > 0
    } catch (Exception ex) {
        ex.printStackTrace()}
    finally {
        sql.close()
    }
}

def getLastCheckedDate(String jdbcUrl, String username, String password, String taxonId) {
    def sql = Sql.newInstance(jdbcUrl, username, password)
    def lastCheckedDate = null

    try {
        def query = "SELECT last_check FROM meta WHERE taxon_id = '${taxonId}'"
        def result = sql.rows(query)

        if (result.size() > 0) {
            // Assuming 'last_check' is a date-like column
            // Adjust the date format pattern based on the actual format in your database
            def dateFormat = new SimpleDateFormat("yyyy-MM-dd") // Adjust the format if needed
            lastCheckedDate = dateFormat.parse(result[0].last_check)
        }
    } catch (Exception ex) {
        ex.printStackTrace()    
    } finally {
        sql.close()
    }

    return lastCheckedDate
}

def insertMetaRecord(String jdbcUrl, String username, String password, String taxonId) {
    def sql = Sql.newInstance(jdbcUrl, username, password)

    try {
        // Get the current date and time
        def currentDate = LocalDateTime.now()
        def dateFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd")
        def formattedDate = currentDate.format(dateFormatter)

        // Execute the SQL INSERT statement
        def insertQuery = "INSERT INTO meta (taxon_id, last_checked_date) VALUES ('${taxonId}', '${formattedDate}')"
        sql.executeUpdate(insertQuery, 'meta_id')
    } catch (Exception ex) {
        ex.printStackTrace()
    } finally {
        sql.close()
    }

}
def updateLastCheckedDate(String jdbcUrl, String username, String password, String taxonId) {
    def sql = Sql.newInstance(jdbcUrl, username, password)

    try {
        // Get the current date and time
        def currentDate = LocalDateTime.now()
        def dateFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd")
        def formattedDate = currentDate.format(dateFormatter)

        // Execute the SQL UPDATE statement
        def updateQuery = "UPDATE meta SET last_checked_date = '${formattedDate}' WHERE taxon_id = '${taxonId}'"
        sql.executeUpdate(updateQuery)
    } catch (Exception ex) {
        ex.printStackTrace()
    }finally {
        sql.close()
    }

}
def build_ncbi_path(gca, assembly_name) {
    final gca_splitted = gca.replaceAll("_","").tokenize(".")[0].split("(?<=\\G.{3})").join('/')
    return  'https://ftp.ncbi.nlm.nih.gov/genomes/all'  + '/' + gca_splitted + '/' + "$gca" +'_' + assembly_name.replaceAll(" ","_") + '/' + "$gca" + '_' + assembly_name.replaceAll(" ","_") + '_genomic.fna.gz'
}

def getPairedFastqsURL(String jdbcUrl, String username, String password, String run_accession) {
    def sql = Sql.newInstance(jdbcUrl, username, password)
    try {
        def query = "SELECT url FROM file INNER JOIN run ON run_id WHERE run_accession = '${run_accession}'"
        def result = sql.rows(query)
    } catch (Exception ex) {
    ex.printStackTrace()    
    } finally {
        sql.close()
    }

    return result
}

def checkFastqc(String jdbcUrl, String username, String password, String run_accession) {
    def sql = Sql.newInstance(jdbcUrl, username, password)
    def query = """ SELECT basic_statistics, per_base_sequence_quality, per_sequence_quality_scores, \
        per_base_sequence_content 
        FROM data_files df 
        INNER JOIN run r on df.run_id =r.run_id 
        WHERE r.run_id= '${run_accession}'
        """
    def qc_status = null 

    try {
        def result = sql.rows(query)
        // Process the results
        results.each { row ->
        def basicStatistics = row.basic_statistics
        def perBaseSequenceQuality = row.per_base_sequence_quality
        def perSequenceQualityScores = row.per_sequence_quality_scores
        def perBaseSequenceContent = row.per_base_sequence_content
        if (basicStatistics=='PASS' && perBaseSequenceQuality='PASS' &&
            perSequenceQualityScores='PASS' && perBaseSequenceContent='PASS') {
            // Execute the SQL UPDATE statement
            def updateQuery = "UPDATE RUN set qc_status = 'QC_PASS' WHERE run_id= '${run_accession}'"
            sql.executeUpdate(updateQuery)
            qc_status = 'QC_PASS'
            }
        else {
            // Execute the SQL UPDATE statement
            def updateQuery = "UPDATE RUN set qc_status = 'QC_FAIL' WHERE run_id= '${run_accession}'"
            sql.executeUpdate(updateQuery)
            qc_status = 'QC_FAIL'
            }
    }
    } catch (Exception ex) {
        ex.printStackTrace()}
    finally {
        sql.close()
    }

    return qc_status
}

def checkOverrepresentedSequences(String jdbcUrl, String username, String password, String run_accession) {
    def sql = Sql.newInstance(jdbcUrl, username, password)
    def query = """ SELECT overrepresented_sequences 
        FROM data_files df 
        INNER JOIN run r on df.run_id =r.run_id 
        WHERE r.run_id= '${run_accession}'
        """
    def overrepresented_sequences = null 

    try {
        def result = sql.rows(query)
        // Process the results
        results.each { row ->
        def OverrepresentedSequences = row.overrepresented_sequences
        
        if (OverrepresentedSequences=='WARN' OR OverrepresentedSequences=='FAIL') {
            overrepresented_sequences = True
            }
        else {
            overrepresented_sequences = False
            }
    }
    } catch (Exception ex) {
        ex.printStackTrace()}
    finally {
        sql.close()
    }
    return overrepresented_sequences
}
def concatString(string1, string2, string3){
    return string1 + '_'+string2 + '_'+string3
}

def calculateIndexBases(genomeFile) {
    def indexBases = Math.min(14, Math.floor((Math.log(genomeFile, 2) / 2) - 1))
    return indexBases
}

def getRunId(String jdbcUrl, String username, String password, String run_accession, String gca, String percentage_mapped) {
    def sql = Sql.newInstance(jdbcUrl, username, password)
    def run_id = null
    try {
        def query = "SELECT run_id  FROM run WHERE run_accession = '${run_accession}'"
        run_id = sql.rows(query)
        return run_id
    } catch (Exception ex) {
        ex.printStackTrace()
    } finally {
        sql.close()
    }
    
}

def updateFastqcStatus(String jdbcUrl, String username, String password, String run_accession) {
    def sql = Sql.newInstance(jdbcUrl, username, password)

    try {
        // Execute the SQL UPDATE statement
        def updateQuery = "UPDATE run SET qc_status = 'ALIGNED' WHERE run_accession = '${run_accession}'"
        sql.executeUpdate(updateQuery)
    } catch (Exception ex) {
        ex.printStackTrace()
    }finally {
        sql.close()
    }

}