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

@Grab('org.codehaus.groovy:groovy-all:2.2.2')
//@Grab(group='org.xerial', module='sqlite-jdbc', version='3.36.0.3')
import groovy.sql.Sql
import groovy.json.JsonOutput

//include { fromQuery } from 'plugin/nf-sqldb'

def make_publish_dir(publish_dir, project, name) {
    list = [publish_dir]
    if (project == "ensembl" && name != "") {
        list = list + [name]
    }
    return list.join("/")
}

def buildMetadata(gca, taxon_id) {
    def db_meta = ["gca": gca, "taxon_id": taxon_id]
    return db_meta
}

def getMetaValue(String dbname, String metaKey) {
    def sql
    println(dbname)
    def driver = 'com.mysql.cj.jdbc.Driver'
    // 'mysql-connector-j-8.0.31' //'org.hsqldb.jdbc.JDBCDriver' //'mysql-connector-j-8.0.31' //'com.mysql.jdbc.Driver' // 'com.mysql.cj.jdbc.Driver'
    def jdbcUrl = "jdbc:mysql://${params.host}:${params.port}/${dbname}"
    sql = Sql.newInstance(jdbcUrl, params.user_w,params.password,driver)
    def query = "SELECT meta_value FROM meta WHERE meta_key = ?"
    def result = sql.rows(query, [metaKey])

    //result= channel.sql.fromQuery(query, db: 'core_db')
    return result
}

/*
def generateMetadataJson(sqlFilePath, jsonName, outputDir) {
    def metadata = []
    println(sqlFilePath.toString())
    new File(sqlFilePath).eachLine { line ->
        // Skip lines that do not contain INSERT INTO statements
        if (line.startsWith("INSERT INTO meta")) {
            // Extract values from the INSERT INTO statement
            def matcher = line =~ /\((\d+), '([^']*)', '([^']*)'\);/
            if (matcher.matches()) {
                def speciesId = matcher[0][1].toInteger()
                def metaKey = matcher[0][2]
                def metaValue = matcher[0][3]

                // Construct metadata entry
                def entry = [species_id: speciesId, meta_key: metaKey, meta_value: metaValue]
                metadata << entry
            }
        }
    }

    // Construct JSON object
    def json = [meta: metadata]
    // Write JSON object to file
    def outputFile = new File(outputDir, jsonName)
    outputFile.withWriter { writer ->
        writer << groovy.json.JsonOutput.toJson(json)
    }
    println(outputFile)
    return outputFile.absolutePath
}
*/
