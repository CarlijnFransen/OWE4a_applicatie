# Author: Joshua Koopmans
# Date: 29-05-2018
# Descr.: This program fills the external database with results yielded from parsed XML files.

from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio import Entrez
import mysql.connector

Entrez.email = 'JJ.Koopmans@student.han.nl'

# Input: -
# Output: -
# Descr.: Calls functions
def main():
    connection, cursor = database_connect()
    #roll_back(cursor)
    xml_find(cursor)
    close(connection, cursor)

# Input: An open database connection
# Output:
# Descr.: Opens 200 .xml files, parses the files, the yielded variables are passed as arguments in the fill function.
def xml_find(cursor):
    E_VALUE_THRESH = 1
    alignment_id = 0
    protein_id = 99
    family_id = 399
    genus_id = 699
    for file_nr in range(200):
        filename = 'result' + str(file_nr+1) + '.xml'
        file = open(filename, 'r')
        blast_record = NCBIXML.read(file)

        print("results"+ str(file_nr+1)+ " / results"+str(200))
        for alignment in blast_record.alignments:
            alignment_id += 1
            protein_id += 1
            family_id += 1
            genus_id += 1
            print("Inserting into database row: %s" % (alignment_id))

            # Isolates accession code from results
            acc_id = str(alignment.title.split("|")[3])

            # Isolates protein name from results
            pre_pre_protein_name = alignment.title.split(";")
            pre_name = pre_pre_protein_name[0].replace(">gi", "").split("|")[4]
            protein_name = pre_name.replace("'", "")

            # Uses accession code to search in protein database
            handle = Entrez.efetch(db="protein", id=acc_id, rettype = "gb", retmode = "txt")
            record = SeqIO.read(handle, "genbank")
            handle.close()

            # Below variables for filling database
            accession_id = record.annotations["accessions"][0]
            organism = record.annotations["organism"]
            try:
                family = record.annotations["taxonomy"][(len(record.annotations["taxonomy"])-2)]
            except IndexError:
                family = "N/A"
            try:
                genus = record.annotations["taxonomy"][(len(record.annotations["taxonomy"])-1)]
            except:
                genus= "N/A"
            try:
                function_list = record.annotations["comment"].replace("\n", "").split(".")

                for item in function_list:
                    if "FUNCTION" in item.upper():
                        protein_function = item[11:]
                    else:
                        protein_function = "N/A"
            except:
                protein_function = "N/A"

            for hsp in alignment.hsps:
                if hsp.expect <= E_VALUE_THRESH:

                    # Below more variables for filling database
                    score = hsp.score
                    query_cover = hsp.align_length
                    e_value = '{:.2e}'.format(hsp.expect)
                    identity_perc = format((hsp.identities/hsp.align_length)*100 , '.2f')
                    positive_perc = format((hsp.positives/hsp.align_length)*100 , '.2f')
                    similarity_perc = float(positive_perc) - float(identity_perc)
            fill(cursor, alignment_id, accession_id, e_value, identity_perc, positive_perc, similarity_perc, query_cover, organism, protein_name, protein_id, family_id, family, genus_id, genus, protein_function, score, file_nr+1)
    print("-"*30)
    print("Total rows inserted: " + str(alignment_id))

# Input: -
# Output: Returns a database connection
# Descr.: Connects to an external database using a host, user, password
def database_connect():

    host = 'ninjajoshuadx.mysql.pythonanywhere-services.com'
    user = 'ninjajoshuadx'
    passwd = 'blaat1234'
    db = 'ninjajoshuadx$project_blok_4'

    connection = mysql.connector.connect(host=host, user=user, password=passwd, database=db)
    cursor = connection.cursor()

    return connection, cursor

# Input: An open database connection, variables
# Output: -
# Descr.: Inserts variables into "family", "genus", "p_function", and "blast_results" tables in database using SQL query, commits
def fill(cursor, alignment_id, accession_id, e_value, identity_perc, positive_perc, similarity_perc, query_cover, organism, protein_name, protein_id, family_id, family_name, genus_id, genus_name, protein_function, bit_score, file_nr):
    try:
        cursor.execute(
            "insert into family(family_id, family_name) values(%s, '%s')"
            % (family_id, family_name))
        cursor.execute("commit")
    except:
        print("ERROR in: family query.")

    try:
        cursor.execute(
            "insert into genus(genus_id, genus_name, family_id) values(%s, '%s', %s)"
            % (genus_id, genus_name, family_id))
        cursor.execute("commit")
    except:
        print("ERROR in: genus query.")
    try:
        cursor.execute(
            "insert into p_function(p_function_id, protein_function) values(%s, '%s')"
            % (protein_id, protein_function))
        cursor.execute("commit")
    except:
        print("ERROR in p_function query.")
    try:
        query="insert into blast_results(result_id, acc_id, bit_score, e_value, identity_perc, similarity_perc, query_cover, protein_name, species, genus_id, p_function_id, query_seq_id) values(%s, '%s', '%s', '%s',%s, %s, %s, '%s', '%s', %s, %s, %s)"
        filled_query = query % (alignment_id, accession_id, bit_score, e_value, identity_perc,similarity_perc, query_cover, protein_name, organism, genus_id, protein_id, file_nr)
        cursor.execute(filled_query)
        cursor.execute("commit")
    except:
        print("ERROR in: " + filled_query)

# Input: An open database connection
# Output: -
# Descr.: Closes database connection
def close(connection, cursor):
    connection.commit()
    cursor.close()
    connection.close()

# Input: An open database connection
# Output: -
# Descr.: Used to revert any commits executed in the database.
def roll_back(cursor):
    cursor.execute("rollback")
main()
