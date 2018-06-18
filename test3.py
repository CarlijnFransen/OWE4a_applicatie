# Author: Joshua Koopmans
# Date: 29-05-2018
# Descr.: This program fills the external database with results yielded from parsed XML files.

from flask import Flask, render_template, request, make_response
import mysql.connector
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez, SeqIO
import re

Entrez.email = "JJ.Koopmans@student.han.nl"

app = Flask(__name__)

@app.route("/")
# Input: Default values for variables
# Output: Renders html templates
# Descr.: Locates key word, parameter in URL. Connects to database, and outputs result in html table. Else return all results.
def index(search_term="", emptyList=[]):
    try:
        if "key_word" in request.args:
            if "parameter" in request.args:
                parameter = request.args.get("parameter")
                order = "DESC"
            else:
                parameter= 1
                order = "ASC"
            search_term = request.args.get("key_word")
            result_list = connect(search_term, parameter, order)
            # zoekwoord goes into module to search database and returns results
            return render_template("results.html", zoekwoord=search_term, List=result_list) #results here (maybe in table))


        #When page loads/no input: Asks for input
        else:
            return make_response(render_template("index.html", zoekwoord=search_term, List=emptyList))
    except:
        print("Error in function 'index' (route = '/') ")
        raise

# Input: -
# Output: Renders html templates
# Descr.: Calls functions
def connect(search_input, parameter, order):
    result_list = []
    host = 'ninjajoshuadx.mysql.pythonanywhere-services.com'
    user = 'ninjajoshuadx'
    passwd = 'blaat1234'
    db = 'ninjajoshuadx$project_blok_4'

    connection = mysql.connector.connect(host=host, user=user, password=passwd, database=db)
    cursor = connection.cursor()

    search = "%" + search_input + "%"

    query = """
        select query_seq_id, ACC_ID, PROTEIN_NAME, PROTEIN_FUNCTION, SPECIES, GENUS_NAME, FAMILY_NAME, e_value, IDENTITY_PERC, SIMILARITY_PERC, query_cover from blast_results
        join p_function on (blast_results.p_function_id =p_function.p_function_id)
        join genus on (blast_results.genus_id = genus.genus_id)
        join family on (genus.family_id = family.family_id)
        natural join query_seq
        where lower(protein_name) like lower('%s')
        or lower(species) like lower('%s')
        or lower(protein_function) like lower('%s')
        or lower(family_name) like lower('%s')
        or lower(genus_name) like lower('%s')
        or lower(acc_id) like lower('%s')
        order by %s %s
        """

    cursor.execute(query % (search, search, search, search, search, search, parameter, order))
    for row in cursor:
        result_list.append(tuple(row))
    return result_list

@app.route("/blast/")
# Input: Default values for variables
# Output: Renders html templates
# Descr.: If a sequence is located in the URL, check if sequence is not empty and if DNA. Else return blast home screen.
def blast(seq="", List=[""]):
    if "sequence" in request.args:
        seq = request.args.get("sequence")

        if seq != "":
            if validate_DNA(seq) == True:
                blast_obj = blast_seq(seq)
                result_list = read_xml(blast_obj)
                return render_template("results.html", zoekwoord=seq, List=result_list)
            else:
                #seq = "Not a DNA-sequence!"
                return render_template("blast.html")
        else:
            #seq = "No input!"
            return render_template("blast.html")

        return render_template("blast_results.html", List=[["fail"]])
    else:

        return render_template("blast.html")

# Input: The inputted sequence
# Output: Returns a boolean indicating if the input is DNA or not
# Descr.: Checks the inputted sequence using regex.
def validate_DNA(sequence):
    regex = "[ATGCN]*"
    check = re.search(regex, sequence.upper())
    if (check.end()-check.start()) == len(sequence):
        return True
    else:
        return False

# Input: The inputted sequence
# Output: Returns a blast_object in .xml file format
# Descr.: Blasts inputted DNA-sequence using qblast.
def blast_seq(seq):
    blast_object = NCBIWWW.qblast("blastx", "nr", str(seq), hitlist_size=10, expect=1e-5)
    return blast_object

# Input: blast_object from qblast
# Output: -
# Descr.: Reads, parses .xml file yielded by qblast. Variables are used to fill a html table.
def read_xml(blast_object):
    Entrez.email = 'JJ.Koopmans@student.han.nl'
    result_list = []

    E_VALUE_THRESH = 1
    blast_record = NCBIXML.read(blast_object)
    alignment_id = 0
    for alignment in blast_record.alignments:
        alignment_id += 1

        # Isolates accession code from results
        acc_id = str(alignment.title.split("|")[3])

        # Uses accession code to search in protein database
        handle = Entrez.efetch(db="protein", id=acc_id, rettype="gb", retmode="txt")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        # Below variables for filling database
        organism = record.annotations["organism"]
        family = record.annotations["taxonomy"][(len(record.annotations["taxonomy"]) - 2)]
        genus = record.annotations["taxonomy"][(len(record.annotations["taxonomy"]) - 1)]

        for hsp in alignment.hsps:
            if hsp.expect <= E_VALUE_THRESH:

                # Below more variables for filling database
                e_value = '{:.2e}'.format(hsp.expect)
                identity_perc = format(hsp.identities / hsp.align_length, '.2%')
                positive_perc = format(hsp.positives / hsp.align_length, '.2%')
            listt=list((acc_id, organism, genus, family, e_value, identity_perc, positive_perc))
            result_list.append(listt)
        return result_list


@app.route("/help/")
# Input: -
# Output: Renders help page template
# Descr.: Renders help page template
def help():
    return render_template("help.html")


if __name__ == '__main__':
    app.run()

