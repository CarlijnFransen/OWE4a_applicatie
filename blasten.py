#auteur; Carlijn Fransen
#datum: 16 mei
#dit programma leest een fasta file in en BLAST alle sequenties die in die fasta file staan, en slaat de resultaten op in een .xml file
'''---------------------------------------------------------------------------------------------------------------------------------------'''
#from Bio import SeqIO
from Bio.Blast import NCBIWWW
#input: fasta file
#output: xml files met BLAST resultaten
def main():
    bestand=read_file()
    headers, seqs = split(bestand)
    blast_and_save(headers, seqs)


def read_file():
    inputfasta = True
    while inputfasta == True:
        try:
            fastafile = input('Geef uw fasta file op: ')
            if '.fa' in fastafile or '.fasta' in fastafile:
                inputfasta= False
                bestand = open(fastafile)                            # leest fasta file in
            else:
                print('Dit is geen .fasta file, geef een .fasta file op.')
                inputfasta= True
        except FileNotFoundError:
            inputfasta= True
            print('Het bestand is niet gevonden, voer een geldige bestandsnaam in.')
    return bestand

def split(bestand):

    headers = []
    seqs = []
    seq = ""

    for line in bestand:                # zorgt dat de headers van de sequentie worden gescheiden
        line=line.strip()
        if ">" in line:
            if seq != "":
                seqs.append(seq)
                seq = ""
            headers.append(line)
        else:
            seq += line.strip()
    seqs.append(seq)
    #print(headers)
    #print(seqs)
    #record = SeqIO.read('fasta.fasta', format='fasta')
    return(headers, seqs)

def blast_and_save(headers, seqs):
    print('Er zijn '+str(len(headers))+ ' sequenties gevonden om te BLASTen.')
    for x in range(len(headers)):

        sequentie=seqs[x]
        sequentie = sequentie.upper()
        y = x+1
        print('Bezig met BLASTen van sequentie '+str(y)+ ' dit kan een aantal minuten duren.')
        result_handle = NCBIWWW.qblast('blastx','nr',sequentie, hitlist_size=10, expect=1e-5)   # er wordt hier geblast, BLOUSUM62 is default,

        print('BLAST-resultaten aan het opslaan... een moment geduld.')
        file_name = 'result'+ str(y)+'.xml'
        with open(file_name, 'w') as out_handle:                                                # de BLAST resultaten worden hier opgeslagen
            out_handle.write(result_handle.read())
        result_handle.close()
        print('De BLAST-resultaten van sequentie '+str(y)+' zijn opgeslagen.')

main()