from openpyxl import load_workbook
#Auteur: Stef van Breemen; Carlijn Fransen
def main():
    excelparse_writetofasta()

def excelparse_writetofasta():
    wb = load_workbook(filename="Course4_dataset_v03.xlsx", read_only=True)     #laad de excel file in
    worksheet = wb["groep13"]
    fasta_file = 'sequence.fasta'
    fasta=open(fasta_file, 'w')             #maakt een .fasta bestand aan, en opent hem in schrijf modus
    x=1
    y=2
    for row in worksheet.rows:
        header = row[0].value               #haalt de header uit de sequentie, replaced de @ met een > en haalt de sequentie op
        header = header.replace("@", ">")
        sequence = row[1].value
        fasta.write(header + '\n')          # schrijft alles weg naar een fasta bestand
        fasta.write(sequence+ '\n')
        fasta.write('\n')
        print('sequentie '+ str(x) +' toegevoegd')
        x+=2
        header = row[3].value               #haalt de header uit de sequentie, replaced de @ met een > en haalt de sequentie op
        header = header.replace("@", ">")
        sequence = row[4].value
        fasta.write(header+'\n')            # schrijft alles weg naar een fasta bestand
        fasta.write(sequence+"\n")
        fasta.write('\n')
        print('sequentie '+ str(y)+' toegevoegd')
        y+=2


main()