# Author: Stef van Breemen
# Date: (last edited: 14-06-2018)
# Descr.: This program fills the external database tables "orientation" and "query_seq" with values using an excel file.

from openpyxl import load_workbook
import mysql.connector

# Input: -
# Output: -
# Descr.: Calls functions
def main():
    connection, cursor = connect()
    fill_orien(cursor)
    excelparse(cursor)
    close(connection, cursor)

# Input: -
# Output: Returns a database connection
# Descr.: Connects to an external database using a host, user, password
def connect():
    host = 'ninjajoshuadx.mysql.pythonanywhere-services.com'
    user = 'ninjajoshuadx'
    passwd = 'blaat1234'
    db = 'ninjajoshuadx$project_blok_4'
    connection = mysql.connector.connect(host=host, user=user, password=passwd, database=db)
    cursor = connection.cursor()

    return connection, cursor

# Input: An open database connection
# Output: -
# Descr.: Parses the .xlsx file. Calls function to fill database.
def excelparse(cursor):
    wb = load_workbook(filename="Course4_dataset_v03.xlsx", read_only=True)
    worksheet = wb["groep13"]
    pk = 1

    for row in worksheet.rows:
        header = row[0].value
        header = header.replace("@", ">")
        sequence = row[1].value
        asciiscore = sum([ord(letter) - 32 for letter in row[2].value])
        fill(cursor, pk, header, sequence, asciiscore, 1)
        pk += 1
        header = row[3].value
        header = header.replace("@", ">")
        sequence = row[4].value
        asciiscore = sum([ord(letter) - 32 for letter in row[5].value])
        fill(cursor, pk, header, sequence, asciiscore, 2)
        pk += 1

# Input: An open database connection
# Output: -
# Descr.: Inserts variables into "orientation" table in database using SQL query, commits
def fill_orien(cursor):
    cursor.execute("insert into orientation(orientation_id, orientation_value) values (1, 'Forward')")
    cursor.execute("commit")
    cursor.execute("insert into orientation(orientation_id, orientation_value) values (2, 'Reverse')")
    cursor.execute("commit")

# Input: An open database connection, variables
# Output: -
# Descr.: Inserts variables into "query_seq" table in database using SQL query, commits
def fill(cursor, pk, header, sequence, asciiscore, orientation):
    cursor.execute(
        "insert into query_seq(query_seq_id, HEADER, sequence, Ascii_score, orientation_id) values(%s, '%s', '%s', %s, %s)"
        % (pk, header, sequence, asciiscore, orientation))
    cursor.execute("commit")

# Input: An open database connection
# Output: -
# Descr.: Closes database connection
def close(connection, cursor):
    connection.commit()
    cursor.close()
    connection.close()


main()



