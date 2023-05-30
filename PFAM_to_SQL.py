import csv
import itertools
import pandas as pd
import mysql.connector as mariadb

if __name__ == "__main__":
    mariadb_connection = mariadb.connect(user="ruslan_gumerov",
                                         password="",
                                         database="swisspfam",
                                         unix_socket="/run/mysqld/mysqld.sock")

    cursor = mariadb_connection.cursor()

    # Prepare the insert statement
    add_data = ("INSERT INTO swisspfam "
                "(organism_code, protein_code, uniprot_code, protein_length, domain_code, domain_start, domain_end) "
                "VALUES (%s, %s, %s, %s, %s, %s, %s)")

    # Open the file and read it in chunks
    with open('/data/ruslan_gumerov/swisspfam', 'r') as file:
        csv_reader = csv.reader(file, delimiter='\t')
        for chunk in iter(lambda: list(itertools.islice(csv_reader, 1000)), []):
            insert_list = []
            for line in chunk:
                if len(line) != 0:
                    line = line[0]
                    if line.startswith(">"):
                        organism_code = line.split(" ")[0].split("_")[1]
                        protein_code = line.split(" ")[0].split("_")[0][1:]
                        uniprot_code = line.split(" ")[-3]
                        protein_length = line.split(" ")[-2]
                    else:
                        if len(line.split(" ")) != 1:
                            domain_code = [x for x in line.split(" ") if x.startswith("PF")][0]
                            domain_start = int(line.split(" ")[-1].split("-")[0])
                            domain_end = int(line.split(" ")[-1].split("-")[1])
                            insert_list.append((organism_code, protein_code, uniprot_code, protein_length, domain_code, domain_start, domain_end))

            cursor.executemany(add_data, insert_list)
            mariadb_connection.commit()

    cursor = mariadb_connection.cursor()

    mnemonics = pd.read_csv("/home/ruslan_gumerov/noisy/mnemonic_taxid_df.csv", index_col = 0)

    for tax in mnemonics["swisspfam_mnemonic"]:
        cursor.execute(f"CREATE TABLE {tax} (organism_code VARCHAR(255),"
                        "protein_code VARCHAR(255), uniprot_code VARCHAR(255),"
                        "protein_length INT, domain_code VARCHAR(255),"
                        "domain_start INT, domain_end INT, domain_length INT);")

        cursor.execute(f"INSERT INTO {tax} SELECT * FROM swisspfam WHERE swisspfam.organism_code = '{tax}';")
        mariadb_connection.commit()
        print(tax)

    cursor.close()
    mariadb_connection.close()