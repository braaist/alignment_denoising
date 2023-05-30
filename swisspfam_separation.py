import pandas as pd
import mysql.connector as mariadb

mariadb_connection = mariadb.connect(user="ruslan_gumerov",
                                     password="",
                                     database="swisspfam",
                                     unix_socket="/run/mysqld/mysqld.sock")

cursor = mariadb_connection.cursor()

mnemonics = pd.read_csv("/home/ruslan_gumerov/noisy/swisspfam_mnemonic_taxid_order_df.csv", index_col = 0)

for tax in mnemonics["swisspfam_mnemonic"]:
    cursor.execute(f"CREATE TABLE {tax} (organism_code VARCHAR(255),"
                    "protein_code VARCHAR(255), uniprot_code VARCHAR(255),"
                    "protein_length INT, domain_code VARCHAR(255),"
                    "domain_start INT, domain_end INT, domain_length INT);")

    cursor.execute(f"INSERT INTO {tax} SELECT * FROM swisspfam WHERE swisspfam.organism_code = '{tax}';")
    mariadb_connection.commit()
    print(tax)

cursor.close()
