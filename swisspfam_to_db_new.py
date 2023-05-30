import mysql.connector as mariadb
import csv
import itertools


if __name__ == "__main__":
    mariadb_connection = mariadb.connect(user="ruslan_gumerov",
                                         password="",
                                         database="swisspfam",
                                         unix_socket="/run/mysqld/mysqld.sock")

    cursor = mariadb_connection.cursor()

    # Prepare the insert statement
    add_data = ("INSERT INTO swisspfam "
                "(organism_code, protein_code, uniprot_code, protein_length, domain_code, domain_start, domain_end, domain_length) "
                "VALUES (%s, %s, %s, %s, %s, %s, %s, %s)")

    cursor.execute(f"CREATE TABLE swisspfam (organism_code VARCHAR(255),"
                        "protein_code VARCHAR(255), uniprot_code VARCHAR(255),"
                        "protein_length INT, domain_code VARCHAR(255),"
                        "domain_start INT, domain_end INT, domain_length INT);")

    mariadb_connection.commit()

    # Open the file and read it in chunks
    with open('/data/ruslan_gumerov/swisspfam', 'r') as file:
        csv_reader = csv.reader(file, delimiter='\t')
        chunk_num = 0
        for chunk in iter(lambda: list(itertools.islice(csv_reader, 10000)), []):
            print(chunk_num)
            chunk_num += 1
            insert_list = []
            broken_protein_list = []
            organism_code_list = set()
            for line in chunk:
                if len(line) != 0:
                    line = line[0]
                    if line.startswith(">"):
                        organism_code = line.split(" ")[0].split("_")[1]
                        protein_code = line.split(" ")[0].split("_")[0][1:]
                        uniprot_code = line.split(" ")[-3]
                        protein_length = line.split(" ")[-2]
                        organism_code_list.add(organism_code)
                    else:
                        if len(line.split(" ")) != 1:
                            domain_code = [x for x in line.split(" ") if x.startswith("PF")][0]
                            domain_start = int(line.split(" ")[-1].split("-")[0])
                            domain_end = int(line.split(" ")[-1].split("-")[1])
                            domain_length = domain_end - domain_start
                            #get only proteins without broken domains
                            if line.split(' ')[-2].split('-')[0].isdigit():
                                broken_protein_list.append((organism_code, protein_code))
                            else:
                                insert_list.append((organism_code, protein_code, 
                                                    uniprot_code, protein_length, 
                                                    domain_code, domain_start, 
                                                    domain_end, domain_length))
            for item in broken_protein_list:
                organism = item[0]
                broken_protein = item[1]
                removal_list = [x for x in insert_list if organism in x and broken_protein in x]
                for item in removal_list:
                    insert_list.remove(item)
            cursor.executemany(add_data, insert_list)
            mariadb_connection.commit()
            
    #Split the table into separated tables for each species
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
