from flask import Flask, render_template, request, redirect, url_for
import httpx
import csv
import os
import pandas as pd
from pymongo import MongoClient
import pubchempy as pcp
import requests
from flask import Blueprint, render_template

form = Blueprint('form', __name__, static_folder='static', template_folder='templates')

# MongoDB Atlas connection string (replace with your own)
connection_string = "mongodb://localhost:27017/"

# Connect to MongoDB Atlas
client = MongoClient(connection_string)
db = client["Drug_Interaction"]
collection = db["interaction_files"]

print("✅ Connected to MongoDB Atlas.")

def get_pubchem_info(drug_name):
    """Fetch PubChem CID and DrugBank ID from PubChem using the drug name."""
    try:
        results = pcp.get_compounds(drug_name, 'name')
        if results:
            compound = results[0]
            pubchem_cid = compound.cid  # PubChem CID
            # Extract DrugBank ID from synonyms
            drugbank_id = None
            if compound.synonyms:
                for synonym in compound.synonyms:
                    if synonym.startswith("DB"):
                        drugbank_id = synonym  # DrugBank ID
                        break
            return pubchem_cid, drugbank_id
    except Exception as e:
        print(f"PubChem lookup failed for {drug_name}: {e}")
    return None, None

def get_chembl_info(drug_name):
    """Fetch ChEMBL ID from ChEMBL using REST API."""
    try:
        url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?pref_name__iexact={drug_name}"
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            data = response.json()
            molecules = data.get("molecules", [])
            if molecules:
                chembl_id = molecules[0].get("molecule_chembl_id")  # ChEMBL ID
                return chembl_id
    except Exception as e:
        print(f"ChEMBL lookup failed for {drug_name}: {e}")
    return None

def get_all_drug_ids(drug_name):
    """Fetch all IDs (PubChem CID, DrugBank ID, and ChEMBL ID) for a given drug name."""
    pubchem_cid, drugbank_id = get_pubchem_info(drug_name)
    chembl_id = get_chembl_info(drug_name)
    
    return {
        "Drug Name": drug_name,
        "PubChem CID": pubchem_cid,
        "DrugBank ID": drugbank_id,
        "ChEMBL ID": chembl_id
    }

def save_to_csv(drug_info, filename="drug_names_ids.csv"):
    """Save drug information to a CSV file."""
    file_exists = os.path.isfile(filename)
    with open(filename, mode="a", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=["Drug Name", "PubChem CID", "DrugBank ID", "ChEMBL ID"])
        if not file_exists:
            writer.writeheader()  # Write header if file doesn't exist
        writer.writerow(drug_info)

def is_drug_in_csv(drug_name, filename="drug_names_ids.csv"):
    """Check if a drug is already in the CSV file."""
    if not os.path.isfile(filename):
        return False
    with open(filename, mode="r", encoding="utf-8") as file:
        reader = csv.DictReader(file)
        for row in reader:
            if row["Drug Name"].lower() == drug_name.lower():
                return True
    return False

def load_similarity_matrix(filename="Drug_data/chem_similarity.csv"):
    """Load the drug similarity matrix from a CSV file."""
    similarity_matrix = pd.read_csv(filename, index_col=0)
    return similarity_matrix

def get_top_similar_drugs(drugbank_id, similarity_matrix, top_n=3):
    """Get the top N similar drugs for a given DrugBank ID."""
    if drugbank_id not in similarity_matrix.index:
        return []
    # Get the similarity scores for the given drug
    similarities = similarity_matrix.loc[drugbank_id]
    # Sort by similarity (descending) and get the top N
    top_drugs = similarities.sort_values(ascending=False).index[1:top_n + 1]  # Skip self-similarity
    return list(top_drugs)

def get_drug_name_from_drugbank_id(drugbank_id):
    """Fetch the drug name from PubChem or ChEMBL using the DrugBank ID."""
    # Try fetching from PubChem
    try:
        results = pcp.get_compounds(drugbank_id, 'name')
        if results:
            compound = results[0]
            if compound.synonyms:
                return compound.synonyms[0]  # Return the first synonym (usually the common name)
    except Exception as e:
        print(f"PubChem lookup failed for {drugbank_id}: {e}")
    
    # Try fetching from ChEMBL
    try:
        url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_chembl_id={drugbank_id}"
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            data = response.json()
            molecules = data.get("molecules", [])
            if molecules:
                return molecules[0].get("pref_name")  # Preferred name
    except Exception as e:
        print(f"ChEMBL lookup failed for {drugbank_id}: {e}")
    
    # If no name is found, return None
    return None

def save_similar_drugs_to_csv(input_drug_name, similar_drugbank_ids, filename="SAME_DRUG.csv"):
    """Save the input drug name and its similar drugs (resolved names) to a CSV file."""
    file_exists = os.path.isfile(filename)
    similar_drug_names = []

    # Resolve drug names for similar DrugBank IDs
    for drugbank_id in similar_drugbank_ids:
        drug_name = get_drug_name_from_drugbank_id(drugbank_id)
        if drug_name:
            similar_drug_names.append(drug_name)
        else:
            print(f"❌ No name found for DrugBank ID {drugbank_id}. Skipping...")

    # Save to CSV
    with open(filename, mode="a", newline="", encoding="utf-8") as file:
        fieldnames = ["Input Drug Name", "Similar Drug 1", "Similar Drug 2", "Similar Drug 3"]
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        if not file_exists:
            writer.writeheader()  # Write header if file doesn't exist
        writer.writerow({
            "Input Drug Name": input_drug_name,
            "Similar Drug 1": similar_drug_names[0] if len(similar_drug_names) > 0 else "",
            "Similar Drug 2": similar_drug_names[1] if len(similar_drug_names) > 1 else "",
            "Similar Drug 3": similar_drug_names[2] if len(similar_drug_names) > 2 else ""
        })

    return similar_drug_names

def generate_combinations_from_same_drug(drug_names, filename="SAME_DRUG.csv"):
    """Generate combinations of drugs based on similar drugs from SAME_DRUG.csv."""
    # Load the SAME_DRUG.csv file
    same_drug_dict = {}
    if not os.path.isfile(filename):
        print(f"❌ {filename} not found. Please ensure the file exists.")
        return []
    
    with open(filename, mode="r", encoding="utf-8") as file:
        reader = csv.DictReader(file)
        for row in reader:
            input_drug = row["Input Drug Name"]
            similar_drugs = [row["Similar Drug 1"], row["Similar Drug 2"], row["Similar Drug 3"]]
            # Remove empty strings (if any)
            similar_drugs = [drug for drug in similar_drugs if drug]
            same_drug_dict[input_drug] = similar_drugs
    
    # Create a list of all groups (input drugs and their similar drugs)
    groups = []
    for drug in drug_names:
        if drug in same_drug_dict:
            group = [drug] + same_drug_dict[drug]  # Input drug + similar drugs
            groups.append(group)
        else:
            print(f"❌ {drug} not found in SAME_DRUG.csv. Skipping...")
    
    # Generate combinations
    combinations = set()  # Use a set to avoid duplicate pairs
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            for drug1 in groups[i]:
                for drug2 in groups[j]:
                    # Ensure no duplicate pairs (order ignored)
                    pair = tuple(sorted([drug1, drug2]))
                    combinations.add(pair)
    
    # Convert the set of tuples to a sorted list of comma-separated strings
    combinations = sorted([f"{pair[0]},{pair[1]}" for pair in combinations])
    
    return combinations

def fetch_drug_data(drug_name):
    """Fetch drug interaction data from PubChem and save it to MongoDB."""
    # Check if the drug data already exists in MongoDB
    existing_document = collection.find_one({"drug_name": drug_name})
    if existing_document:
        print(f"✅ Data for {drug_name} already exists in MongoDB. Skipping...")
        return
    
    url = f"https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={{%22download%22:%22*%22,%22collection%22:%22drugbankddi%22,%22order%22:[%22cid2,asc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22pubchem_name_%5E{drug_name}%24_drugbankddi%22,%22where%22:{{%22ands%22:[{{%22name%22:%22%5E{drug_name}%24%22}}]}}}}"

    try:
        response = httpx.get(url)
        response.raise_for_status()

        # Save the CSV content to MongoDB
        file_name = f"{drug_name}_response.csv"
        document = {
            "drug_name": drug_name,
            "file_name": file_name,
            "content": response.text  # Save the CSV content as text
        }
        collection.insert_one(document)
        print(f"File saved to MongoDB as '{file_name}'")

    except httpx.HTTPStatusError as e:
        print(f"HTTP error occurred: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# def search_interactions_for_drugs(drugs_to_check):
#     """Search for interactions between the given drugs using data from MongoDB."""
#     interactions = set()

#     # Fetch data for all drugs in the combination
#     for drug_name in drugs_to_check:
#         fetch_drug_data(drug_name.strip())

#     # Search for interactions between the drugs
#     for drug_name in drugs_to_check:
#         file_name = f"{drug_name.strip()}_response.csv"
#         document = collection.find_one({"file_name": file_name})
#         if document:
#             content = document["content"]
#             # Parse the CSV content
#             reader = csv.reader(content.splitlines())
#             for row in reader:
#                 drug1 = row[3]
#                 drug2 = row[5]
#                 interaction = row[6]
                
#                 if drug1 in drugs_to_check and drug2 in drugs_to_check:
#                     interaction_tuple = tuple(sorted([drug1, drug2]))
#                     if interaction_tuple not in interactions:
#                         interactions.add(interaction_tuple)
#                         print(f"{drug1} - {drug2}: {interaction}")
#         else:
#             print(f"File for {drug_name.strip()} not found in MongoDB.")
    
#     if not interactions:
#         print("No interactions found between the drugs.")
#     return interactions

def search_interactions_for_drugs(drugs_to_check):
    """Search for interactions between the given drugs using data from MongoDB."""
    interactions = set()

    # Fetch data for all drugs in the combination
    for drug_name in drugs_to_check:
        fetch_drug_data(drug_name.strip())

    # Search for interactions between the drugs
    for drug_name in drugs_to_check:
        file_name = f"{drug_name.strip()}_response.csv"
        document = collection.find_one({"file_name": file_name})
        if document:
            content = document["content"]
            # Parse the CSV content
            reader = csv.DictReader(content.splitlines())
            for row in reader:
                drug1 = row["name"]
                drug2 = row["name2"]
                description = row["descr"]
                
                if drug1 in drugs_to_check and drug2 in drugs_to_check:
                    interaction_tuple = tuple(sorted([drug1, drug2]))
                    if interaction_tuple not in interactions:
                        interactions.add((interaction_tuple, description))
                        print(f"{drug1} - {drug2}: {description}")
        else:
            print(f"File for {drug_name.strip()} not found in MongoDB.")
    
    if not interactions:
        print("No interactions found between the drugs.")
    return interactions

def are_similar_drugs_present(input_drug_name, filename="SAME_DRUG.csv"):
    """Check if similar drugs for the input drug are already present in SAME_DRUG.csv."""
    if not os.path.isfile(filename):
        return False  # File doesn't exist, so similar drugs are not present
    
    with open(filename, mode="r", encoding="utf-8") as file:
        reader = csv.DictReader(file)
        for row in reader:
            if row["Input Drug Name"].lower() == input_drug_name.lower():
                # Check if at least one similar drug is present
                if row["Similar Drug 1"] or row["Similar Drug 2"] or row["Similar Drug 3"]:
                    return True  # Similar drugs are already present
        return False  # Input drug not found or no similar drugs present

def get_similar_drugs_from_csv(input_drug_name, filename="SAME_DRUG.csv"):
    """Retrieve similar drugs for a given input drug from SAME_DRUG.csv."""
    similar_drugs = []
    if not os.path.isfile(filename):
        return similar_drugs
    
    with open(filename, mode="r", encoding="utf-8") as file:
        reader = csv.DictReader(file)
        for row in reader:
            if row["Input Drug Name"].lower() == input_drug_name.lower():
                similar_drugs = [row["Similar Drug 1"], row["Similar Drug 2"], row["Similar Drug 3"]]
                # Remove empty strings (if any)
                similar_drugs = [drug for drug in similar_drugs if drug]
                break
    return similar_drugs        

# @form.route("/", methods=["GET", "POST"])
# def index():
#     if request.method == "POST":
#         # Step 1: Take user input
#         input_drug_names = request.form.get("drug_names").strip().split()
        
#         # Step 2: Check if similar drugs are already present in SAME_DRUG.csv
#         similarity_matrix = load_similarity_matrix()
#         similar_drugs_info = {}
#         for drug_name in input_drug_names:
#             if are_similar_drugs_present(drug_name):  # Check if similar drugs are already present
#                 print(f"✅ Similar drugs for {drug_name} already exist in SAME_DRUG.csv. Skipping...")
#             else:
#                 print(f"❌ Similar drugs for {drug_name} not found in SAME_DRUG.csv. Fetching...")
#                 # Fetch DrugBank ID for the input drug
#                 drug_info = get_all_drug_ids(drug_name)
#                 drugbank_id = drug_info["DrugBank ID"]
#                 if drugbank_id:
#                     # Fetch similar DrugBank IDs from the similarity matrix
#                     similar_drugbank_ids = get_top_similar_drugs(drugbank_id, similarity_matrix)
#                     similar_drug_names = []
                    
#                     # Resolve DrugBank IDs to drug names
#                     for drugbank_id in similar_drugbank_ids:
#                         drug_name_resolved = get_drug_name_from_drugbank_id(drugbank_id)
#                         if drug_name_resolved:
#                             similar_drug_names.append(drug_name_resolved)
#                             print(f"✅ Resolved name for DrugBank ID {drugbank_id}: {drug_name_resolved}")
#                         else:
#                             print(f"❌ No name resolved for DrugBank ID {drugbank_id}. Skipping...")
                    
#                     # Save similar drug names to SAME_DRUG.csv
#                     if similar_drug_names:  # Only save if at least one similar drug name was resolved
#                         save_similar_drugs_to_csv(drug_name, similar_drug_names)
#                         similar_drugs_info[drug_name] = similar_drug_names
#                         print(f"✅ Similar drugs for {drug_name}: {similar_drug_names}")
#                     else:
#                         print(f"❌ No similar drugs resolved for {drug_name}. Skipping...")
#                 else:
#                     print(f"❌ No DrugBank ID found for {drug_name}.")
        
#         # Step 3: Generate combinations
#         combinations = generate_combinations_from_same_drug(input_drug_names)
#         if combinations:
#             print("\n✅ Generated Combinations and Interactions:")
#             interactions = []
#             for combo in combinations:
#                 print(f"\n🔍 Checking interactions for: {combo}")
#                 drugs_to_check = [drug.strip() for drug in combo.split(",")]
#                 interaction_results = search_interactions_for_drugs(drugs_to_check)
#                 if interaction_results:
#                     interactions.append((combo, interaction_results))
#             return render_template("results.html", input_drug_names=input_drug_names, similar_drugs_info=similar_drugs_info, combinations=combinations, interactions=interactions)
#         else:
#             print("❌ No combinations generated.")
#             return render_template("results.html", input_drug_names=input_drug_names, similar_drugs_info=similar_drugs_info, combinations=[], interactions=[])
    
#     return render_template("form.html")

#finailllll---------------------------------------------------------------down

# @form.route("/", methods=["GET", "POST"])
# def index():
#     if request.method == "POST":
#         # Step 1: Take user input
#         input_drug_names = request.form.get("drug_names").strip().split()
        
#         # Step 2: Check if similar drugs are already present in SAME_DRUG.csv
#         similar_drugs_info = {}
#         for drug_name in input_drug_names:
#             if are_similar_drugs_present(drug_name):  # Check if similar drugs are already present
#                 print(f"✅ Similar drugs for {drug_name} already exist in SAME_DRUG.csv. Skipping...")
#                 similar_drugs_info[drug_name] = get_similar_drugs_from_csv(drug_name)
#             else:
#                 print(f"❌ Similar drugs for {drug_name} not found in SAME_DRUG.csv. Fetching...")
#                 # Fetch DrugBank ID for the input drug
#                 drug_info = get_all_drug_ids(drug_name)
#                 drugbank_id = drug_info["DrugBank ID"]
#                 if drugbank_id:
#                     # Fetch similar DrugBank IDs from the similarity matrix
#                     similarity_matrix = load_similarity_matrix()
#                     similar_drugbank_ids = get_top_similar_drugs(drugbank_id, similarity_matrix)
#                     similar_drug_names = []
                    
#                     # Resolve DrugBank IDs to drug names
#                     for drugbank_id in similar_drugbank_ids:
#                         drug_name_resolved = get_drug_name_from_drugbank_id(drugbank_id)
#                         if drug_name_resolved:
#                             similar_drug_names.append(drug_name_resolved)
#                         else:
#                             print(f"❌ No name resolved for DrugBank ID {drugbank_id}. Skipping...")
                    
#                     # Save similar drug names to SAME_DRUG.csv
#                     if similar_drug_names:  # Only save if at least one similar drug name was resolved
#                         save_similar_drugs_to_csv(drug_name, similar_drug_names)
#                         similar_drugs_info[drug_name] = similar_drug_names
#                         print(f"✅ Similar drugs for {drug_name}: {similar_drug_names}")
#                     else:
#                         print(f"❌ No similar drugs resolved for {drug_name}. Skipping...")
#                 else:
#                     print(f"❌ No DrugBank ID found for {drug_name}.")
        
#         # Step 3: Generate combinations
#         combinations = generate_combinations_from_same_drug(input_drug_names)
#         if combinations:
#             print("\n✅ Generated Combinations and Interactions:")
#             interactions = []
#             for combo in combinations:
#                 print(f"\n🔍 Checking interactions for: {combo}")
#                 drugs_to_check = [drug.strip() for drug in combo.split(",")]
#                 interaction_results = search_interactions_for_drugs(drugs_to_check)
#                 if interaction_results:
#                     interactions.append((combo, interaction_results))
#             return render_template("results.html", input_drug_names=input_drug_names, similar_drugs_info=similar_drugs_info, combinations=combinations, interactions=interactions)
#         else:
#             print("❌ No combinations generated.")
#             return render_template("results.html", input_drug_names=input_drug_names, similar_drugs_info=similar_drugs_info, combinations=[], interactions=[])
    
#     return render_template("form.html")


@form.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        # Step 1: Take user input
        input_drug_names = request.form.get("drug_names").strip().split()
        
        # Step 2: Check if similar drugs are already present in SAME_DRUG.csv
        similar_drugs_info = {}
        for drug_name in input_drug_names:
            if are_similar_drugs_present(drug_name):  # Check if similar drugs are already present
                print(f"✅ Similar drugs for {drug_name} already exist in SAME_DRUG.csv. Skipping...")
                similar_drugs_info[drug_name] = get_similar_drugs_from_csv(drug_name)
            else:
                print(f"❌ Similar drugs for {drug_name} not found in SAME_DRUG.csv. Fetching...")
                # Fetch DrugBank ID for the input drug
                drug_info = get_all_drug_ids(drug_name)
                drugbank_id = drug_info["DrugBank ID"]
                if drugbank_id:
                    # Fetch similar DrugBank IDs from the similarity matrix
                    similarity_matrix = load_similarity_matrix()
                    similar_drugbank_ids = get_top_similar_drugs(drugbank_id, similarity_matrix)
                    similar_drug_names = []
                    
                    # Resolve DrugBank IDs to drug names
                    for drugbank_id in similar_drugbank_ids:
                        drug_name_resolved = get_drug_name_from_drugbank_id(drugbank_id)
                        if drug_name_resolved:
                            similar_drug_names.append(drug_name_resolved)
                        else:
                            print(f"❌ No name resolved for DrugBank ID {drugbank_id}. Skipping...")
                    
                    # Save similar drug names to SAME_DRUG.csv
                    if similar_drug_names:  # Only save if at least one similar drug name was resolved
                        save_similar_drugs_to_csv(drug_name, similar_drug_names)
                        similar_drugs_info[drug_name] = similar_drug_names
                        print(f"✅ Similar drugs for {drug_name}: {similar_drug_names}")
                    else:
                        print(f"❌ No similar drugs resolved for {drug_name}. Skipping...")
                else:
                    print(f"❌ No DrugBank ID found for {drug_name}.")
        
        # Step 3: Generate combinations including original input drugs
        if all(drug_name in similar_drugs_info for drug_name in input_drug_names):
            # Include original input drugs in the list of drugs to combine
            all_drugs = input_drug_names.copy()
            for drug_name in input_drug_names:
                all_drugs.extend(similar_drugs_info[drug_name])
            
            # Generate all possible combinations
            combinations = set()
            for i in range(len(all_drugs)):
                for j in range(i + 1, len(all_drugs)):
                    pair = tuple(sorted([all_drugs[i], all_drugs[j]]))
                    combinations.add(pair)
            
            # Convert the set of tuples to a sorted list of comma-separated strings
            combinations = sorted([f"{pair[0]},{pair[1]}" for pair in combinations])
            
            if combinations:
                print("\n✅ Generated Combinations and Interactions:")
                interactions = []
                for combo in combinations:
                    print(f"\n🔍 Checking interactions for: {combo}")
                    drugs_to_check = [drug.strip() for drug in combo.split(",")]
                    interaction_results = search_interactions_for_drugs(drugs_to_check)
                    if interaction_results:
                        interactions.append((combo, interaction_results))
                return render_template("results.html", input_drug_names=input_drug_names, similar_drugs_info=similar_drugs_info, combinations=combinations, interactions=interactions)
            else:
                print("❌ No combinations generated.")
                return render_template("results.html", input_drug_names=input_drug_names, similar_drugs_info=similar_drugs_info, combinations=[], interactions=[])
        else:
            print("❌ Unable to generate combinations. Some input drugs do not have similar drugs resolved.")
            return render_template("results.html", input_drug_names=input_drug_names, similar_drugs_info=similar_drugs_info, combinations=[], interactions=[])
    
    return render_template("form.html")
#---------------------------------------------------------------------------------------------------------------------------------------------------
# from flask import Flask, render_template, request, redirect, url_for
# import httpx
# import csv
# import os
# import pandas as pd
# from pymongo import MongoClient
# import pubchempy as pcp
# import requests
# from flask import Blueprint, render_template



# form = Blueprint('form', __name__, static_folder='static', template_folder='templates')

# # MongoDB Atlas connection string (replace with your own)
# connection_string = "mongodb://localhost:27017/"

# # Connect to MongoDB Atlas
# client = MongoClient(connection_string)
# db = client["Drug_Interaction"]
# collection = db["interaction_files"]

# print("✅ Connected to MongoDB Atlas.")

# def get_pubchem_info(drug_name):
#     """Fetch PubChem CID and DrugBank ID from PubChem using the drug name."""
#     try:
#         results = pcp.get_compounds(drug_name, 'name')
#         if results:
#             compound = results[0]
#             pubchem_cid = compound.cid  # PubChem CID
#             # Extract DrugBank ID from synonyms
#             drugbank_id = None
#             if compound.synonyms:
#                 for synonym in compound.synonyms:
#                     if synonym.startswith("DB"):
#                         drugbank_id = synonym  # DrugBank ID
#                         break
#             return pubchem_cid, drugbank_id
#     except Exception as e:
#         print(f"PubChem lookup failed for {drug_name}: {e}")
#     return None, None

# def get_chembl_info(drug_name):
#     """Fetch ChEMBL ID from ChEMBL using REST API."""
#     try:
#         url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?pref_name__iexact={drug_name}"
#         response = requests.get(url, timeout=5)
#         if response.status_code == 200:
#             data = response.json()
#             molecules = data.get("molecules", [])
#             if molecules:
#                 chembl_id = molecules[0].get("molecule_chembl_id")  # ChEMBL ID
#                 return chembl_id
#     except Exception as e:
#         print(f"ChEMBL lookup failed for {drug_name}: {e}")
#     return None

# def get_all_drug_ids(drug_name):
#     """Fetch all IDs (PubChem CID, DrugBank ID, and ChEMBL ID) for a given drug name."""
#     pubchem_cid, drugbank_id = get_pubchem_info(drug_name)
#     chembl_id = get_chembl_info(drug_name)
    
#     return {
#         "Drug Name": drug_name,
#         "PubChem CID": pubchem_cid,
#         "DrugBank ID": drugbank_id,
#         "ChEMBL ID": chembl_id
#     }

# def save_to_csv(drug_info, filename="drug_names_ids.csv"):
#     """Save drug information to a CSV file."""
#     file_exists = os.path.isfile(filename)
#     with open(filename, mode="a", newline="", encoding="utf-8") as file:
#         writer = csv.DictWriter(file, fieldnames=["Drug Name", "PubChem CID", "DrugBank ID", "ChEMBL ID"])
#         if not file_exists:
#             writer.writeheader()  # Write header if file doesn't exist
#         writer.writerow(drug_info)

# def is_drug_in_csv(drug_name, filename="drug_names_ids.csv"):
#     """Check if a drug is already in the CSV file."""
#     if not os.path.isfile(filename):
#         return False
#     with open(filename, mode="r", encoding="utf-8") as file:
#         reader = csv.DictReader(file)
#         for row in reader:
#             if row["Drug Name"].lower() == drug_name.lower():
#                 return True
#     return False

# def load_similarity_matrix(filename="Drug_data\chem_similarity.csv"):
#     """Load the drug similarity matrix from a CSV file."""
#     similarity_matrix = pd.read_csv(filename, index_col=0)
#     return similarity_matrix

# def get_top_similar_drugs(drugbank_id, similarity_matrix, top_n=3):
#     """Get the top N similar drugs for a given DrugBank ID."""
#     if drugbank_id not in similarity_matrix.index:
#         return []
#     # Get the similarity scores for the given drug
#     similarities = similarity_matrix.loc[drugbank_id]
#     # Sort by similarity (descending) and get the top N
#     top_drugs = similarities.sort_values(ascending=False).index[1:top_n + 1]  # Skip self-similarity
#     return list(top_drugs)

# def get_drug_name_from_drugbank_id(drugbank_id):
#     """Fetch the drug name from PubChem or ChEMBL using the DrugBank ID."""
#     # Try fetching from PubChem
#     try:
#         results = pcp.get_compounds(drugbank_id, 'name')
#         if results:
#             compound = results[0]
#             if compound.synonyms:
#                 return compound.synonyms[0]  # Return the first synonym (usually the common name)
#     except Exception as e:
#         print(f"PubChem lookup failed for {drugbank_id}: {e}")
    
#     # Try fetching from ChEMBL
#     try:
#         url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_chembl_id={drugbank_id}"
#         response = requests.get(url, timeout=5)
#         if response.status_code == 200:
#             data = response.json()
#             molecules = data.get("molecules", [])
#             if molecules:
#                 return molecules[0].get("pref_name")  # Preferred name
#     except Exception as e:
#         print(f"ChEMBL lookup failed for {drugbank_id}: {e}")
    
#     # If no name is found, return None
#     return None

# def save_similar_drugs_to_csv(input_drug_name, similar_drugbank_ids, filename="SAME_DRUG.csv"):
#     """Save the input drug name and its similar drugs (resolved names) to a CSV file."""
#     file_exists = os.path.isfile(filename)
#     similar_drug_names = []

#     # Resolve drug names for similar DrugBank IDs
#     for drugbank_id in similar_drugbank_ids:
#         drug_name = get_drug_name_from_drugbank_id(drugbank_id)
#         if drug_name:
#             similar_drug_names.append(drug_name)
#         else:
#             print(f"❌ No name found for DrugBank ID {drugbank_id}. Skipping...")

#     # Save to CSV
#     with open(filename, mode="a", newline="", encoding="utf-8") as file:
#         fieldnames = ["Input Drug Name", "Similar Drug 1", "Similar Drug 2", "Similar Drug 3"]
#         writer = csv.DictWriter(file, fieldnames=fieldnames)
#         if not file_exists:
#             writer.writeheader()  # Write header if file doesn't exist
#         writer.writerow({
#             "Input Drug Name": input_drug_name,
#             "Similar Drug 1": similar_drug_names[0] if len(similar_drug_names) > 0 else "",
#             "Similar Drug 2": similar_drug_names[1] if len(similar_drug_names) > 1 else "",
#             "Similar Drug 3": similar_drug_names[2] if len(similar_drug_names) > 2 else ""
#         })

#     return similar_drug_names

# def generate_combinations_from_same_drug(drug_names, filename="SAME_DRUG.csv"):
#     """Generate combinations of drugs based on similar drugs from SAME_DRUG.csv."""
#     # Load the SAME_DRUG.csv file
#     same_drug_dict = {}
#     if not os.path.isfile(filename):
#         print(f"❌ {filename} not found. Please ensure the file exists.")
#         return []
    
#     with open(filename, mode="r", encoding="utf-8") as file:
#         reader = csv.DictReader(file)
#         for row in reader:
#             input_drug = row["Input Drug Name"]
#             similar_drugs = [row["Similar Drug 1"], row["Similar Drug 2"], row["Similar Drug 3"]]
#             # Remove empty strings (if any)
#             similar_drugs = [drug for drug in similar_drugs if drug]
#             same_drug_dict[input_drug] = similar_drugs
    
#     # Create a list of all groups (input drugs and their similar drugs)
#     groups = []
#     for drug in drug_names:
#         if drug in same_drug_dict:
#             group = [drug] + same_drug_dict[drug]  # Input drug + similar drugs
#             groups.append(group)
#         else:
#             print(f"❌ {drug} not found in SAME_DRUG.csv. Skipping...")
    
#     # Generate combinations
#     combinations = set()  # Use a set to avoid duplicate pairs
#     for i in range(len(groups)):
#         for j in range(i + 1, len(groups)):
#             for drug1 in groups[i]:
#                 for drug2 in groups[j]:
#                     # Ensure no duplicate pairs (order ignored)
#                     pair = tuple(sorted([drug1, drug2]))
#                     combinations.add(pair)
    
#     # Convert the set of tuples to a sorted list of comma-separated strings
#     combinations = sorted([f"{pair[0]},{pair[1]}" for pair in combinations])
    
#     return combinations

# def fetch_drug_data(drug_name):
#     """Fetch drug interaction data from PubChem and save it to MongoDB."""
#     url = f"https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={{%22download%22:%22*%22,%22collection%22:%22drugbankddi%22,%22order%22:[%22cid2,asc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22pubchem_name_%5E{drug_name}%24_drugbankddi%22,%22where%22:{{%22ands%22:[{{%22name%22:%22%5E{drug_name}%24%22}}]}}}}"

#     try:
#         response = httpx.get(url)
#         response.raise_for_status()

#         # Save the CSV content to MongoDB
#         file_name = f"{drug_name}_response.csv"
#         document = {
#             "drug_name": drug_name,
#             "file_name": file_name,
#             "content": response.text  # Save the CSV content as text
#         }
#         collection.insert_one(document)
#         print(f"File saved to MongoDB as '{file_name}'")

#     except httpx.HTTPStatusError as e:
#         print(f"HTTP error occurred: {e}")
#     except Exception as e:
#         print(f"An unexpected error occurred: {e}")

# def search_interactions_for_drugs(drugs_to_check):
#     """Search for interactions between the given drugs using data from MongoDB."""
#     interactions = set()

#     # Fetch data for all drugs in the combination
#     for drug_name in drugs_to_check:
#         fetch_drug_data(drug_name.strip())

#     # Search for interactions between the drugs
#     for drug_name in drugs_to_check:
#         file_name = f"{drug_name.strip()}_response.csv"
#         document = collection.find_one({"file_name": file_name})
#         if document:
#             content = document["content"]
#             # Parse the CSV content
#             reader = csv.reader(content.splitlines())
#             for row in reader:
#                 drug1 = row[3]
#                 drug2 = row[5]
#                 interaction = row[6]
                
#                 if drug1 in drugs_to_check and drug2 in drugs_to_check:
#                     interaction_tuple = tuple(sorted([drug1, drug2]))
#                     if interaction_tuple not in interactions:
#                         interactions.add(interaction_tuple)
#                         print(f"{drug1} - {drug2}: {interaction}")
#         else:
#             print(f"File for {drug_name.strip()} not found in MongoDB.")
    
#     if not interactions:
#         print("No interactions found between the drugs.")
#     return interactions

# @form.route("/", methods=["GET", "POST"])
# def index():
#     if request.method == "POST":
#         # Step 1: Take user input
#         input_drug_names = request.form.get("drug_names").strip().split()
        
#         # Step 2: Fetch IDs and save to CSV
#         for drug_name in input_drug_names:
#             if not is_drug_in_csv(drug_name):
#                 drug_info = get_all_drug_ids(drug_name)
#                 if drug_info["PubChem CID"] or drug_info["DrugBank ID"] or drug_info["ChEMBL ID"]:
#                     save_to_csv(drug_info)
#                     print(f"✅ Information for {drug_name} saved to drug_names_ids.csv.")
#                 else:
#                     print(f"❌ No information found for {drug_name}.")
#             else:
#                 print(f"✅ {drug_name} already exists in drug_names_ids.csv. Skipping...")
        
#         # Step 3: Fetch similar drugs, resolve their names, and save to SAME_DRUG.csv
#         similarity_matrix = load_similarity_matrix()
#         similar_drugs_info = {}
#         for drug_name in input_drug_names:
#             drug_info = get_all_drug_ids(drug_name)
#             drugbank_id = drug_info["DrugBank ID"]
#             if drugbank_id:
#                 similar_drugbank_ids = get_top_similar_drugs(drugbank_id, similarity_matrix)
#                 similar_drug_names = save_similar_drugs_to_csv(drug_name, similar_drugbank_ids)
#                 similar_drugs_info[drug_name] = similar_drug_names
#                 print(f"✅ Similar drugs for {drug_name}: {similar_drug_names}")
#             else:
#                 print(f"❌ No DrugBank ID found for {drug_name}.")
        
#         # Step 4: Generate combinations
#         combinations = generate_combinations_from_same_drug(input_drug_names)
#         if combinations:
#             print("\n✅ Generated Combinations and Interactions:")
#             interactions = []
#             for combo in combinations:
#                 print(f"\n🔍 Checking interactions for: {combo}")
#                 drugs_to_check = [drug.strip() for drug in combo.split(",")]
#                 interaction_results = search_interactions_for_drugs(drugs_to_check)
#                 if interaction_results:
#                     interactions.append((combo, interaction_results))
#             return render_template("results.html", input_drug_names=input_drug_names, similar_drugs_info=similar_drugs_info, combinations=combinations, interactions=interactions)
#         else:
#             print("❌ No combinations generated.")
#             return render_template("results.html", input_drug_names=input_drug_names, similar_drugs_info=similar_drugs_info, combinations=[], interactions=[])
    
#     return render_template("form.html")

# from flask import Flask, render_template, request, redirect, url_for
# import httpx
# import csv
# import os
# import pandas as pd
# from pymongo import MongoClient
# import pubchempy as pcp
# import requests
# from flask import Blueprint, render_template

# form = Blueprint('form', __name__, static_folder='static', template_folder='templates')

# # MongoDB Atlas connection string (replace with your own)
# connection_string = "mongodb://localhost:27017/"

# # Connect to MongoDB Atlas
# client = MongoClient(connection_string)
# db = client["Drug_Interaction"]
# collection = db["interaction_files"]

# print("✅ Connected to MongoDB Atlas.")

# def get_pubchem_info(drug_name):
#     """Fetch PubChem CID and DrugBank ID from PubChem using the drug name."""
#     try:
#         results = pcp.get_compounds(drug_name, 'name')
#         if results:
#             compound = results[0]
#             pubchem_cid = compound.cid  # PubChem CID
#             # Extract DrugBank ID from synonyms
#             drugbank_id = None
#             if compound.synonyms:
#                 for synonym in compound.synonyms:
#                     if synonym.startswith("DB"):
#                         drugbank_id = synonym  # DrugBank ID
#                         break
#             return pubchem_cid, drugbank_id
#     except Exception as e:
#         print(f"PubChem lookup failed for {drug_name}: {e}")
#     return None, None

# def get_chembl_info(drug_name):
#     """Fetch ChEMBL ID from ChEMBL using REST API."""
#     try:
#         url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?pref_name__iexact={drug_name}"
#         response = requests.get(url, timeout=5)
#         if response.status_code == 200:
#             data = response.json()
#             molecules = data.get("molecules", [])
#             if molecules:
#                 chembl_id = molecules[0].get("molecule_chembl_id")  # ChEMBL ID
#                 return chembl_id
#     except Exception as e:
#         print(f"ChEMBL lookup failed for {drug_name}: {e}")
#     return None

# def get_all_drug_ids(drug_name):
#     """Fetch all IDs (PubChem CID, DrugBank ID, and ChEMBL ID) for a given drug name."""
#     pubchem_cid, drugbank_id = get_pubchem_info(drug_name)
#     chembl_id = get_chembl_info(drug_name)
    
#     return {
#         "Drug Name": drug_name,
#         "PubChem CID": pubchem_cid,
#         "DrugBank ID": drugbank_id,
#         "ChEMBL ID": chembl_id
#     }

# def save_to_csv(drug_info, filename="drug_names_ids.csv"):
#     """Save drug information to a CSV file."""
#     file_exists = os.path.isfile(filename)
#     with open(filename, mode="a", newline="", encoding="utf-8") as file:
#         writer = csv.DictWriter(file, fieldnames=["Drug Name", "PubChem CID", "DrugBank ID", "ChEMBL ID"])
#         if not file_exists:
#             writer.writeheader()  # Write header if file doesn't exist
#         writer.writerow(drug_info)

# def is_drug_in_csv(drug_name, filename="drug_names_ids.csv"):
#     """Check if a drug is already in the CSV file."""
#     if not os.path.isfile(filename):
#         return False
#     with open(filename, mode="r", encoding="utf-8") as file:
#         reader = csv.DictReader(file)
#         for row in reader:
#             if row["Drug Name"].lower() == drug_name.lower():
#                 return True
#     return False

# def load_similarity_matrix(filename="Drug_data/chem_similarity.csv"):
#     """Load the drug similarity matrix from a CSV file."""
#     similarity_matrix = pd.read_csv(filename, index_col=0)
#     return similarity_matrix

# def get_top_similar_drugs(drugbank_id, similarity_matrix, top_n=3):
#     """Get the top N similar drugs for a given DrugBank ID."""
#     if drugbank_id not in similarity_matrix.index:
#         return []
#     # Get the similarity scores for the given drug
#     similarities = similarity_matrix.loc[drugbank_id]
#     # Sort by similarity (descending) and get the top N
#     top_drugs = similarities.sort_values(ascending=False).index[1:top_n + 1]  # Skip self-similarity
#     return list(top_drugs)

# def get_drug_name_from_drugbank_id(drugbank_id):
#     """Fetch the drug name from PubChem or ChEMBL using the DrugBank ID."""
#     # Try fetching from PubChem
#     try:
#         results = pcp.get_compounds(drugbank_id, 'name')
#         if results:
#             compound = results[0]
#             if compound.synonyms:
#                 return compound.synonyms[0]  # Return the first synonym (usually the common name)
#     except Exception as e:
#         print(f"PubChem lookup failed for {drugbank_id}: {e}")
    
#     # Try fetching from ChEMBL
#     try:
#         url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_chembl_id={drugbank_id}"
#         response = requests.get(url, timeout=5)
#         if response.status_code == 200:
#             data = response.json()
#             molecules = data.get("molecules", [])
#             if molecules:
#                 return molecules[0].get("pref_name")  # Preferred name
#     except Exception as e:
#         print(f"ChEMBL lookup failed for {drugbank_id}: {e}")
    
#     # If no name is found, return None
#     return None

# def save_similar_drugs_to_csv(input_drug_name, similar_drugbank_ids, filename="SAME_DRUG.csv"):
#     """Save the input drug name and its similar drugs (resolved names) to a CSV file."""
#     file_exists = os.path.isfile(filename)
#     similar_drug_names = []

#     # Resolve drug names for similar DrugBank IDs
#     for drugbank_id in similar_drugbank_ids:
#         drug_name = get_drug_name_from_drugbank_id(drugbank_id)
#         if drug_name:
#             similar_drug_names.append(drug_name)
#         else:
#             print(f"❌ No name found for DrugBank ID {drugbank_id}. Skipping...")

#     # Save to CSV
#     with open(filename, mode="a", newline="", encoding="utf-8") as file:
#         fieldnames = ["Input Drug Name", "Similar Drug 1", "Similar Drug 2", "Similar Drug 3"]
#         writer = csv.DictWriter(file, fieldnames=fieldnames)
#         if not file_exists:
#             writer.writeheader()  # Write header if file doesn't exist
#         writer.writerow({
#             "Input Drug Name": input_drug_name,
#             "Similar Drug 1": similar_drug_names[0] if len(similar_drug_names) > 0 else "",
#             "Similar Drug 2": similar_drug_names[1] if len(similar_drug_names) > 1 else "",
#             "Similar Drug 3": similar_drug_names[2] if len(similar_drug_names) > 2 else ""
#         })

#     return similar_drug_names

# def generate_combinations_from_same_drug(drug_names, filename="SAME_DRUG.csv"):
#     """Generate combinations of drugs based on similar drugs from SAME_DRUG.csv."""
#     # Load the SAME_DRUG.csv file
#     same_drug_dict = {}
#     if not os.path.isfile(filename):
#         print(f"❌ {filename} not found. Please ensure the file exists.")
#         return []
    
#     with open(filename, mode="r", encoding="utf-8") as file:
#         reader = csv.DictReader(file)
#         for row in reader:
#             input_drug = row["Input Drug Name"]
#             similar_drugs = [row["Similar Drug 1"], row["Similar Drug 2"], row["Similar Drug 3"]]
#             # Remove empty strings (if any)
#             similar_drugs = [drug for drug in similar_drugs if drug]
#             same_drug_dict[input_drug] = similar_drugs
    
#     # Create a list of all groups (input drugs and their similar drugs)
#     groups = []
#     for drug in drug_names:
#         if drug in same_drug_dict:
#             group = [drug] + same_drug_dict[drug]  # Input drug + similar drugs
#             groups.append(group)
#         else:
#             print(f"❌ {drug} not found in SAME_DRUG.csv. Skipping...")
    
#     # Generate combinations
#     combinations = set()  # Use a set to avoid duplicate pairs
#     for i in range(len(groups)):
#         for j in range(i + 1, len(groups)):
#             for drug1 in groups[i]:
#                 for drug2 in groups[j]:
#                     # Ensure no duplicate pairs (order ignored)
#                     pair = tuple(sorted([drug1, drug2]))
#                     combinations.add(pair)
    
#     # Convert the set of tuples to a sorted list of comma-separated strings
#     combinations = sorted([f"{pair[0]},{pair[1]}" for pair in combinations])
    
#     return combinations

# def fetch_drug_data(drug_name):
#     """Fetch drug interaction data from PubChem and save it to MongoDB."""
#     # Check if the drug data already exists in MongoDB
#     existing_document = collection.find_one({"drug_name": drug_name})
#     if existing_document:
#         print(f"✅ Data for {drug_name} already exists in MongoDB. Skipping...")
#         return
    
#     url = f"https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={{%22download%22:%22*%22,%22collection%22:%22drugbankddi%22,%22order%22:[%22cid2,asc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22pubchem_name_%5E{drug_name}%24_drugbankddi%22,%22where%22:{{%22ands%22:[{{%22name%22:%22%5E{drug_name}%24%22}}]}}}}"

#     try:
#         response = httpx.get(url)
#         response.raise_for_status()

#         # Save the CSV content to MongoDB
#         file_name = f"{drug_name}_response.csv"
#         document = {
#             "drug_name": drug_name,
#             "file_name": file_name,
#             "content": response.text  # Save the CSV content as text
#         }
#         collection.insert_one(document)
#         print(f"File saved to MongoDB as '{file_name}'")

#     except httpx.HTTPStatusError as e:
#         print(f"HTTP error occurred: {e}")
#     except Exception as e:
#         print(f"An unexpected error occurred: {e}")

# def search_interactions_for_drugs(drugs_to_check):
#     """Search for interactions between the given drugs using data from MongoDB."""
#     interactions = set()

#     # Fetch data for all drugs in the combination
#     for drug_name in drugs_to_check:
#         fetch_drug_data(drug_name.strip())

#     # Search for interactions between the drugs
#     for drug_name in drugs_to_check:
#         file_name = f"{drug_name.strip()}_response.csv"
#         document = collection.find_one({"file_name": file_name})
#         if document:
#             content = document["content"]
#             # Parse the CSV content
#             reader = csv.reader(content.splitlines())
#             for row in reader:
#                 drug1 = row[3]
#                 drug2 = row[5]
#                 interaction = row[6]
                
#                 if drug1 in drugs_to_check and drug2 in drugs_to_check:
#                     interaction_tuple = tuple(sorted([drug1, drug2]))
#                     if interaction_tuple not in interactions:
#                         interactions.add(interaction_tuple)
#                         print(f"{drug1} - {drug2}: {interaction}")
#         else:
#             print(f"File for {drug_name.strip()} not found in MongoDB.")
    
#     if not interactions:
#         print("No interactions found between the drugs.")
#     return interactions

# @form.route("/", methods=["GET", "POST"])
# def index():
#     if request.method == "POST":
#         # Step 1: Take user input
#         input_drug_names = request.form.get("drug_names").strip().split()
        
#         # Step 2: Fetch IDs and save to CSV
#         for drug_name in input_drug_names:
#             if not is_drug_in_csv(drug_name):
#                 drug_info = get_all_drug_ids(drug_name)
#                 if drug_info["PubChem CID"] or drug_info["DrugBank ID"] or drug_info["ChEMBL ID"]:
#                     save_to_csv(drug_info)
#                     print(f"✅ Information for {drug_name} saved to drug_names_ids.csv.")
#                 else:
#                     print(f"❌ No information found for {drug_name}.")
#             else:
#                 print(f"✅ {drug_name} already exists in drug_names_ids.csv. Skipping...")
        
#         # Step 3: Fetch similar drugs, resolve their names, and save to SAME_DRUG.csv
#         similarity_matrix = load_similarity_matrix()
#         similar_drugs_info = {}
#         for drug_name in input_drug_names:
#             drug_info = get_all_drug_ids(drug_name)
#             drugbank_id = drug_info["DrugBank ID"]
#             if drugbank_id:
#                 similar_drugbank_ids = get_top_similar_drugs(drugbank_id, similarity_matrix)
#                 similar_drug_names = save_similar_drugs_to_csv(drug_name, similar_drugbank_ids)
#                 similar_drugs_info[drug_name] = similar_drug_names
#                 print(f"✅ Similar drugs for {drug_name}: {similar_drug_names}")
#             else:
#                 print(f"❌ No DrugBank ID found for {drug_name}.")
        
#         # Step 4: Generate combinations
#         combinations = generate_combinations_from_same_drug(input_drug_names)
#         if combinations:
#             print("\n✅ Generated Combinations and Interactions:")
#             interactions = []
#             for combo in combinations:
#                 print(f"\n🔍 Checking interactions for: {combo}")
#                 drugs_to_check = [drug.strip() for drug in combo.split(",")]
#                 interaction_results = search_interactions_for_drugs(drugs_to_check)
#                 if interaction_results:
#                     interactions.append((combo, interaction_results))
#             return render_template("results.html", input_drug_names=input_drug_names, similar_drugs_info=similar_drugs_info, combinations=combinations, interactions=interactions)
#         else:
#             print("❌ No combinations generated.")
#             return render_template("results.html", input_drug_names=input_drug_names, similar_drugs_info=similar_drugs_info, combinations=[], interactions=[])
    
#     return render_template("form.html")

# from flask import Flask, render_template, request
# import httpx
# import csv
# import os
# import pandas as pd
# from pymongo import MongoClient
# import pubchempy as pcp
# import requests
# from flask import Blueprint

# # Initialize Flask Blueprint
# form = Blueprint('form', __name__, static_folder='static', template_folder='templates')

# # MongoDB connection
# connection_string = "mongodb://localhost:27017/"
# client = MongoClient(connection_string)
# db = client["Drug_Interaction"]
# collection = db["interaction_files"]

# print("✅ Connected to MongoDB.")

# # Helper Functions

# def get_pubchem_info(drug_name):
#     """Fetch PubChem CID and DrugBank ID from PubChem using the drug name."""
#     try:
#         results = pcp.get_compounds(drug_name, 'name')
#         if results:
#             compound = results[0]
#             pubchem_cid = compound.cid
#             drugbank_id = None
#             if compound.synonyms:
#                 for synonym in compound.synonyms:
#                     if synonym.startswith("DB"):
#                         drugbank_id = synonym
#                         break
#             return pubchem_cid, drugbank_id
#     except Exception as e:
#         print(f"PubChem lookup failed for {drug_name}: {e}")
#     return None, None

# def get_chembl_info(drug_name):
#     """Fetch ChEMBL ID from ChEMBL using REST API."""
#     try:
#         url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?pref_name__iexact={drug_name}"
#         response = requests.get(url, timeout=5)
#         if response.status_code == 200:
#             data = response.json()
#             molecules = data.get("molecules", [])
#             if molecules:
#                 return molecules[0].get("molecule_chembl_id")
#     except Exception as e:
#         print(f"ChEMBL lookup failed for {drug_name}: {e}")
#     return None

# def get_all_drug_ids(drug_name):
#     """Fetch all IDs (PubChem CID, DrugBank ID, and ChEMBL ID) for a given drug name."""
#     pubchem_cid, drugbank_id = get_pubchem_info(drug_name)
#     chembl_id = get_chembl_info(drug_name)
#     return {
#         "Drug Name": drug_name,
#         "PubChem CID": pubchem_cid,
#         "DrugBank ID": drugbank_id,
#         "ChEMBL ID": chembl_id
#     }

# def save_to_csv(drug_info, filename="drug_names_ids.csv"):
#     """Save drug information to a CSV file."""
#     file_exists = os.path.isfile(filename)
#     with open(filename, mode="a", newline="", encoding="utf-8") as file:
#         writer = csv.DictWriter(file, fieldnames=["Drug Name", "PubChem CID", "DrugBank ID", "ChEMBL ID"])
#         if not file_exists:
#             writer.writeheader()
#         writer.writerow(drug_info)

# def is_drug_in_csv(drug_name, filename="drug_names_ids.csv"):
#     """Check if a drug is already in the CSV file."""
#     if not os.path.isfile(filename):
#         return False
#     with open(filename, mode="r", encoding="utf-8") as file:
#         reader = csv.DictReader(file)
#         for row in reader:
#             if row["Drug Name"].lower() == drug_name.lower():
#                 return True
#     return False

# def load_similarity_matrix(filename="Drug_data/chem_similarity.csv"):
#     """Load the drug similarity matrix from a CSV file."""
#     return pd.read_csv(filename, index_col=0)

# def get_top_similar_drugs(drugbank_id, similarity_matrix, top_n=3):
#     """Get the top N similar drugs for a given DrugBank ID."""
#     if drugbank_id not in similarity_matrix.index:
#         return []
#     similarities = similarity_matrix.loc[drugbank_id]
#     top_drugs = similarities.sort_values(ascending=False).index[1:top_n + 1]
#     return list(top_drugs)

# def get_drug_name_from_drugbank_id(drugbank_id):
#     """Fetch the drug name from PubChem or ChEMBL using the DrugBank ID."""
#     try:
#         results = pcp.get_compounds(drugbank_id, 'name')
#         if results:
#             compound = results[0]
#             if compound.synonyms:
#                 return compound.synonyms[0]
#     except Exception as e:
#         print(f"PubChem lookup failed for {drugbank_id}: {e}")
#     try:
#         url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_chembl_id={drugbank_id}"
#         response = requests.get(url, timeout=5)
#         if response.status_code == 200:
#             data = response.json()
#             molecules = data.get("molecules", [])
#             if molecules:
#                 return molecules[0].get("pref_name")
#     except Exception as e:
#         print(f"ChEMBL lookup failed for {drugbank_id}: {e}")
#     return None

# def save_similar_drugs_to_csv(input_drug_name, similar_drugbank_ids, filename="SAME_DRUG.csv"):
#     """Save the input drug name and its similar drugs (resolved names) to a CSV file."""
#     file_exists = os.path.isfile(filename)
#     similar_drug_names = []
#     for drugbank_id in similar_drugbank_ids:
#         drug_name = get_drug_name_from_drugbank_id(drugbank_id)
#         if drug_name:
#             similar_drug_names.append(drug_name)
#         else:
#             print(f"❌ No name found for DrugBank ID {drugbank_id}. Skipping...")
#     with open(filename, mode="a", newline="", encoding="utf-8") as file:
#         fieldnames = ["Input Drug Name", "Similar Drug 1", "Similar Drug 2", "Similar Drug 3"]
#         writer = csv.DictWriter(file, fieldnames=fieldnames)
#         if not file_exists:
#             writer.writeheader()
#         writer.writerow({
#             "Input Drug Name": input_drug_name,
#             "Similar Drug 1": similar_drug_names[0] if len(similar_drug_names) > 0 else "",
#             "Similar Drug 2": similar_drug_names[1] if len(similar_drug_names) > 1 else "",
#             "Similar Drug 3": similar_drug_names[2] if len(similar_drug_names) > 2 else ""
#         })
#     return similar_drug_names



# def is_drug_in_same_drug_csv(drug_name, filename="SAME_DRUG.csv"):
#     """Check if a drug is already in SAME_DRUG.csv."""
#     if not os.path.isfile(filename):
#         return False
#     with open(filename, mode="r", encoding="utf-8") as file:
#         reader = csv.DictReader(file)
#         for row in reader:
#             if row["Input Drug Name"].lower() == drug_name.lower():
#                 return True
#     return False

# def get_similar_drugs_from_csv(drug_name, filename="SAME_DRUG.csv"):
#     """Get similar drugs for a given drug from SAME_DRUG.csv."""
#     similar_drugs = []
#     if not os.path.isfile(filename):
#         return similar_drugs
#     with open(filename, mode="r", encoding="utf-8") as file:
#         reader = csv.DictReader(file)
#         for row in reader:
#             if row["Input Drug Name"].lower() == drug_name.lower():
#                 similar_drugs = [row["Similar Drug 1"], row["Similar Drug 2"], row["Similar Drug 3"]]
#                 similar_drugs = [drug for drug in similar_drugs if drug]
#                 break
#     return similar_drugs



# def generate_combinations_from_same_drug(drug_names, filename="SAME_DRUG.csv"):
#     """Generate combinations of drugs based on similar drugs from SAME_DRUG.csv."""
#     same_drug_dict = {}
#     if not os.path.isfile(filename):
#         print(f"❌ {filename} not found. Please ensure the file exists.")
#         return []
#     with open(filename, mode="r", encoding="utf-8") as file:
#         reader = csv.DictReader(file)
#         for row in reader:
#             input_drug = row["Input Drug Name"]
#             similar_drugs = [row["Similar Drug 1"], row["Similar Drug 2"], row["Similar Drug 3"]]
#             similar_drugs = [drug for drug in similar_drugs if drug]
#             same_drug_dict[input_drug] = similar_drugs
#     groups = []
#     for drug in drug_names:
#         if drug in same_drug_dict:
#             group = [drug] + same_drug_dict[drug]
#             groups.append(group)
#         else:
#             print(f"❌ {drug} not found in SAME_DRUG.csv. Skipping...")
#     combinations = set()
#     for i in range(len(groups)):
#         for j in range(i + 1, len(groups)):
#             for drug1 in groups[i]:
#                 for drug2 in groups[j]:
#                     pair = tuple(sorted([drug1, drug2]))
#                     combinations.add(pair)
#     return sorted([f"{pair[0]},{pair[1]}" for pair in combinations])

# def fetch_drug_data(drug_name):
#     """Fetch drug interaction data from PubChem and save it to MongoDB."""
#     existing_document = collection.find_one({"drug_name": drug_name})
#     if existing_document:
#         print(f"✅ Data for {drug_name} already exists in MongoDB. Skipping...")
#         return
#     url = f"https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={{%22download%22:%22*%22,%22collection%22:%22drugbankddi%22,%22order%22:[%22cid2,asc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22pubchem_name_%5E{drug_name}%24_drugbankddi%22,%22where%22:{{%22ands%22:[{{%22name%22:%22%5E{drug_name}%24%22}}]}}}}"
#     try:
#         response = httpx.get(url)
#         response.raise_for_status()
#         file_name = f"{drug_name}_response.csv"
#         document = {
#             "drug_name": drug_name,
#             "file_name": file_name,
#             "content": response.text
#         }
#         collection.insert_one(document)
#         print(f"File saved to MongoDB as '{file_name}'")
#     except httpx.HTTPStatusError as e:
#         print(f"HTTP error occurred: {e}")
#     except Exception as e:
#         print(f"An unexpected error occurred: {e}")

# def search_interactions_for_drugs(drugs_to_check):
#     """Search for interactions between the given drugs using data from MongoDB."""
#     interactions = []
#     for drug_name in drugs_to_check:
#         fetch_drug_data(drug_name.strip())
#     for drug_name in drugs_to_check:
#         file_name = f"{drug_name.strip()}_response.csv"
#         document = collection.find_one({"file_name": file_name})
#         if document:
#             content = document["content"]
#             reader = csv.reader(content.splitlines())
#             for row in reader:
#                 drug1 = row[3]
#                 drug2 = row[5]
#                 interaction_description = row[6]
#                 if drug1 in drugs_to_check and drug2 in drugs_to_check:
#                     interaction_tuple = tuple(sorted([drug1, drug2]))
#                     # Avoid duplicate interactions
#                     if not any(interaction_tuple == existing_tuple for existing_tuple, _ in interactions):
#                         interactions.append((interaction_tuple, interaction_description))
#                         print(f"{drug1} - {drug2}: {interaction_description}")
#         else:
#             print(f"File for {drug_name.strip()} not found in MongoDB.")
#     if not interactions:
#         print("No interactions found between the drugs.")
#     return interactions

# # Flask Routes

# @form.route("/", methods=["GET", "POST"])
# def index():
#     if request.method == "POST":
#         input_drug_names = request.form.get("drug_names").strip().split()
#         for drug_name in input_drug_names:
#             if not is_drug_in_csv(drug_name):
#                 drug_info = get_all_drug_ids(drug_name)
#                 if drug_info["PubChem CID"] or drug_info["DrugBank ID"] or drug_info["ChEMBL ID"]:
#                     save_to_csv(drug_info)
#                     print(f"✅ Information for {drug_name} saved to drug_names_ids.csv.")
#                 else:
#                     print(f"❌ No information found for {drug_name}.")
#             else:
#                 print(f"✅ {drug_name} already exists in drug_names_ids.csv. Skipping...")
#         similar_drugs_info = {}
#         for drug_name in input_drug_names:
#             if is_drug_in_same_drug_csv(drug_name):
#                 print(f"✅ Similar drugs for {drug_name} already exist in SAME_DRUG.csv. Skipping lookup...")
#                 similar_drugs = get_similar_drugs_from_csv(drug_name)
#                 similar_drugs_info[drug_name] = similar_drugs
#             else:
#                 print(f"❌ Similar drugs for {drug_name} not found in SAME_DRUG.csv. Fetching...")
#                 similarity_matrix = load_similarity_matrix()
#                 drug_info = get_all_drug_ids(drug_name)
#                 drugbank_id = drug_info["DrugBank ID"]
#                 if drugbank_id:
#                     similar_drugbank_ids = get_top_similar_drugs(drugbank_id, similarity_matrix)
#                     similar_drug_names = save_similar_drugs_to_csv(drug_name, similar_drugbank_ids)
#                     similar_drugs_info[drug_name] = similar_drug_names
#                     print(f"✅ Similar drugs for {drug_name}: {similar_drug_names}")
#                 else:
#                     print(f"❌ No DrugBank ID found for {drug_name}.")
#         combinations = generate_combinations_from_same_drug(input_drug_names)
#         if combinations:
#             print("\n✅ Generated Combinations and Interactions:")
#             interactions = []
#             for combo in combinations:
#                 print(f"\n🔍 Checking interactions for: {combo}")
#                 drugs_to_check = [drug.strip() for drug in combo.split(",")]
#                 interaction_results = search_interactions_for_drugs(drugs_to_check)
#                 if interaction_results:
#                     interactions.append((combo, interaction_results))
#             return render_template("results.html", input_drug_names=input_drug_names, similar_drugs_info=similar_drugs_info, combinations=combinations, interactions=interactions)
#         else:
#             print("❌ No combinations generated.")
#             return render_template("results.html", input_drug_names=input_drug_names, similar_drugs_info=similar_drugs_info, combinations=[], interactions=[])
#     return render_template("form.html")



# from flask import Flask, render_template, request
# import httpx
# import csv
# import os
# import pandas as pd
# from pymongo import MongoClient
# import pubchempy as pcp
# import requests
# from flask import Blueprint

# # Initialize Flask Blueprint
# form = Blueprint('form', __name__, static_folder='static', template_folder='templates')

# # MongoDB connection
# connection_string = "mongodb://localhost:27017/"
# client = MongoClient(connection_string)
# db = client["Drug_Interaction"]
# collection = db["interaction_files"]

# print("✅ Connected to MongoDB.")

# # Helper Functions

# def get_pubchem_info(drug_name):
#     """Fetch PubChem CID and DrugBank ID from PubChem using the drug name."""
#     try:
#         results = pcp.get_compounds(drug_name, 'name')
#         if results:
#             compound = results[0]
#             pubchem_cid = compound.cid
#             drugbank_id = None
#             if compound.synonyms:
#                 for synonym in compound.synonyms:
#                     if synonym.startswith("DB"):
#                         drugbank_id = synonym
#                         break
#             return pubchem_cid, drugbank_id
#     except Exception as e:
#         print(f"PubChem lookup failed for {drug_name}: {e}")
#     return None, None

# def get_chembl_info(drug_name):
#     """Fetch ChEMBL ID from ChEMBL using REST API."""
#     try:
#         url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?pref_name__iexact={drug_name}"
#         response = requests.get(url, timeout=5)
#         if response.status_code == 200:
#             data = response.json()
#             molecules = data.get("molecules", [])
#             if molecules:
#                 return molecules[0].get("molecule_chembl_id")
#     except Exception as e:
#         print(f"ChEMBL lookup failed for {drug_name}: {e}")
#     return None

# def get_all_drug_ids(drug_name):
#     """Fetch all IDs (PubChem CID, DrugBank ID, and ChEMBL ID) for a given drug name."""
#     pubchem_cid, drugbank_id = get_pubchem_info(drug_name)
#     chembl_id = get_chembl_info(drug_name)
#     return {
#         "Drug Name": drug_name,
#         "PubChem CID": pubchem_cid,
#         "DrugBank ID": drugbank_id,
#         "ChEMBL ID": chembl_id
#     }

# def save_to_csv(drug_info, filename="drug_names_ids.csv"):
#     """Save drug information to a CSV file."""
#     file_exists = os.path.isfile(filename)
#     with open(filename, mode="a", newline="", encoding="utf-8") as file:
#         writer = csv.DictWriter(file, fieldnames=["Drug Name", "PubChem CID", "DrugBank ID", "ChEMBL ID"])
#         if not file_exists:
#             writer.writeheader()
#         writer.writerow(drug_info)

# def is_drug_in_csv(drug_name, filename="drug_names_ids.csv"):
#     """Check if a drug is already in the CSV file."""
#     if not os.path.isfile(filename):
#         return False
#     with open(filename, mode="r", encoding="utf-8") as file:
#         reader = csv.DictReader(file)
#         for row in reader:
#             if row["Drug Name"].lower() == drug_name.lower():
#                 return True
#     return False

# def load_similarity_matrix(filename="Drug_data/chem_similarity.csv"):
#     """Load the drug similarity matrix from a CSV file."""
#     return pd.read_csv(filename, index_col=0)

# def get_top_similar_drugs(drugbank_id, similarity_matrix, top_n=3):
#     """Get the top N similar drugs for a given DrugBank ID."""
#     if drugbank_id not in similarity_matrix.index:
#         return []
#     similarities = similarity_matrix.loc[drugbank_id]
#     top_drugs = similarities.sort_values(ascending=False).index[1:top_n + 1]
#     return list(top_drugs)

# def get_drug_name_from_drugbank_id(drugbank_id):
#     """Fetch the drug name from PubChem or ChEMBL using the DrugBank ID."""
#     try:
#         results = pcp.get_compounds(drugbank_id, 'name')
#         if results:
#             compound = results[0]
#             if compound.synonyms:
#                 return compound.synonyms[0]
#     except Exception as e:
#         print(f"PubChem lookup failed for {drugbank_id}: {e}")
#     try:
#         url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_chembl_id={drugbank_id}"
#         response = requests.get(url, timeout=5)
#         if response.status_code == 200:
#             data = response.json()
#             molecules = data.get("molecules", [])
#             if molecules:
#                 return molecules[0].get("pref_name")
#     except Exception as e:
#         print(f"ChEMBL lookup failed for {drugbank_id}: {e}")
#     return None

# def save_similar_drugs_to_csv(input_drug_name, similar_drugbank_ids, filename="SAME_DRUG.csv"):
#     """Save the input drug name and its similar drugs (resolved names) to a CSV file."""
#     file_exists = os.path.isfile(filename)
#     similar_drug_names = []
#     for drugbank_id in similar_drugbank_ids:
#         drug_name = get_drug_name_from_drugbank_id(drugbank_id)
#         if drug_name:
#             similar_drug_names.append(drug_name)
#         else:
#             print(f"❌ No name found for DrugBank ID {drugbank_id}. Skipping...")
#     with open(filename, mode="a", newline="", encoding="utf-8") as file:
#         fieldnames = ["Input Drug Name", "Similar Drug 1", "Similar Drug 2", "Similar Drug 3"]
#         writer = csv.DictWriter(file, fieldnames=fieldnames)
#         if not file_exists:
#             writer.writeheader()
#         writer.writerow({
#             "Input Drug Name": input_drug_name,
#             "Similar Drug 1": similar_drug_names[0] if len(similar_drug_names) > 0 else "",
#             "Similar Drug 2": similar_drug_names[1] if len(similar_drug_names) > 1 else "",
#             "Similar Drug 3": similar_drug_names[2] if len(similar_drug_names) > 2 else ""
#         })
#     return similar_drug_names

# def is_drug_in_same_drug_csv(drug_name, filename="SAME_DRUG.csv"):
#     """Check if a drug is already in SAME_DRUG.csv."""
#     if not os.path.isfile(filename):
#         return False
#     with open(filename, mode="r", encoding="utf-8") as file:
#         reader = csv.DictReader(file)
#         for row in reader:
#             if row["Input Drug Name"].lower() == drug_name.lower():
#                 return True
#     return False

# def get_similar_drugs_from_csv(drug_name, filename="SAME_DRUG.csv"):
#     """Get similar drugs for a given drug from SAME_DRUG.csv."""
#     similar_drugs = []
#     if not os.path.isfile(filename):
#         return similar_drugs
#     with open(filename, mode="r", encoding="utf-8") as file:
#         reader = csv.DictReader(file)
#         for row in reader:
#             if row["Input Drug Name"].lower() == drug_name.lower():
#                 similar_drugs = [row["Similar Drug 1"], row["Similar Drug 2"], row["Similar Drug 3"]]
#                 similar_drugs = [drug for drug in similar_drugs if drug]
#                 break
#     return similar_drugs

# def generate_combinations_from_same_drug(drug_names, filename="SAME_DRUG.csv"):
#     """Generate combinations of drugs based on similar drugs from SAME_DRUG.csv."""
#     same_drug_dict = {}
#     if not os.path.isfile(filename):
#         print(f"❌ {filename} not found. Please ensure the file exists.")
#         return []
#     with open(filename, mode="r", encoding="utf-8") as file:
#         reader = csv.DictReader(file)
#         for row in reader:
#             input_drug = row["Input Drug Name"]
#             similar_drugs = [row["Similar Drug 1"], row["Similar Drug 2"], row["Similar Drug 3"]]
#             similar_drugs = [drug for drug in similar_drugs if drug]
#             same_drug_dict[input_drug] = similar_drugs
#     groups = []
#     for drug in drug_names:
#         if drug in same_drug_dict:
#             group = [drug] + same_drug_dict[drug]
#             groups.append(group)
#         else:
#             print(f"❌ {drug} not found in SAME_DRUG.csv. Skipping...")
#     combinations = set()
#     for i in range(len(groups)):
#         for j in range(i + 1, len(groups)):
#             for drug1 in groups[i]:
#                 for drug2 in groups[j]:
#                     pair = tuple(sorted([drug1, drug2]))
#                     combinations.add(pair)
#     return sorted([f"{pair[0]},{pair[1]}" for pair in combinations])

# def fetch_drug_data(drug_name):
#     """Fetch drug interaction data from PubChem and save it to MongoDB."""
#     existing_document = collection.find_one({"drug_name": drug_name})
#     if existing_document:
#         print(f"✅ Data for {drug_name} already exists in MongoDB. Skipping...")
#         return
#     url = f"https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={{%22download%22:%22*%22,%22collection%22:%22drugbankddi%22,%22order%22:[%22cid2,asc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22pubchem_name_%5E{drug_name}%24_drugbankddi%22,%22where%22:{{%22ands%22:[{{%22name%22:%22%5E{drug_name}%24%22}}]}}}}"
#     try:
#         response = httpx.get(url)
#         response.raise_for_status()
#         file_name = f"{drug_name}_response.csv"
#         document = {
#             "drug_name": drug_name,
#             "file_name": file_name,
#             "content": response.text
#         }
#         collection.insert_one(document)
#         print(f"File saved to MongoDB as '{file_name}'")
#     except httpx.HTTPStatusError as e:
#         print(f"HTTP error occurred: {e}")
#     except Exception as e:
#         print(f"An unexpected error occurred: {e}")

# def search_interactions_for_drugs(drugs_to_check):
#     """Search for interactions between the given drugs using data from MongoDB."""
#     interactions = set()
#     for drug_name in drugs_to_check:
#         fetch_drug_data(drug_name.strip())
#     for drug_name in drugs_to_check:
#         file_name = f"{drug_name.strip()}_response.csv"
#         document = collection.find_one({"file_name": file_name})
#         if document:
#             content = document["content"]
#             reader = csv.reader(content.splitlines())
#             for row in reader:
#                 drug1 = row[3]
#                 drug2 = row[5]
#                 interaction = row[6]
#                 if drug1 in drugs_to_check and drug2 in drugs_to_check:
#                     interaction_tuple = tuple(sorted([drug1, drug2]))
#                     if interaction_tuple not in interactions:
#                         interactions.add(interaction_tuple)
#                         print(f"{drug1} - {drug2}: {interaction}")
#         else:
#             print(f"File for {drug_name.strip()} not found in MongoDB.")
#     if not interactions:
#         print("No interactions found between the drugs.")
#     return interactions

# # Flask Routes

# @form.route("/", methods=["GET", "POST"])
# def index():
#     if request.method == "POST":
#         input_drug_names = request.form.get("drug_names").strip().split()
#         for drug_name in input_drug_names:
#             if not is_drug_in_csv(drug_name):
#                 drug_info = get_all_drug_ids(drug_name)
#                 if drug_info["PubChem CID"] or drug_info["DrugBank ID"] or drug_info["ChEMBL ID"]:
#                     save_to_csv(drug_info)
#                     print(f"✅ Information for {drug_name} saved to drug_names_ids.csv.")
#                 else:
#                     print(f"❌ No information found for {drug_name}.")
#             else:
#                 print(f"✅ {drug_name} already exists in drug_names_ids.csv. Skipping...")
#         similar_drugs_info = {}
#         for drug_name in input_drug_names:
#             if is_drug_in_same_drug_csv(drug_name):
#                 print(f"✅ Similar drugs for {drug_name} already exist in SAME_DRUG.csv. Skipping lookup...")
#                 similar_drugs = get_similar_drugs_from_csv(drug_name)
#                 similar_drugs_info[drug_name] = similar_drugs
#             else:
#                 print(f"❌ Similar drugs for {drug_name} not found in SAME_DRUG.csv. Fetching...")
#                 similarity_matrix = load_similarity_matrix()
#                 drug_info = get_all_drug_ids(drug_name)
#                 drugbank_id = drug_info["DrugBank ID"]
#                 if drugbank_id:
#                     similar_drugbank_ids = get_top_similar_drugs(drugbank_id, similarity_matrix)
#                     similar_drug_names = save_similar_drugs_to_csv(drug_name, similar_drugbank_ids)
#                     similar_drugs_info[drug_name] = similar_drug_names
#                     print(f"✅ Similar drugs for {drug_name}: {similar_drug_names}")
#                 else:
#                     print(f"❌ No DrugBank ID found for {drug_name}.")
#         combinations = generate_combinations_from_same_drug(input_drug_names)
#         if combinations:
#             print("\n✅ Generated Combinations and Interactions:")
#             interactions = []
#             for combo in combinations:
#                 print(f"\n🔍 Checking interactions for: {combo}")
#                 drugs_to_check = [drug.strip() for drug in combo.split(",")]
#                 interaction_results = search_interactions_for_drugs(drugs_to_check)
#                 if interaction_results:
#                     interactions.append((combo, interaction_results))
#             return render_template("results.html", input_drug_names=input_drug_names, similar_drugs_info=similar_drugs_info, combinations=combinations, interactions=interactions)
#         else:
#             print("❌ No combinations generated.")
#             return render_template("results.html", input_drug_names=input_drug_names, similar_drugs_info=similar_drugs_info, combinations=[], interactions=[])
#     return render_template("form.html")


#---------------------------------------------------------------------------------------------------------------------------------------------------