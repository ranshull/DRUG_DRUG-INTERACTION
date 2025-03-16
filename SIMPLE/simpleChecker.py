from flask import Flask, render_template, request, redirect, url_for
from flask import Flask, render_template, request, jsonify
import httpx
import csv
import os
from flask import Blueprint, render_template

simpleChecker = Blueprint('simpleChecker', __name__, static_folder='static', template_folder='templates')


DRUG_DATA_FOLDER = "Drug_data"

def fetch_drug_data(drug_name):
    """Fetch interaction data for a drug from PubChem API."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={{%22download%22:%22*%22,%22collection%22:%22drugbankddi%22,%22order%22:[%22cid2,asc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22pubchem_name_%5E{drug_name}%24_drugbankddi%22,%22where%22:{{%22ands%22:[{{%22name%22:%22%5E{drug_name}%24%22}}]}}}}"

    try:
        response = httpx.get(url)
        response.raise_for_status()
        # filename = f"{drug_name}_response.csv"
        filename = os.path.join(DRUG_DATA_FOLDER, f"{drug_name}_response.csv")

        # Save API response as a CSV file
        with open(filename, "wb") as file:
            file.write(response.content)
        return filename

    except httpx.HTTPStatusError as e:
        return f"HTTP error: {e}"
    except Exception as e:
        return f"Error: {e}"

def search_interactions_for_drugs(drugs_to_check):
    """Check interactions among multiple drugs."""
    interactions = []

    # Fetch data for each drug
    for drug_name in drugs_to_check:
        fetch_drug_data(drug_name)

    # Read saved files and search for interactions
    for drug_name in drugs_to_check:
        # filename = f"{drug_name}_response.csv"
        filename = os.path.join(DRUG_DATA_FOLDER, f"{drug_name}_response.csv")
        if os.path.exists(filename):
            with open(filename, newline='', encoding='utf-8') as file:
                reader = csv.reader(file)
                for row in reader:
                    if len(row) < 7:
                        continue
                    drug1, drug2, interaction = row[3], row[5], row[6]
                    if drug1 in drugs_to_check and drug2 in drugs_to_check:
                        interactions.append(f"{drug1} - {drug2}: {interaction}")
    
    return interactions if interactions else ["No interactions found."]

@simpleChecker .route("/")
def index():
    return render_template("simpleChecker.html")

@simpleChecker .route("/check_interaction", methods=["POST"])
def check_interaction():
    """Handle AJAX requests for drug interaction checks."""
    drugs = request.json.get("drugs")
    if not drugs:
        return jsonify({"error": "No drugs provided."}), 400

    drugs_list = [drug.strip() for drug in drugs.split()]
    results = search_interactions_for_drugs(drugs_list)
    return jsonify({"interactions": results})

