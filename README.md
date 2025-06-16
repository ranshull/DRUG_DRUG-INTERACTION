# MedSafe — Drug Interaction Finder

**Author:** Anshul
**Field:** Computer Science / Data Science

## 🔬 What is MedSafe?

**MedSafe** is a smart web-based tool designed to detect potentially harmful drug interactions using chemical structural similarity. It provides quick, reliable results that are accessible to researchers, students, healthcare professionals, and patients.

🔗 **Live Demo:** [https://drug-drug-interaction.onrender.com/](https://drug-drug-interaction.onrender.com/)

---

## 🚀 Features

* Perform interaction analysis using SMILES-based structural similarity
* Built-in checker and analyzer for drug combinations
* Uses molecular fingerprints and Jaccard similarity
* Integrated with DrugBank and PubChem datasets
* Fully deployed and accessible via the web (hosted on Render)
* Built with Python, Flask, and MongoDB

**Note:** Machine learning models were tested in notebooks but are not used in the deployed version. The live system uses similarity-based logic (no ML in production).

---

## 📈 Goals and Future Scope

### 🎯 Current Goals

* Help users avoid risky drug interactions
* Deliver fast and reliable interaction results
* Maintain an intuitive interface for both professionals and learners
* Ensure scalability and performance

### 🔮 Future Enhancements

* Expand beyond structural similarity to include:

  * Biological pathways
  * Enzyme and target interactions
  * Pharmacodynamic/pharmacokinetic data

---

## 🧪 How to Use

### 1. Checker Button  

Input drugs separated by spaces  
Example:  
Atenolol Phenytoin Xanthine  
Omeprazole Clopidogrel  
Warfarin Phenytoin  

### 2. Analyze Interaction Button

Input drugs separated by commas for a detailed report  
Example:  
Atenolol, Phenytoin, Xanthine  
Omeprazole, Clopidogrel  
Warfarin, Phenytoin  

---

## 🛠️ Tech Stack

* **Backend:** Python, Flask
* **Database:** MongoDB
* **Frontend:** HTML/CSS/JavaScript (Flask-based)
* **Deployment:** Render
* **Data Sources:** DrugBank, PubChem


