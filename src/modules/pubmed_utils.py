import json
from Bio import Entrez

from spaCy_utils import process_text


def configure_entrez(email, api_key):
    """Configures Entrez credentials."""
    Entrez.email = email
    Entrez.api_key = api_key

def search_pubmed(query, max_results=10):
    """Searches PubMed articles and returns a list of IDs."""
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def fetch_article_details(id_list):
    """Fetches abstracts and keywords of articles by ID."""
    if not id_list:
        return []
    
    ids = ",".join(id_list)
    handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="xml")
    records = Entrez.read(handle)
    handle.close()

    results = []
    for article in records["PubmedArticle"]:
        title = article["MedlineCitation"]["Article"].get("ArticleTitle", "No Title Available")
        
        # Converts abstract to string
        abstract_data = article["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", ["No Abstract"])
        abstract = " ".join(abstract_data) if isinstance(abstract_data, list) else str(abstract_data)

        keywords = article["MedlineCitation"].get("KeywordList", [])
        keywords = [kw for sublist in keywords for kw in sublist] if keywords else ["No Keywords"]

        # Processar com spaCy
        spacy_results = process_text(abstract)

        results.append({
            "title": title, 
            "abstract": abstract, 
            "keywords": keywords,
            "spacy_entities": spacy_results["entities"],
            "spacy_matched_terms": spacy_results["matched_terms"]
        })
    
    return results

def save_results_to_json(articles, filename="pubmed_results.json"):
    """Saves the results in a JSON file."""
    with open(filename, mode='w', encoding='utf-8') as file:
        json.dump(articles, file, ensure_ascii=False, indent=4)
    print(f"Results saved in {filename}")
