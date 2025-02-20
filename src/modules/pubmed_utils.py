import json
from Bio import Entrez

def configure_entrez(email, api_key):
    """Configures Entrez credentials (PubMed)."""
    Entrez.email = email
    Entrez.api_key = api_key

def search_pubmed(query, max_results=10):
    """Searches for articles on PubMed and returns a list of IDs."""
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
        abstract = article["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", ["No Abstract"])[0]
        keywords = article["MedlineCitation"].get("KeywordList", [])
        keywords = [kw for sublist in keywords for kw in sublist] if keywords else ["No Keywords"]

        results.append({"title": title, "abstract": abstract, "keywords": keywords})
    
    return results

def save_results_to_json(articles, filename="pubmed_results.json"):
    """Saves the results to a JSON file."""
    with open(filename, mode='w', encoding='utf-8') as file:
        json.dump(articles, file, ensure_ascii=False, indent=4)
    print(f"Results saved in {filename}")
