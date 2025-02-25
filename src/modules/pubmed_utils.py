import json
from Bio import Entrez

from modules.spaCy_utils import process_text


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
    """Fetches article details including title, abstract, authors, keywords, journal, and DOI."""
    if not id_list:
        return []

    ids = ",".join(id_list)
    handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="xml")
    records = Entrez.read(handle)
    handle.close()

    results = []
    for article in records["PubmedArticle"]:
        medline = article["MedlineCitation"]
        article_data = medline["Article"]

        title = article_data.get("ArticleTitle", "No Title Available")
        abstract_data = article_data.get("Abstract", {}).get("AbstractText", ["No Abstract"])
        abstract = " ".join(abstract_data) if isinstance(abstract_data, list) else str(abstract_data)
        keywords = medline.get("KeywordList", [])
        keywords = [kw for sublist in keywords for kw in sublist] if keywords else ["No Keywords"]
        
        # Extracting authors
        authors = []
        if "AuthorList" in article_data:
            for author in article_data["AuthorList"]:
                if "LastName" in author and "ForeName" in author:
                    authors.append(f"{author['ForeName']} {author['LastName']}")

        # Extracting journal name and DOI
        journal = article_data.get("Journal", {}).get("Title", "No Journal Info")
        doi = "No DOI"
        if "ELocationID" in article_data:
            for eloc in article_data["ELocationID"]:
                if eloc.attributes.get("EIdType") == "doi":
                    doi = eloc.lower()

        # spaCy processing
        spacy_results = process_text(abstract)

        results.append({
            "title": title, 
            "abstract": abstract, 
            "keywords": keywords,
            "authors": authors,
            "journal": journal,
            "doi": doi,
            "spacy_entities": spacy_results["entities"],
            "spacy_matched_terms": spacy_results["matched_terms"]
        })

    return results


def save_results_to_json(articles, filename="pubmed_results.json"):
    """Saves the results in a JSON file."""
    with open(filename, mode='w', encoding='utf-8') as file:
        json.dump(articles, file, ensure_ascii=False, indent=4)
    print(f"Results saved in {filename}")
