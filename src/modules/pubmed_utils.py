import json
import os
from Bio import Entrez
from pymongo import MongoClient
from tqdm import tqdm
from modules.mongoDB_utils import configure_mongoDB_connection, save_to_mongo
from modules.spaCy_utils import process_text

def configure_entrez(email, api_key):
    """Configure the Entrez module with the provided email and API key."""
    Entrez.email = email
    Entrez.api_key = api_key

def configure_pubmed(email, api_key):
    # PubMed configuration
    email = os.getenv("EMAIL")
    api_key = os.getenv("API_KEY_PUBMED")

    if not email or not api_key:
        raise ValueError("Missing EMAIL or API_KEY_PUBMED in environment variables.")
    
    configure_entrez(email, api_key)

def fetch_papers(query, num_results=20):
    """Fetch articles from the PubMed API and return results in JSON format."""
    handle = Entrez.esearch(db="pubmed", term=query, retmax=num_results)
    record = Entrez.read(handle)
    handle.close()
    
    id_list = record.get("IdList", [])
    if not id_list:
        return json.dumps({"results": []}, indent=4)

    articles = fetch_article_details(id_list)
    return json.dumps({"results": articles}, ensure_ascii=False, indent=4)


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


def search_pubmed(query, num_results):
    """Fetch articles from the PubMed API and save them to MongoDB."""
    
    collection = configure_mongoDB_connection()
    articles_json = fetch_papers(query, num_results)
    
    # Parse JSON string into Python dictionary
    articles_dict = json.loads(articles_json)

    # Extract the list of articles
    articles = articles_dict.get("results", [])

    save_to_mongo(articles, collection, "PubMed")
    return articles
