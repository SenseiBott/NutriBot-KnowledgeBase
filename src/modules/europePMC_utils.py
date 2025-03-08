import requests
import json
from pymongo import MongoClient
from tqdm import tqdm

from modules.mongoDB_utils import configure_mongoDB_connection

def fetch_papers(query, num_results):
    """Fetch articles from the Europe PMC API."""
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    params = {
        "query": query,
        "resultType": "core",
        "pageSize": num_results,
        "format": "json"
    }
    
    response = requests.get(url, params=params)
    if response.status_code == 200:
        return response.json().get("resultList", {}).get("result", [])
    else:
        print(f"Error fetching data: {response.status_code}")
        return []


def save_to_mongo(papers, collection):
    """Save articles to MongoDB."""
    if not papers:
        print("No articles to save.")
        return
    
    for paper in tqdm(papers, desc="Saving articles to MongoDB"):
        doc = {
            "title": paper.get("title", ""),
            "authors": paper.get("authorString", ""),
            "year": int(paper.get("pubYear", 0)),
            "source": "Europe PMC",
            "abstract": paper.get("abstractText", ""),
            "keywords": paper.get("keywordList", {}).get("keyword", []),
            "doi": paper.get("doi", ""),
            "paper_id": paper.get("id", ""),
            "last_updated": paper.get("firstPublicationDate", ""),
        }
        collection.insert_one(doc)
    
    print("All articles have been successfully saved to MongoDB!")

def search_europe_pmc(query, num_results):
    """Fetch articles from the Europe PMC API and save them to MongoDB."""
    collection = configure_mongoDB_connection()
    papers = fetch_papers(query, num_results)
    save_to_mongo(papers, collection)
    return papers
