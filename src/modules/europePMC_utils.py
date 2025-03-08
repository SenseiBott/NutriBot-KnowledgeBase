import requests
from tqdm import tqdm

from modules.mongoDB_utils import configure_mongoDB_connection, save_to_mongo

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

def search_europe_pmc(query, num_results):
    """Fetch articles from the Europe PMC API and save them to MongoDB."""
    collection = configure_mongoDB_connection()
    papers = fetch_papers(query, num_results)
    save_to_mongo(papers, collection, "Europe PMC")
    return papers
