import requests
import json
from datetime import datetime

def search_semantic_scholar(query, max_results=20):
    """Search for articles on Semantic Scholar."""
    url = "https://api.semanticscholar.org/graph/v1/paper/search"
    params = {
        "query": query,
        "limit": max_results,
        "fields": "title,abstract,authors,year,citationCount,referenceCount"
    }
    
    response = requests.get(url, params=params)
    if response.status_code == 200:
        return response.json()
    else:
        print(f"Error: {response.status_code}")
        return None

def save_results_to_json(articles, filename="semanticscholar_results.json"):
    """Saves the results to a JSON file."""
    with open(filename, "w", encoding="utf-8") as file:
        json.dump(articles, file, ensure_ascii=False, indent=4)
    print(f"Results saved in {filename}")