import requests
import json

def search_europe_pmc(query, max_results=20):
    """Search for articles in Europe PMC."""
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    params = {
        "query": query,
        "resultType": "core",
        "pageSize": max_results,
        "format": "json" 
    }
    
    response = requests.get(url, params=params)
    if response.status_code == 200:
        return response.json()
    else:
        print(f"Error: {response.status_code}")
        return None

def save_results_to_json(articles, filename="europe_pmc_results.json"):
    """Saves the results to a JSON file."""
    with open(filename, "w", encoding="utf-8") as file:
        json.dump(articles, file, ensure_ascii=False, indent=4)
    print(f"Results saved in {filename}")
