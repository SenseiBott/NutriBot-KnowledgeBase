import sys
import os
from dotenv import load_dotenv

# Add the 'modules' directory to the import path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'modules')))

from pubmed_utils import configure_entrez, search_pubmed, fetch_article_details, save_results_to_json
from europePMC_utils import search_europe_pmc, save_results_to_json as save_europe_results
from datetime import datetime

# Load environment variables from the .env file
load_dotenv()

def main():
    """Main function to fetch and save articles from PubMed and Europe PMC."""
    
    # PubMed configuration
    email = os.getenv("EMAIL")
    api_key = os.getenv("API_KEY_PUBMED")

    if not email or not api_key:
        raise ValueError("Missing EMAIL or API_KEY_PUBMED in environment variables.")
    
    configure_entrez(email, api_key)

    # Define the search query
    query = '("dietary supplements"[MeSH] OR "nutritional supplements") AND ("disease prevention"[MeSH] OR "health benefits")'
    max_articles = 20

    # Search for articles on PubMed
    print("🔎 Searching for articles on PubMed...")
    article_ids = search_pubmed(query, max_results=max_articles)
    articles = fetch_article_details(article_ids)
    pubmed_filename = f"data/articles_PubMed_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.json"
    save_results_to_json(articles, filename=pubmed_filename)

    # Search for articles on Europe PMC
    print("🔎 Searching for articles on Europe PMC...")
    europe_pmc_articles = search_europe_pmc(query)
    if europe_pmc_articles:
        europe_filename = f"data/articles_EuropePMC_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.json"
        save_europe_results(europe_pmc_articles, filename=europe_filename)

if __name__ == "__main__":
    main()
