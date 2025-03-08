from dotenv import load_dotenv
from modules.europePMC_utils import search_europe_pmc
from modules.pubmed_utils import search_pubmed

def main():
    """Main function to fetch and save articles from PubMed and Europe PMC."""
    
    # Load environment variables from the .env file
    load_dotenv()
    
    # Define the search query
    query = '("supplement AND disease prevention")'
    max_articles = 20

    # Search for articles on PubMed
    print("🔎 Searching for articles on PubMed...")
    pubmed_articles = search_pubmed(query, max_articles)

    # Search for articles on Europe PMC
    print("🔎 Searching for articles on Europe PMC...")
    europe_pmc_articles = search_europe_pmc(query, max_articles)
    

if __name__ == "__main__":
    main()
