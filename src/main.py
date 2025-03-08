from dotenv import load_dotenv
from modules.europePMC_utils import search_europe_pmc
from modules.pubmed_utils import search_pubmed
from modules.semanticscholar_utils import search_semanticscholar

def search_and_print(source_name, search_function, query, max_articles):
    """Generic function to search for articles from a given source."""
    print(f"🔎 Searching for articles on {source_name}...")
    return search_function(query, max_articles)

def main():
    """Main function to fetch and save articles from multiple sources."""
    
    # Load environment variables from the .env file
    load_dotenv()
    
    # Define the search query
    query = '("supplement AND disease prevention")'
    max_articles = 1

    # Define sources and their corresponding functions
    sources = {
        "PubMed": search_pubmed,
        "Europe PMC": search_europe_pmc,
        "Semantic Scholar": search_semanticscholar
    }

    # Perform searches for each source
    results = {source: search_and_print(source, func, query, max_articles) for source, func in sources.items()}

if __name__ == "__main__":
    main()
