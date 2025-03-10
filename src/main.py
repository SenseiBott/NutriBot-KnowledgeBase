from dotenv import load_dotenv
from modules.europePMC_utils import search_europe_pmc
from modules.pubmed_utils import search_pubmed
from modules.semanticscholar_utils import search_semanticscholar
from modules.wikipedia_utils import search_wikipedia 

def search_and_print(source_name, search_function, query, max_articles=None, year_range=None):
    """Generic function to search for articles from a given source."""
    print(f"🔎 Searching for articles on {source_name}...")
    
    if source_name == "Wikipedia":
        return search_function(query)
    elif year_range:
        return search_function(query, max_articles, year_range)
    else:
        return search_function(query, max_articles)


def main():
    """Main function to fetch and save articles from multiple sources."""
    
    # Load environment variables from the .env file
    load_dotenv()
    
    # Define the search query
    query = '("supplement AND disease prevention")'
    max_articles = 1
    year_range = (2020, 2025) # 
    
    # Define sources and their corresponding functions
    sources = {
        "PubMed": search_pubmed,
        "Europe PMC": search_europe_pmc,
        "Semantic Scholar": search_semanticscholar,
        "Wikipedia": search_wikipedia  # Nova fonte
    }

    
    results = {}
    for source, func in sources.items():
        if source == "Wikipedia":
            results[source] = search_and_print(source, func, query)
        else:
            results[source] = search_and_print(source, func, query, max_articles, year_range)

    if "Wikipedia" in results and results["Wikipedia"]:
        wiki_result = results["Wikipedia"]
        print(f"\n📖 Wikipedia: {wiki_result['title']}")
        print(f"🔗 URL: {wiki_result['url']}")
        print(f"📝 Resumo: {wiki_result['summary']}")
    else:
        print("\n❌ Nenhum resultado encontrado na Wikipedia.")

if __name__ == "__main__":
    main()