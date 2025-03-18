from modules.europePMC_utils import search_europe_pmc
from modules.googleScholar_utils import search_google_scholar
from modules.pubmed_utils import search_pubmed
from modules.semanticscholar_utils import search_semanticscholar
from modules.wikipedia_utils import search_wikipedia 

def search_and_print(source, func, query, max_articles=1, year_range=(2020, 2025)):
    """Performs the search and prints the results."""
    print(f"🔎 Searching on {source}...")
    results = func(query) if source == "Wikipedia" else func(query, max_articles, year_range)
    
    if source == "Wikipedia" and results:
        print(f"\n📖 Wikipedia: {results['title']}")
        print(f"🔗 URL: {results['url']}")
        print(f"📝 Summary: {results['summary']}")
    
    return results

def main():
    """Runs searches on multiple sources."""
    query = '("dietary supplements" AND "disease prevention" AND "pharmaceuticals")'
    
    sources = {
        "PubMed": search_pubmed,
        "Europe PMC": search_europe_pmc,
        "Semantic Scholar": search_semanticscholar,
        "Wikipedia": search_wikipedia,
        #"Google Scholar": search_google_scholar
    }

    results = {src: search_and_print(src, func, query) for src, func in sources.items()}

    if not results.get("Wikipedia"):
        print("\n❌ No results found on Wikipedia.")

if __name__ == "__main__":
    main()
