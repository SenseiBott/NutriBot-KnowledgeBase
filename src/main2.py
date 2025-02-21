import sys
import os
from dotenv import load_dotenv
from datetime import datetime

# Add the 'modules' directory to the import path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'modules')))

from semanticscholar_utils import search_semantic_scholar, save_results_to_json as save_semantic_results

# Load environment variables from the .env file
load_dotenv()

def main():
    """Main function to fetch and save articles from Semantic Scholar."""
    
    # Define the search query
    query = '("dietary supplements" OR "nutritional supplements") AND ("disease prevention" OR "health benefits")'
    max_articles = 20

    # Search for articles on Semantic Scholar
    print("🔎 Searching for articles on Semantic Scholar...")
    semantic_scholar_articles = search_semantic_scholar(query, max_results=max_articles)
    if semantic_scholar_articles:
        semantic_filename = f"data/articles_SemanticScholar_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.json"
        save_semantic_results(semantic_scholar_articles, filename=semantic_filename)

if __name__ == "__main__":
    main()