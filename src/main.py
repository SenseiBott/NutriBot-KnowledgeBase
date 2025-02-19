import os
from dotenv import load_dotenv
from pubmed_utils import configure_entrez, save_results_to_json, search_pubmed, fetch_article_details, print_article_results

load_dotenv()

def main():
    """Main function to search and display PubMed articles."""
    email = os.getenv("EMAIL")
    api_key = os.getenv("API_KEY_PUBMED")

    if not email or not api_key:
        raise ValueError("Missing EMAIL or API_KEY_PUBMED in environment variables.")

    configure_entrez(email, api_key)

    query = '("dietary supplements" OR "nutritional supplements") AND ("disease prevention" OR "health benefits")'
    max_articles = 10

    print(f"Searching for articles with query: {query}")
    article_ids = search_pubmed(query, max_results=max_articles)
    print(f"Found article IDs: {article_ids}")

    articles = fetch_article_details(article_ids)
    print_article_results(articles)
    save_results_to_json(articles, filename="articles.json")

if __name__ == "__main__":
    main()
