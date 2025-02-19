import json
from Bio import Entrez

def configure_entrez(email, api_key):
    """Configure Entrez credentials."""
    Entrez.email = email
    Entrez.api_key = api_key

def search_pubmed(query, max_results=10):
    """Search for articles on PubMed and return a list of IDs."""
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def fetch_article_details(id_list):
    """Fetch abstracts and keywords of articles by ID."""
    if not id_list:
        return []
    
    ids = ",".join(id_list)
    handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="xml")
    records = Entrez.read(handle)
    handle.close()

    results = []
    for article in records["PubmedArticle"]:
        title = article["MedlineCitation"]["Article"]["ArticleTitle"]
        abstract = article["MedlineCitation"]["Article"].get("Abstract", {}).get("AbstractText", ["N/A"])[0]
        keywords = article["MedlineCitation"].get("KeywordList", [])
        keywords = [kw for sublist in keywords for kw in sublist] if keywords else ["N/A"]
        
        results.append({"title": title, "abstract": abstract, "keywords": keywords})
    
    return results

def print_article_results(articles):
    """Displays the results in a formatted way."""
    if not articles:
        print("No articles found.")
        return

    for idx, article in enumerate(articles):
        print("\n===================================")
        print(f"Title {idx+1}: {article['title']}")
        print(f"Abstract: {article['abstract']}")
        print(f"Keywords: {', '.join(article['keywords'])}")

def save_results_to_json(articles, filename="articles.json"):
    """Save the article results to a JSON file."""
    with open(filename, mode='w', encoding='utf-8') as file:
        json.dump(articles, file, ensure_ascii=False, indent=4)

    print(f"Articles saved in {filename}")
