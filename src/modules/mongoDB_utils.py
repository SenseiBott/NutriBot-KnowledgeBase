# Configure MongoDB connection
from pymongo import MongoClient
from tqdm import tqdm

def configure_mongoDB_connection():
    """Configure MongoDB connection."""
    client = MongoClient("mongodb://localhost:27017/")  
    db = client["supplement_agent"]  
    collection = db["papers"] 
    return collection


def save_to_mongo(papers,collection, source):
    
    """Save articles to MongoDB."""
    if not papers:
        print("No articles to save.")
        return
    
    if source == "PubMed":
        for paper in tqdm(papers, desc="Saving PubMed articles to MongoDB"):
            doc = {
                "title": paper.get("title", ""),
                "authors": paper.get("authors", []),
                "year": int(paper.get("year", 0)),
                "source": "PubMed",
                "abstract": paper.get("abstract", ""),
                "keywords": paper.get("keywords", []),
                "doi": paper.get("doi", ""),
                "journal": paper.get("journal", ""),
                "last_updated": paper.get("last_updated", ""),
                "spacy_entities": paper.get("spacy_entities", []),
                "spacy_matched_terms": paper.get("spacy_matched_terms", [])
            }
            collection.insert_one(doc)
    else:
        if source == "Europe PMC":
            for paper in tqdm(papers, desc="Saving EuropePMC articles to MongoDB"):
                doc = {
                    "title": paper.get("title", ""),
                    "authors": paper.get("authorString", ""),
                    "year": int(paper.get("pubYear", 0)),
                    "source": "Europe PMC",
                    "abstract": paper.get("abstractText", ""),
                    "keywords": paper.get("keywordList", {}).get("keyword", []),
                    "doi": paper.get("doi", ""),
                    "last_updated": paper.get("firstPublicationDate", ""),
                }
                collection.insert_one(doc)

    print("All articles have been successfully saved to MongoDB!")
    