from pymongo import MongoClient
from tqdm import tqdm
from modules.spaCy_utils import process_text

def configure_mongoDB_connection():
    """Configure MongoDB connection."""
    client = MongoClient("mongodb://localhost:27017/")  
    db = client["supplement_agent"]  
    collection = db["papers"] 
    return collection


def save_to_mongo(papers, source):
    """Save articles to MongoDB."""
    if not papers:
        print("No articles to save.")
        return
    else:
        collection = configure_mongoDB_connection()
        
        if source == "PubMed":
            for paper in tqdm(papers, desc="Saving PubMed articles to MongoDB"):
                year = paper.get("year", 0)
                if year == "No Year Available":
                    year = 0

                doc = {
                    "title": paper.get("title", ""),
                    "authors": paper.get("authors", []),
                    "year": int(paper.get(year, 0)),
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
        
        elif source == "Europe PMC":
            for paper in tqdm(papers, desc="Saving EuropePMC articles to MongoDB"):
                
                # Process Text
                if paper.get("abstractText", ""):
                    abstract = paper.get("abstractText", "")
                    spacy_results = process_text(abstract)
                else:
                    abstract = ""

                # Transform authors into a list of strings
                authors = paper.get("authorList", {}).get("author", [])
                authors = [f"{author.get('firstName', '')} {author.get('lastName', '')}" for author in authors]

                doc = {
                    "title": paper.get("title", ""),
                    "authors": authors,
                    "year": int(paper.get("pubYear", 0) or 0),
                    "source": "Europe PMC",
                    "abstract": paper.get("abstractText", ""),
                    "keywords": paper.get("keywordList", {}).get("keyword", []),
                    "doi": paper.get("doi", ""),
                    "last_updated": paper.get("firstPublicationDate", ""),
                    "spacy_entities": spacy_results.get("entities", []),
                    "spacy_matched_terms": spacy_results.get("matched_terms", [])
                }
                collection.insert_one(doc)

        elif source == "Semantic Scholar":
            for paper in tqdm(papers, desc="Saving Semantic Scholar articles to MongoDB"):
                
                # Process Text
                if paper.get("abstract", ""):
                    abstract = paper.get("abstract", "")
                    spacy_results = process_text(abstract)
                else:
                    abstract = ""

                doc = {
                    "title": paper.get("title", ""),
                    "authors": [author.get("name", "") for author in paper.get("authors", [])],
                    "year": int(paper.get("year", 0) or 0),
                    "source": "Semantic Scholar", 
                    "abstract": paper.get("abstract", ""),
                    "keywords": [],  
                    "doi": paper.get("externalIds", {}).get("DOI", ""), 
                    "journal": paper.get("journal", {}).get("name", ""),
                    "last_updated": "",
                    "spacy_entities": spacy_results.get("entities", []),
                    "spacy_matched_terms": spacy_results.get("matched_terms", [])
                }
                collection.insert_one(doc)
        
        elif source == "GoogleScholar":
            for paper in tqdm(papers, desc="Saving Google Scholar articles to MongoDB"):

                # Process Text
                if paper.get("abstract", ""):
                    abstract = paper.get("abstract", "")
                    spacy_results = process_text(abstract)
                else:
                    abstract = ""

                authors = paper.get("authors", "No Authors")
                if isinstance(authors, list):
                    authors = ", ".join(authors)  # Combine authors into a single string

                doc = {
                    "title": paper.get("title", ""),
                    "authors": authors,
                    "year": int(paper.get("year", 0) or 0),
                    "source": "Google Scholar", 
                    "abstract": paper.get("abstract", ""),
                    "keywords": paper.get("keywords", []),  
                    "doi": paper.get("doi", ""),
                    "journal": paper.get("journal", ""),
                    "last_updated": "",
                    "spacy_entities": spacy_results.get("entities", []),
                    "spacy_matched_terms": spacy_results.get("matched_terms", [])
                }
                collection.insert_one(doc)

    print("All articles have been successfully saved to MongoDB!")
