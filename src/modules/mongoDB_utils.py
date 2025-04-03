import os
import uuid
from pymongo import MongoClient
from tqdm import tqdm
from modules.spaCy_utils import process_text
import os
from pinecone import Pinecone, ServerlessSpec

def configure_mongoDB_connection():
    """Configure MongoDB connection."""
    client = MongoClient("mongodb://localhost:27017/")  
    db = client["supplement_agent"]  
    collection = db["papers"] 
    return collection


def configure_pinecone_connection():
    pinecone_api_key = os.getenv("PINECONE_API_KEY")
    pc = Pinecone(api_key=pinecone_api_key)
    index_name = "papers"
    if index_name not in pc.list_indexes().names():
        pc.create_index(
            name=index_name,
            dimension=384,  # Match BGE-small-en
            metric='cosine',
            spec=ServerlessSpec(cloud='aws', region='us-east-1')
        )
    return pc.Index(index_name)

# Gerar um ID único
def generate_unique_id():
    """Generate a unique ID for each paper."""
    return str(uuid.uuid4())


def save_to_mongo_and_pinecone(papers, source):
    """Save articles to MongoDB."""
    if not papers:
        print("No articles to save.")
        return
    else:
        collection = configure_mongoDB_connection()
        index = configure_pinecone_connection()
        
        if source == "PubMed":
            for paper in tqdm(papers, desc="Saving PubMed articles to MongoDB"):
                year = paper.get("year", 0)
                if year == "No Year Available":
                    year = 0

                # Processar o abstract, se disponível
                abstract = paper.get("abstract", "")
                spacy_results = process_text(abstract) if abstract else {
                    "entities": [], "matched_terms": {}, "chunks": [], "embeddings": []
                }

                paper_id = generate_unique_id()  # Gerar um ID único para o artigo  
                
                doc = {
                    "paper_id": paper_id,
                    "title": paper.get("title", ""),
                    "authors": paper.get("authors", []),
                    "year": int(year),
                    "source": "PubMed",
                    "abstract": abstract,
                    "keywords": paper.get("keywords", []),
                    "doi": paper.get("doi", ""),
                    "journal": paper.get("journal", ""),
                    "last_updated": paper.get("last_updated", ""),
                    "spacy_entities": spacy_results["entities"],
                    "spacy_matched_terms": spacy_results["matched_terms"],
                    "chunks": spacy_results["chunks"],
                }
                collection.insert_one(doc)
                
                metadata = {
                    "title": paper.get("title", ""),
                    "authors": paper.get("authors", []),
                    "year": int(year),
                    "source": "PubMed",
                    "abstract": abstract,
                    "keywords": paper.get("keywords", []),
                    "doi": paper.get("doi", ""),
                    "journal": paper.get("journal", ""),
                    "last_updated": paper.get("last_updated", ""),
                }
                
                # Save embeddings and metadata to Pinecone
                index.upsert(vectors=[(paper_id, spacy_results["embeddings"], metadata)])
        
        elif source == "Europe PMC":
            for paper in tqdm(papers, desc="Saving EuropePMC articles to MongoDB"):
                # Processar o abstract, se disponível
                abstract = paper.get("abstractText", "")
                spacy_results = process_text(abstract) if abstract else {
                    "entities": [], "matched_terms": {}, "chunks": [], "embeddings": []
                }

                # Transformar autores em lista de strings
                authors = paper.get("authorList", {}).get("author", [])
                authors = [f"{author.get('firstName', '')} {author.get('lastName', '')}" for author in authors]
                
                paper_id = generate_unique_id()  # Gerar um ID único para o artigo  
                   
                doc = {
                    "paper_id": paper_id,
                    "title": paper.get("title", ""),
                    "authors": authors,
                    "year": int(paper.get("pubYear", 0) or 0),
                    "source": "Europe PMC",
                    "abstract": abstract,
                    "keywords": paper.get("keywordList", {}).get("keyword", []),
                    "doi": paper.get("doi", ""),
                    "last_updated": paper.get("firstPublicationDate", ""),
                    "spacy_entities": spacy_results["entities"],
                    "spacy_matched_terms": spacy_results["matched_terms"],
                    "chunks": spacy_results["chunks"],
                }
                collection.insert_one(doc)
                
                # Save embeddings to Pinecone
                index.upsert(vectors=[(paper_id, spacy_results["embeddings"])])

        elif source == "Semantic Scholar":
            for paper in tqdm(papers, desc="Saving Semantic Scholar articles to MongoDB"):
                # Processar o abstract, se disponível
                abstract = paper.get("abstract", "")
                spacy_results = process_text(abstract) if abstract else {
                    "entities": [], "matched_terms": {}, "chunks": [], "embeddings": []
                }

                paper_id = generate_unique_id()  # Gerar um ID único para o artigo
                
                doc = {
                    "paper_id": paper_id,
                    "title": paper.get("title", ""),
                    "authors": [author.get("name", "") for author in paper.get("authors", [])],
                    "year": int(paper.get("year", 0) or 0),
                    "source": "Semantic Scholar",
                    "abstract": abstract,
                    "keywords": [],
                    "doi": paper.get("externalIds", {}).get("DOI", ""),
                    "journal": paper.get("journal", {}).get("name", "") if paper.get("journal") else "",
                    "last_updated": "",
                    "spacy_entities": spacy_results["entities"],
                    "spacy_matched_terms": spacy_results["matched_terms"],
                    "chunks": spacy_results["chunks"],
                }
                collection.insert_one(doc)
                
                # Save embeddings to Pinecone
                index.upsert(vectors=[(paper_id, spacy_results["embeddings"])])

        
        elif source == "GoogleScholar":
            for paper in tqdm(papers, desc="Saving Google Scholar articles to MongoDB"):
                # Processar o abstract, se disponível
                abstract = paper.get("abstract", "")
                spacy_results = process_text(abstract) if abstract else {
                    "entities": [], "matched_terms": {}, "chunks": [], "embeddings": []
                }

                authors = paper.get("authors", "No Authors")
                if isinstance(authors, list):
                    authors = ", ".join(authors)  # Combinar autores em uma string

                paper_id = generate_unique_id()  # Gerar um ID único para o artigo
                
                doc = {
                    "paper_id": paper_id,
                    "title": paper.get("title", ""),
                    "authors": authors,
                    "year": int(paper.get("year", 0) or 0),
                    "source": "Google Scholar",
                    "abstract": abstract,
                    "keywords": paper.get("keywords", []),
                    "doi": paper.get("doi", ""),
                    "journal": paper.get("journal", ""),
                    "last_updated": "",
                    "spacy_entities": spacy_results["entities"],
                    "spacy_matched_terms": spacy_results["matched_terms"],
                    "chunks": spacy_results["chunks"], 
                }
                collection.insert_one(doc)
                
                # Save embeddings to Pinecone
                index.upsert(vectors=[(paper_id, spacy_results["embeddings"])])

    print("All articles have been successfully saved to MongoDB!")