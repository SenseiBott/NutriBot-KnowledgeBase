import os
import uuid
import numpy as np
from tqdm import tqdm
from modules.spaCy_utils import process_text
from pinecone import Pinecone, ServerlessSpec

def configure_pinecone_connection():
    """Configure Pinecone connection."""
    pinecone_api_key = os.getenv("PINECONE_API_KEY")
    print(f"PINECONE_API_KEY: {pinecone_api_key}")  # Debugging line to check API key
    if not pinecone_api_key:
        raise ValueError("PINECONE_API_KEY não está definida.")
    pc = Pinecone(api_key=pinecone_api_key)
    index_name = "project" 
    if index_name not in pc.list_indexes().names():
        pc.create_index(
            name=index_name,
            dimension=1024,  
            metric='cosine',
            spec=ServerlessSpec(cloud='aws', region='us-east-1')
        )
    return pc.Index(index_name)

def generate_unique_id():
    """Generate a unique ID for each paper."""
    return str(uuid.uuid4())

def extract_paper_attributes(paper, source):
    """Extract paper attributes based on the source API."""
    if source == "PubMed":
        year = paper.get("year", 0)
        if year == "No Year Available":
            year = 0
        return {
            "title": paper.get("title", ""),
            "authors": paper.get("authors", []),
            "year": int(year),
            "source": "PubMed",
            "abstract": paper.get("abstract", ""),
            "keywords": paper.get("keywords", []),
            "doi": paper.get("doi", ""),
            "journal": paper.get("journal", ""),
            "last_updated": paper.get("last_updated", "")
        }
    elif source == "Europe PMC":
        authors = paper.get("authorList", {}).get("author", [])
        authors = [f"{author.get('firstName', '')} {author.get('lastName', '')}" for author in authors]
        return {
            "title": paper.get("title", ""),
            "authors": authors,
            "year": int(paper.get("pubYear", 0) or 0),
            "source": "Europe PMC",
            "abstract": paper.get("abstractText", ""),
            "keywords": paper.get("keywordList", {}).get("keyword", []),
            "doi": paper.get("doi", ""),
            "journal": "",
            "last_updated": paper.get("firstPublicationDate", "")
        }
    elif source == "Semantic Scholar":
        return {
            "title": paper.get("title", ""),
            "authors": [author.get("name", "") for author in paper.get("authors", [])],
            "year": int(paper.get("year", 0) or 0),
            "source": "Semantic Scholar",
            "abstract": paper.get("abstract", ""),
            "keywords": [],
            "doi": paper.get("externalIds", {}).get("DOI", ""),
            "journal": paper.get("journal", {}).get("name", "") if paper.get("journal") else "",
            "last_updated": ""
        }
    elif source == "Google Scholar":
        authors = paper.get("authors", "No Authors")
        if isinstance(authors, list):
            authors = ", ".join(authors)  # Combine authors into a string
        return {
            "title": paper.get("title", ""),
            "authors": authors,
            "year": int(paper.get("year", 0) or 0),
            "source": "Google Scholar",
            "abstract": paper.get("abstract", ""),
            "keywords": paper.get("keywords", []),
            "doi": paper.get("doi", ""),
            "journal": paper.get("journal", ""),
            "last_updated": ""
        }
    elif source == "level1":
        # Extract data for trusted sources (hierarchy level 1)
        return {
            "title": paper.get("title", ""),
            "source": paper.get("source", "Unknown Source"),
            "link": paper.get("link", ""),
            "content": paper.get("content", ""),
            "scraped_at": paper.get("scraped_at", "")
        }
    else:
        raise ValueError(f"Unsupported source: {source}")

def save_paper_to_mongo_and_pinecone(paper, source, index):
    """Save a single paper to Pinecone with the new vector structure."""
    # Extract paper attributes
    paper_data = extract_paper_attributes(paper, source)
    
    # Process the abstract with spaCy
    if source == "level1":
        content = paper_data["content"]
        # Process the content with spaCy
        spacy_results = process_text(content) if content else {
            "entities": [], "matched_terms": {}, "chunks": [], "embeddings": np.zeros((0, 1024))
        }
    else:
        abstract = paper_data["abstract"]
        # Process the abstract with spaCy
        spacy_results = process_text(abstract) if abstract else {
            "entities": [], "matched_terms": {}, "chunks": [], "embeddings": np.zeros((0, 1024))
        }

    # Generate a unique paper ID
    paper_id = generate_unique_id()
    
    # Save embeddings and chunk text to Pinecone with the new structure
    embeddings = spacy_results["embeddings"]  # Shape: [n_chunks, 384]
    chunks = spacy_results["chunks"]
    if len(embeddings) != len(chunks):
        print(f"Warning: Mismatch between embeddings ({len(embeddings)}) and chunks ({len(chunks)}) for paper {paper_id}")
        return

    vectors = []
    for i, (embedding, chunk) in enumerate(zip(embeddings, chunks)):
        if embedding.shape != (1024,):
            print(f"Invalid embedding shape for chunk {i} of paper {paper_id}: {embedding.shape}")
            continue
        chunk_id = f"{paper_id}_chunk_{i}"

        # Define hierarchy level based on source
        hierarchy_level = 1 if source == "level1" else 2

        # New metadata structure
        if source == "level1":
            metadata = {
                "text": chunk,
                "title": paper_data["title"],
                "link": paper_data["link"],
                "source": paper_data["source"],  # Added source information for level1 (NIH, FDA, etc.)
                "topic": paper_data["title"].split()[0] if paper_data["title"] else "",  # Use first word of title as topic
                "hierarchy": hierarchy_level
            }
        else:
            metadata = {
                "text": chunk,  # 'text' replaces 'chunk_text'
                "title": paper_data["title"],
                "link": paper_data.get("doi", ""),  # Use DOI as link, or empty string if not available
                "year": str(paper_data["year"]) if paper_data["year"] else "",  # Convert year to string, handle None/0
                "topic": paper_data.get("keywords", [])[0] if paper_data.get("keywords") else "",  # Use first keyword as topic, or empty string
                "hierarchy": 2  # Static value for now, as no hierarchical level is provided
            }
        vectors.append({
            "id": chunk_id,
            "values": embedding.tolist(),
            "metadata": metadata
        })

    # Upsert vectors to Pinecone in batch with namespace "ns1"
    if vectors:
        index.upsert(vectors=vectors, namespace="ns1")

def save_to_mongo_and_pinecone(papers, source):
    """Save articles to Pinecone."""
    if not papers:
        print("No articles to save.")
        return
    
    # Configure connections
    index = configure_pinecone_connection()
    
    # Process each paper
    for paper in tqdm(papers, desc=f"Saving {source} articles"):
        save_paper_to_mongo_and_pinecone(paper, source, index)
    
    print(f"All {source} articles have been successfully saved!")