import spacy
from spacy.matcher import PhraseMatcher
import re

# Load the SciBERT model from spaCy
nlp = spacy.load("en_core_sci_scibert")

# Function to normalize text
def normalize_text(text):
    """
    Normalizes text: removes punctuation, extra spaces, and converts to lowercase.
    """
    text = text.lower()  # Convert to lowercase
    text = re.sub(r'[^\w\s]', '', text)  # Remove punctuation
    text = re.sub(r'\s+', ' ', text)  # Remove extra spaces
    return text.strip()

# Function to load pharmaceutical terms from a .txt file
def load_pharmaceuticals_from_file(filename):
    """Loads pharmaceutical names from a text file and normalizes the terms."""
    with open(filename, 'r', encoding='utf-8') as file:
        pharmaceuticals = file.readlines()
    return [normalize_text(pharmaceutical) for pharmaceutical in pharmaceuticals]

# Load and normalize the list of pharmaceuticals
pharmaceutical_terms = load_pharmaceuticals_from_file("data/drugs.txt")

# Lists of other biomedical terms (diseases, supplements, and medical concepts)
disease_terms = [
    "cardiovascular diseases", "cancer", "diabetes", "hypertension", "obesity", 
    "asthma", "alzheimer's disease", "parkinson's disease", "rheumatoid arthritis", 
    "osteoporosis", "stroke", "copd", "depression", "anxiety", "insomnia", "migraine", 
    "epilepsy", "ibd", "hepatitis", "tuberculosis", "malaria", "hiv/aids", "mental health disorders", "disease prevention"
]

supplement_terms = [
    "vitamin c", "vitamin d", "vitamin e", "vitamin b12", "calcium", "magnesium", 
    "omega-3 fatty acids", "probiotics", "fiber", "zinc", "iron", "turmeric", 
    "coenzyme q10", "glucosamine", "chondroitin", "ginseng", "echinacea", "green tea extract", 
    "garcinia cambogia", "spirulina", "aloe vera", "ashwagandha", "l-carnitine", "resveratrol", 
    "melatonin", "elderberry", "fish oil", "flaxseed oil", "beta-glucan", "biotin", "collagen", "nutritional supplements"
]

medical_concept_terms = [
    "immune function", "oxidative stress", "antioxidant properties", "cellular health", "immune response"
]

# Normalize all terms
disease_terms = [normalize_text(term) for term in disease_terms]
supplement_terms = [normalize_text(term) for term in supplement_terms]
medical_concept_terms = [normalize_text(term) for term in medical_concept_terms]

# Function to create a PhraseMatcher
def create_matcher(nlp, terms):
    """Creates a PhraseMatcher with normalized terms."""
    matcher = PhraseMatcher(nlp.vocab, attr="LOWER")  # Case-insensitive matching
    patterns = [nlp.make_doc(term) for term in terms]
    matcher.add("TERM_MATCHER", patterns)
    return matcher

# Create matchers for diseases, supplements, pharmaceuticals, and medical concepts
disease_matcher = create_matcher(nlp, disease_terms)
supplement_matcher = create_matcher(nlp, supplement_terms)
pharmaceutical_matcher = create_matcher(nlp, pharmaceutical_terms)
medical_concept_matcher = create_matcher(nlp, medical_concept_terms)

# Function to process the text and extract entities and corresponding terms
def process_text(text):
    """Processes the text and extracts entities and corresponding terms."""
    # Normalize the input text
    normalized_text = normalize_text(text)
    doc = nlp(normalized_text)

    # Extract named entities (NER)
    entities = [(ent.text, ent.label_) for ent in doc.ents]

    # Extract matched terms
    disease_matches = [doc[start:end].text for _, start, end in disease_matcher(doc)]
    supplement_matches = [doc[start:end].text for _, start, end in supplement_matcher(doc)]
    pharmaceutical_matches = [doc[start:end].text for _, start, end in pharmaceutical_matcher(doc)]
    medical_concept_matches = [doc[start:end].text for _, start, end in medical_concept_matcher(doc)]

    # Categorize entities with priority
    categorized_entities = []
    for ent_text, ent_label in entities:
        if ent_text in pharmaceutical_matches:
            categorized_entities.append((ent_text, "PHARMACEUTICAL"))
        elif ent_text in disease_matches:
            categorized_entities.append((ent_text, "DISEASE"))
        elif ent_text in supplement_matches:
            categorized_entities.append((ent_text, "SUPPLEMENT"))
        elif ent_text in medical_concept_matches:
            categorized_entities.append((ent_text, "MEDICAL_CONCEPT"))
        else:
            # Ignore generic terms like "essential", "patients", etc.
            continue

    return {
        "entities": categorized_entities,
        "matched_terms": {
            "diseases": disease_matches,
            "supplements": supplement_matches,
            "pharmaceuticals": pharmaceutical_matches,
            "medical_concepts": medical_concept_matches
        }
    }
