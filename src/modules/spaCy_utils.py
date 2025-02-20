import spacy
from spacy.matcher import PhraseMatcher

# Load NLP spaCy model
nlp = spacy.load("en_core_web_trf")

# List of terms (diseases, supplements, pharmaceuticals) 
disease_terms = [
    "Cardiovascular diseases", "Cancer", "Diabetes", "Hypertension", "Obesity", 
    "Asthma", "Alzheimer's disease", "Parkinson's disease", "Rheumatoid arthritis", 
    "Osteoporosis", "Stroke", "COPD", "Depression", "Anxiety", "Insomnia", "Migraine", 
    "Epilepsy", "IBD", "Hepatitis", "Tuberculosis", "Malaria", "HIV/AIDS", "Mental health disorders", "disease prevention"
]

supplement_terms = [
    "Vitamin C", "Vitamin D", "Vitamin E", "Vitamin B12", "Calcium", "Magnesium", 
    "Omega-3 fatty acids", "Probiotics", "Fiber", "Zinc", "Iron", "Turmeric", 
    "Coenzyme Q10", "Glucosamine", "Chondroitin", "Ginseng", "Echinacea", "Green tea extract", 
    "Garcinia Cambogia", "Spirulina", "Aloe vera", "Ashwagandha", "L-carnitine", "Resveratrol", 
    "Melatonin", "Elderberry", "Fish oil", "Flaxseed oil", "Beta-glucan", "Biotin", "Collagen", "nutritional supplements"
]

pharmaceutical_terms = [
    "Acetaminophen", "Ibuprofen", "Aspirin", "Penicillin", "Amoxicillin", "Caffeine", 
    "Morphine", "Oxycodone", "Diazepam", "Lorazepam", "Atorvastatin", "Simvastatin", "Metformin",
    "Insulin", "Levothyroxine", "Prednisone", "Azithromycin", "Fluoxetine", "Sertraline", "pharmaceutical drugs"
]

# Combine all terms into one list
all_terms = disease_terms + supplement_terms + pharmaceutical_terms


def process_text(text):
    """Process text and extract entities and matched terms."""
    doc = nlp(text)
    
    # Extract entities
    entities = [(ent.text, ent.label_) for ent in doc.ents]
    
    # Match specific terms
    matcher = PhraseMatcher(nlp.vocab)
    patterns = [nlp.make_doc(term) for term in all_terms]
    matcher.add("SUPPLEMENT_AND_DISEASE", patterns)
    
    matches = matcher(doc)
    match_terms = [doc[start:end].text for _, start, end in matches]
    
    return {"entities": entities, "matched_terms": match_terms}
