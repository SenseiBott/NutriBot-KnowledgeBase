import json
import os
from dotenv import load_dotenv
from tqdm import tqdm
from modules.pinecone_utils import save_to_mongo_and_pinecone
from pathlib import Path 

def load_trusted_data(json_file_path):
    """Load trusted data from JSON file."""
    try:
        with open(json_file_path, 'r', encoding='utf-8') as file:
            data = json.load(file)
            return data.get("supplements", [])
    except FileNotFoundError:
        print(f"Error: File {json_file_path} not found.")
        return []
    except json.JSONDecodeError:
        print(f"Error: File {json_file_path} contains invalid JSON.")
        return []
    except Exception as e:
        print(f"Error loading trusted data: {str(e)}")
        return []

def process_trusted_data(json_file_path):
    """Process trusted data from JSON file and save to Pinecone."""
    load_dotenv()  # Load environment variables
    
    print(f"Loading trusted data from {json_file_path}...")
    supplements = load_trusted_data(json_file_path)
    
    if not supplements:
        print("No trusted data found or error loading data.")
        return []
    
    print(f"Found {len(supplements)} supplements from trusted sources.")
    
    # Save data to Pinecone with level1 source identifier
    save_to_mongo_and_pinecone(supplements, "level1")
    
    return supplements

def save_results_to_json(data, filename="processed_trusted_data.json"):
    """Saves the processed trusted data to a JSON file."""
    with open(filename, "w", encoding="utf-8") as file:
        json.dump(data, file, ensure_ascii=False, indent=4)
    print(f"Results saved in {filename}")

def main():
    
    script_dir = Path(__file__).parent
    
    json_file_path = script_dir / "trusted_data" / "supplements_data.json"

    supplements = process_trusted_data(json_file_path)
    

if __name__ == "__main__":
    main()