# Configure MongoDB connection
from pymongo import MongoClient

def configure_mongoDB_connection():
    """Configure MongoDB connection."""
    client = MongoClient("mongodb://localhost:27017/")  
    db = client["supplement_agent"]  
    collection = db["papers"] 
    return collection


