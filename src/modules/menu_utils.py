from rich.console import Console
from rich.prompt import Prompt
from rich.table import Table

# Lista de queries relevantes sobre prevenção farmacológica
VITAMIN_QUERIES = [
    '(TITLE:preventive OR ABSTRACT:preventive) AND pharmacology AND ("chronic disease" OR "chronic diseases")',
    '"primary prevention" AND (pharmaceuticals OR drugs OR medications)',
    'statins AND ("disease prevention" OR "cardiovascular prevention")',
    'metformin AND ("cancer prevention" OR "oncoprevention")',
    'aspirin AND ("primary prevention" OR "prophylactic use")',
    '("immunomodulatory drugs" OR "immunosuppressants") AND ("infection prevention" OR "prophylaxis")',
    'nutraceuticals AND pharmaceuticals AND ("disease prevention" OR "preventive health")',
    '"drug repurposing" AND ("disease prevention" OR prophylaxis)',
    'antihypertensives AND ("preventive use" OR "primary prevention")',
    'chemoprevention AND ("randomized controlled trial" OR RCT)'
]

def display_vitamin_queries():
    """Displays the list of vitamin queries in a formatted table."""
    console = Console()
    table = Table(title="🧪 Preventive Pharmacology Research Queries", show_header=True, header_style="bold magenta")
    table.add_column("ID", style="dim", width=3)
    table.add_column("Query", style="cyan")
    
    for i, query in enumerate(VITAMIN_QUERIES, 1):
        table.add_row(str(i), query)
    
    console.print(table)
    return len(VITAMIN_QUERIES)

def select_query():
    """Allows user to select a query from the list."""
    console = Console()
    total_queries = display_vitamin_queries()
    
    while True:
        try:
            choice = int(Prompt.ask(f"[bold white]Select a query (1-{total_queries}) or 0 for custom query[/bold white]"))
            if choice == 0:
                return Prompt.ask("[bold white]Enter your custom query:[/bold white]")
            elif 1 <= choice <= total_queries:
                return VITAMIN_QUERIES[choice - 1]
            else:
                console.print(f"[bold red]❌ Please enter a number between 0 and {total_queries}[/bold red]")
        except ValueError:
            console.print("[bold red]❌ Please enter a valid number[/bold red]")

def display_enhanced_menu():
    """Displays the enhanced menu with vitamin queries option."""
    console = Console()
    
    menu_table = Table(title="🔬 Enhanced Research Tool", show_header=True, header_style="bold green")
    menu_table.add_column("Option", style="bold blue", width=8)
    menu_table.add_column("Description", style="white")
    
    menu_options = [
        ("1", "PubMed Search"),
        ("2", "Europe PMC Search"),
        ("3", "Semantic Scholar Search"),
        ("4", "Wikipedia Search"),
        ("5", "Google Scholar Search"),
        ("6", "Search All Sources"),
        ("7", "🧪 Google Scholar Batch (Single Query)"),
        ("8", "📋 View Research Queries List"),
        ("9", "🚀 Execute ALL Research Queries (Batch Mode)"),
        ("T", "🔧 Test Google Scholar Connection"),
        ("Q", "Quit")
    ]
    
    for option, description in menu_options:
        menu_table.add_row(option, description)
    
    console.print(menu_table)
    return Prompt.ask("\n[bold white]Choose an option[/bold white]", default="1")

def display_source_selection_menu():
    """Displays the source selection menu for batch searches."""
    console = Console()
    
    source_table = Table(title="Choose Source for Batch Search", show_header=True, header_style="bold cyan")
    source_table.add_column("Option", style="bold blue", width=8)
    source_table.add_column("Source", style="white")
    
    source_options = [
        ("1", "Google Scholar"),
        ("2", "PubMed"),
        ("3", "Europe PMC"),
        ("4", "Semantic Scholar")
    ]
    
    for opt, src in source_options:
        source_table.add_row(opt, src)
    
    console.print(source_table)
    return source_options

def get_source_map():
    """Returns the mapping of menu choices to source names."""
    return {
        "1": "google_scholar",
        "2": "pubmed", 
        "3": "europe_pmc",
        "4": "semantic_scholar"
    }

def get_sources_dict():
    """Returns the sources dictionary for the main menu."""
    from modules.pubmed_utils import search_pubmed
    from modules.europePMC_utils import search_europe_pmc
    from modules.semanticscholar_utils import search_semanticscholar
    from modules.wikipedia_utils import search_wikipedia
    from modules.googleScholar_utils import search_google_scholar
    
    return {
        "1": ("PubMed", search_pubmed),
        "2": ("Europe PMC", search_europe_pmc),
        "3": ("Semantic Scholar", search_semanticscholar),
        "4": ("Wikipedia", search_wikipedia),
        "5": ("Google Scholar", search_google_scholar),
        "6": ("All Sources", None)
    }

def display_batch_confirmation(num_queries, selected_source, num_docs_per_query):
    """Displays batch operation confirmation details."""
    console = Console()
    
    console.print(f"\n[bold yellow]⚠️  About to execute {num_queries} queries on {selected_source.replace('_', ' ').title()}[/bold yellow]")
    console.print(f"[bold yellow]⚠️  This will retrieve up to {num_queries * num_docs_per_query} documents total[/bold yellow]")
    console.print(f"[bold yellow]⚠️  Estimated time: {num_queries * 2} seconds (with API delays)[/bold yellow]")
    
    return Prompt.ask("[bold white]Continue? (y/n)[/bold white]", default="n")

def display_menu():
    """Legacy function name for backward compatibility."""
    return display_enhanced_menu()