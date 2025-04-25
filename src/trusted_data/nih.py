import time
from datetime import datetime
import json
import re
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from bs4 import BeautifulSoup
import traceback

def format_url_name(name):
    formatted = re.sub(r'[^a-zA-Z0-9\s]', '', name.lower())
    formatted = formatted.replace(' ', '')
    return formatted

def dismiss_modal(driver):
    try:
        no_thanks_button = WebDriverWait(driver, 5).until(
            EC.element_to_be_clickable((By.ID, "prefix-dismissButton"))
        )
        no_thanks_button.click()
    except:
        pass


def check_page_exists(driver, timeout=5):
    try:
        WebDriverWait(driver, timeout).until(
            EC.presence_of_element_located((By.TAG_NAME, "body"))
        )
        
        # Verificar expressões específicas que indicam que a página não existe
        error_indicators = [
            "Page not found", 
            "404", 
            "The requested page could not be found",
            "We can't find the page"
        ]
        
        for indicator in error_indicators:
            if indicator in driver.page_source:
                return False
        
        # Verificar a presença de elementos que devem existir em páginas válidas
        try:
            article = WebDriverWait(driver, 3).until(
                EC.presence_of_element_located((By.TAG_NAME, "article"))
            )
            return True
        except:
            try:
                fact_sheet = WebDriverWait(driver, 3).until(
                    EC.presence_of_element_located((By.ID, "fact-sheet"))
                )
                return True
            except:
                h1_elements = driver.find_elements(By.TAG_NAME, "h1")
                if h1_elements:
                    return True
                return False
    except:
        return False
    
def clean_text(text):
    """Remove citações e \n """
    cleaned_text = re.sub(r'\[\s*\d+(?:\s*,\s*\d+)*\s*\]', '', text)
    cleaned_text = re.sub(r'\s+', ' ', cleaned_text).strip()
    cleaned_text = cleaned_text.replace('\n', ' ')
    return cleaned_text

def process_table(table_element):
    caption = table_element.find('caption')
    table_caption = caption.get_text().strip() if caption else "Tabela"
    
    # Obter linhas do cabeçalho
    headers = []
    thead = table_element.find('thead')
    if thead:
        for th in thead.find_all('th'):
            headers.append(clean_text(th.get_text().strip()))
    
    # Obter linhas de dados
    rows = []
    tbody = table_element.find('tbody')
    if tbody:
        tr_elements = tbody.find_all('tr')
    else:
        tr_elements = table_element.find_all('tr')
    
    for tr in tr_elements:
        row_data = [clean_text(td.get_text().strip()) for td in tr.find_all('td')]
        if not row_data:
            continue
        rows.append(row_data)
    
    text_output = [clean_text(table_caption)]
        
    # Converter cada linha em texto natural
    for row in rows:
        if len(row) >= 2:
            item_name = row[0]
            value_parts = []
            for i, cell in enumerate(row[1:], start=1):
                if i < len(headers):
                    value_parts.append(f"{headers[i]}: {cell}")
                else:
                    value_parts.append(cell)
            sentence = f"{item_name}: {', '.join(value_parts)}"
            text_output.append(sentence)
        else:
            text_output.append(f"• {row[0]}")
    
    return " ".join(text_output)

def extract_section_content(section_header):
    """Extrai o conteúdo de uma seção específica e retorna como texto formatado"""
    formatted_content = []
    current = section_header.find_next_sibling()
    
    if current and current.name == 'div':
        section_container = current
        elements = section_container.find_all(['p', 'table', 'ul', 'ol', 'h3'])
    else:
        elements = []
        while current and current.name != 'h2':
            if current.name in ['p', 'table', 'ul', 'ol', 'h3']:
                elements.append(current)
            current = current.find_next_sibling()
    
    for elem in elements:
        if elem.name == 'h3':
            formatted_content.append(f"{clean_text(elem.get_text().strip())}")
        elif elem.name == 'p':
            text = clean_text(elem.get_text())
            if text:
                formatted_content.append(text)
        elif elem.name == 'table':
            formatted_content.append(process_table(elem))
        elif elem.name in ['ul', 'ol']:
            list_items = []
            for li in elem.find_all('li'):
                marker = "• " if elem.name == 'ul' else "1. "
                list_items.append(f"{marker}{clean_text(li.get_text())}")
            formatted_content.extend(list_items)
    
    return " ".join(formatted_content)

def scrape_article_content(driver):
    """Extrai o conteúdo do artigo excluindo referências e disclaimers"""
    soup = BeautifulSoup(driver.page_source, 'html.parser')
    
    title_element = soup.find('h1', class_='fsTitle')
    title = title_element.get_text().strip() if title_element else "Sem título"
    
    article = soup.find('article') or soup.find('div', id='fact-sheet')
    
    if not article:
        print("Não foi possível encontrar o conteúdo do artigo.")
        return None
    
    for section_id in ['ref', 'disc', 'divCitations', 'divDisclaimer']:
        reference_sections = article.find_all(['section', 'div'], id=lambda x: x and section_id in x)
        for section in reference_sections:
            section.extract()
    
    content_data = {
        "title": title,
        "content": []  
    }
    
    section_headers = article.find_all('h2')
    
    for header in section_headers:
        section_id = header.get('id')
        section_title = clean_text(header.get_text().strip())
        
        if (section_id and (section_id in ['ref', 'disc'])) or \
           ('reference' in section_title.lower()) or \
           ('disclaimer' in section_title.lower()) or \
           ('table of contents' in section_title.lower()):
            continue
        
        section_content = f"{section_title} "
        
        section_content += extract_section_content(header)
        
        content_data["content"].append(section_content)
    
    content_data["full_content"] = " ".join(content_data["content"])
    
    return content_data

def scrape_supplement_page(url, supplement_name, driver):
    """Faz o scraping de uma página de suplemento"""
    print(f"\nNavegando para a URL: {url}")
    
    try:
        driver.get(url)
        final_url = driver.current_url

        if not check_page_exists(driver):
            print(f"A página para {supplement_name} não existe.")
            return None
        
        dismiss_modal(driver)
        
        try:
            WebDriverWait(driver, 10).until(
                EC.presence_of_element_located((By.TAG_NAME, "article"))
            )
        except:
            WebDriverWait(driver, 10).until(
                EC.presence_of_element_located((By.ID, "fact-sheet"))
            )
        
        content = scrape_article_content(driver)
        
        if content:
            simplified_content = {
                "title": supplement_name,
                "link": final_url,
                "source" : "NIH",
                "scraped_at": datetime.now().isoformat(),  
                "content": content["full_content"]
            }
            return simplified_content
        else:
            print(f"Não foi possível extrair o conteúdo para {supplement_name}")
            return None
            
    except Exception as e:
        print(f"Erro durante o scraping: {str(e)}")
        traceback.print_exc()
        return None

import os

def process_supplements(supplement_file, output_file="supplements_data.json"):
    with open(supplement_file, 'r', encoding='utf-8') as f:
        data = json.load(f)

    if os.path.exists(output_file):
        with open(output_file, 'r', encoding='utf-8') as f:
            saved_data = json.load(f)
    else:
        saved_data = {"supplements": []}

    processed_names = {s['title'] for s in saved_data.get("supplements", [])}

    chrome_options = Options()
    chrome_options.add_argument("--headless")
    chrome_options.add_argument("--disable-gpu")
    chrome_options.add_argument("--window-size=1920,1080")

    driver = webdriver.Chrome(options=chrome_options)

    supplements = data.get("supplement_terms", [])
    success_count = 0
    fail_count = 0

    try:
        for i, supplement in enumerate(supplements):
            if supplement in processed_names:
                print(f"\n[{i+1}/{len(supplements)}] Já processado: {supplement}.")
                continue

            print(f"\n[{i+1}/{len(supplements)}] A processar: {supplement}")
            url_name = format_url_name(supplement)
            url = f"https://ods.od.nih.gov/factsheets/{url_name}-HealthProfessional"

            content = scrape_supplement_page(url, supplement, driver)

            if content:
                saved_data["supplements"].append(content)
                with open(output_file, 'w', encoding='utf-8') as f:
                    json.dump(saved_data, f, ensure_ascii=False, indent=2)

                success_count += 1
                print(f"Scraping feito para {supplement}.")
            else:
                fail_count += 1
                print(f"Falha ao fazer o scraping de {supplement}.")

            time.sleep(2)

    except Exception as e:
        print(f"Erro ao processar suplementos: {str(e)}")
        traceback.print_exc()
    finally:
        driver.quit()

if __name__ == "__main__":
    supplement_file = "../terms/supplement.json"
    output_file = "supplements_data.json"
    
    process_supplements(supplement_file, output_file)