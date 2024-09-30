import re
import time
import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.wait import WebDriverWait
import os
###################################################
# This file describes a simple pipeline for downloading THE FIRST result of a HADDOCK run; 
# It takes as entry a list of links from the HADDOCK_LINKS file; 
# Using Seleinum it scrapes the first file from the page; 
###################################################


STRUCTURES_PATH = '' #Define the path to store your structures;
HADDOCK_LINKS = '' #Define the path to you file containing your result links - provided via email or in the submission page after submission; 

WAIT_TIME=64800 #18 hours -> Wait time high to allow a docking run to finish, this time can be modified to your use case;
WAIT_TIME02=70 # This wait time relates to the PDB dropdown button - sometimes takes a while to load (or just breaks PAY ATTENTION HERE)
def retrieve_best(driver, url):
    driver.get(url)
    print("Waiting for page to load...")
    time.sleep(5)
    table = WebDriverWait(driver, WAIT_TIME).until(EC.presence_of_element_located((By.XPATH, '//*[@id="content"]/div[2]')))
    content = table.text.split("\n")[9].split(" ")[1]
    #print(f'ZDOCK: {content}')
    target = WebDriverWait(driver, WAIT_TIME).until(EC.presence_of_element_located((By.XPATH, '//*[@id="content"]/div[2]')))
    driver.execute_script("arguments[0].scrollIntoView();", target)
    structure = WebDriverWait(driver, WAIT_TIME).until(EC.presence_of_element_located((By.XPATH, '//*[@id="dropdownMenuLink"]')))
    driver.execute_script("arguments[0].click();", structure)
    print("Found PDB dropdown!")
    time.sleep(2)
    print("Now waiting for Options...")
    pdb = WebDriverWait(driver, WAIT_TIME02).until(EC.element_to_be_clickable((By.XPATH, '//*[@id="content"]/div[2]/table[2]/tbody/tr[1]/td[3]/div/div/li[1]/a')))
    driver.execute_script("arguments[0].click();", pdb)
    time.sleep(3)
    print("PDB downloaded successfully")
    if content == '': content = 999 #In case of failure to retrieve ZDOCK value, set it to an ABSURD value; 
    else: pass
    return content

docking_dict = {"ENTRY_ID":[], "ZDOCK_VALUE":[]}


with open(f"{HADDOCK_LINKS}", "r") as op:
    data = op.read().split("\n")[:-1]
data.reverse()

driver = webdriver.Chrome()
for link  in data:
    idx = link.split("-")[1]
    print(f"Currently attempting link: {link} || and PDB:{idx}")
    content = retrieve_best(driver, link)
    os.system(f'mv ~/Downloads/*.pdb ~/{STRUCTURES_PATH}/{idx}.pdb')
    docking_dict["ENTRY_ID"].append(idx)
    docking_dict["ZDOCK_VALUE"].append(content)
    print(f">{idx} {content}")
#Could use some modification to store idx and zdock values in a csv file
driver.close()
    