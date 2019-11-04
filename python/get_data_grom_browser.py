"""Functions to web-scrap TCGA top mutations across age gaps"""

import time
import pandas as pd
import numpy as np
import os
from splinter import Browser

def getTable(browser):
    table = browser.find_by_xpath('//table[2]').first.outer_html
    df = pd.read_html(table)[0]
    return df

if __name__ == '__main__':
    os.environ["PATH"] += os.pathsep + '../'
    browser = Browser()
    # Visit main website
    browser.visit(
        'http://insulatordb.uthsc.edu/browse_new.php')

    # Select the specie
    sel_specie = browser.find_by_name('Species_type').first
    sel_specie.select("hg19")
    time.sleep(1)

    # Select cell type
    sel_cell = browser.find_by_xpath('//select[@name="cell_symbol"]//option[@value="MCF-7"]').first
    sel_cell.click()
    time.sleep(1)

    # Select cell chr
    sel_chr = browser.find_by_xpath('//select[@name="chr_symbol"]//option[@value="ALL"]').first
    sel_chr.click()
    time.sleep(3)

    dflist = []

    for i in range(221):
        print("Iteration: " + str(i) + "/n")
        df = getTable(browser)
        print(df.shape)

        while df.shape[0] != 1000:
            df = getTable(browser)
            print(df.shape)
            time.sleep(2)
        dflist.append(df)

        nlink = browser.find_by_xpath('//input[@type="image"][@src="right.gif"]').first
        nlink.click()
        time.sleep(5)

    df = getTable(browser)
    print(df.shape)
    dflist.append(df)

    dfresult = pd.concat(dflist)
    dfresult.to_csv('../data/ctcfbs.csv', index=False)

    browser.quit()
