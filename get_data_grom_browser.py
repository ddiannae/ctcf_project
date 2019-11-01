"""Functions to web-scrap TCGA top mutations across age gaps"""

import time
import pandas as pd
import numpy as np
import os
from splinter import Browser


def gdc_explore(browser, cancer, age_from, age_to):
    """Go to GDC explore section"""
    # Revisit each time

    browser.visit(
        'https://portal.gdc.cancer.gov/exploration?searchTableTab=genes')
    # Open Cancer Type Menu (and leave it open)
    agg = browser.find_by_css('div[class="test-term-aggregation css-kfzjep"]')
    button = agg.find_by_xpath('./div[1]/div[6]/div')
    button.click()
    time.sleep(1)

    # Set the cancer primary site
    tag_id = browser.find_by_css('input[id="input-Primary Site-{}"]'.format(
        cancer))
    button = tag_id.find_by_xpath('..')
    button.click()
    time.sleep(1)
    # Set age gap in the form [x, y]
    # Age from
    x = browser.find_by_css(
        'input[id="from-cases.diagnoses.age_at_diagnosis"]').first
    x.fill(age_from)
    time.sleep(1)
    # Age to
    y = browser.find_by_css(
        'input[id="to-cases.diagnoses.age_at_diagnosis"]').first
    y.fill(age_to)
    parent = y.find_by_xpath('..')
    go = parent.find_by_css('a').first
    go.click()
    time.sleep(13)
    drop = browser.find_by_css(
        'span[class="test-page-size-selection dropdown"]')
    if not len(drop):
        return pd.DataFrame()
    button = drop.find_by_xpath('./span[1]/div[1]')
    button.click()
    time.sleep(1)
    drop = browser.find_by_css(
        'span[class="test-page-size-selection dropdown"]')
    # 100 entries
    mark = drop.find_by_xpath('./div[1]/div[6]').first
    mark.click()
    time.sleep(8)

    table = browser.find_by_id('genes-table').first.outer_html
    df = pd.read_html(table)[0]

    return df


if __name__ == '__main__':
    os.environ["PATH"] += os.pathsep + '/home/dianae/Workspace/regulacion_transcripcion'
    browser = Browser()
    # Visit main website
    browser.visit(
        'http://insulatordb.uthsc.edu/browse_new.php')

    # Select the specie
    sel_specie = browser.find_by_name('Species_type').first
    sel_specie.select("hg19")
    time.sleep(1)

    # Select cell chr
    sel_chr = browser.find_by_xpath('//select[@name="chr_symbol"]//option[@value="ALL"]').first
    sel_chr.click()
    time.sleep(1)

    # Select cell type
    sel_cell = browser.find_by_xpath('//select[@name="cell_symbol"]//option[@value="MCF-7"]').first
    sel_cell.click()
    time.sleep(1)

    #Set 10000 records per page
    inpp = browser.find_by_name("my_perpage").first
    inpp.fill("10000")
    time.sleep(10)

    # Cancers
    #cancers = np.loadtxt(
    #    '/Users/rdora/tcga/cancer_types.txt', dtype='U', delimiter='\n')
    #ages = [(a, a + 9) for a in range(10, 81, 10)]
    #for cancer in cancers[65:]:
    #    name = cancer.lower()
    #    name = name.replace(',', '')
    #    name = "-".join(name.split())
    #    cancer = '-'.join(cancer.split())
    #    print(cancer)
    #    for age in ages:
    #        age_name = str(age[0]) + '-' + str(age[1])
    #        print(age_name)
    #        table = gdc_explore(browser, cancer, *age)
    #        if table.shape[0]:
    #            table.to_csv('/Users/rdora/tcga/tables/{}_{}.csv'.format(
    #                name, age_name), index=False)
    #        else:
    #            print('There is no people aged {} to {} in {} Cancer'.format(
    #                age[0], age[1], cancer))
    browser.quit()
