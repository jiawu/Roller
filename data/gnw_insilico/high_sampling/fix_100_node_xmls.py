from lxml import html, etree
import os
import sys

def read_broken_xml(filepath):
    """
    Read the file in, convert it to an element tree
    :param filepath:
    :return:
    """
    with open(filepath, 'r') as fn:
        xmldata = ''.join(fn.readlines())

    # Read as an html
    element_data = html.fromstring(xmldata)

    # Convert to string, then read from string. This prevents weird headers being added
    element = etree.fromstring(etree.tostring(element_data))

    tree = etree.ElementTree(element)
    return tree

def change_case(tree):
    """
    Change to lower camel case. Operates in place
    :param tree:
    :return:
    """
    for e in tree.iter():
        name = (e.tag[e.tag.index('}') + 1:])
        newname = ""
        if 'of' in name:
            splits = name.split('of')
            newname = 'listOf' + splits[-1].title()
        elif 'reference' in name:
            newname = name.replace('reference', 'Reference')
            if 'modifier' in name:
                newname = newname.replace('species', 'Species')
        elif 'law' in name:
            newname = name.replace('law', 'Law')

        if newname:
            e.tag = e.tag.replace(name, newname)


if __name__ == '__main__':
    directory = './Ecoli100/'
    files = os.listdir(directory)
    for file in files:
        if 'xml' in file:
            filepath = directory + file
            try:
                newtree = read_broken_xml(filepath)
            except:
                newtree = etree.parse(filepath)

            change_case(newtree)
            newtree.write(filepath, xml_declaration=True, encoding='UTF-8')

