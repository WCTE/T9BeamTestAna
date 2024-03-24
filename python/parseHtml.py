#!/usr/bin/python

#https://docs.python.org/3/library/html.parser.html
from html.parser import HTMLParser

class MyHTMLParser(HTMLParser):
    def handle_starttag(self, tag, attrs):
        print("Encountered a start tag:", tag)

    def handle_endtag(self, tag):
        print("Encountered an end tag :", tag)

    def handle_data(self, data):
        print("Encountered some data  :", data)

parser = MyHTMLParser()
infile = open('share/wcte-daq.html')
for xline in infile.readlines():
    line = xline[:-1]
    parser.feed(line)


    
    
