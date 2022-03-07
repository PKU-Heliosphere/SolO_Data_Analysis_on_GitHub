from bs4 import BeautifulSoup
import requests
import os

# preparation
_dir = "https://spdf.gsfc.nasa.gov/pub/data/"
_sc = "psp"
_year = "2020"
_inst = "fields" # fields sweap
_payload = ""
_level = "l2"
_data = "mag_rtn"
_URL = _dir + _sc + "/" + _inst + "/" + _payload + "/" + _level + "/" + _data + "/" + _year + "/"
out_dir = "/Users/psr/Documents/cdfData/" + _sc + "/" + _inst + "/" + _payload + "/" + _level + "/" + _data + "/" + _year + "/"
os.makedirs(os.path.dirname(out_dir), exist_ok=True)

# functional
req = requests.get(_URL)
soup = BeautifulSoup(req.text, features="html.parser")
urls = []
names = []
links = []
for i, link in enumerate(soup.findAll('a')):
    _FULLURL = _URL + link.get('href')
    if _FULLURL.endswith('.cdf'):
        urls.append(_FULLURL)
        names.append(soup.select('a')[i].attrs['href'])
for i, (name, url) in enumerate(zip(names, urls)):
    if i >= 0:
        print(url)
        r = requests.get(url)
        with open(out_dir + name.split('/')[-1], "wb") as f:
            f.write(r.content)
