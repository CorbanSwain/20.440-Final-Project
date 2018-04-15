import requests

if __name__ == "__main__":

    resp = requests.get('https://www.ncbi.nlm.nih.gov/pubmed/',
                        params={'term': 'PAXIP1-AS1'})

    file = 'text.xml'

    with open(file, 'w') as f:
        f.write(resp.text)


