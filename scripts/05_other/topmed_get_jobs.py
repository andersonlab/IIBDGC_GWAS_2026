import requests
import json

# imputation server url
url = 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2'
token = 'YOUR_API_TOKEN';

# add token to header (see authentication)
headers = {'X-Auth-Token' : token }

# get all jobs
r = requests.get(url + "/jobs", headers=headers)
if r.status_code != 200:
    raise Exception('GET /jobs/ {}'.format(r.status_code))

# print all jobs
for job in r.json():
    print('{} [{}]'.format(job['id'], job['state']))
