# submit data to topmed:

import requests
import json
import sys

# set study:
study = sys.argv[1]
print(study)

# imputation server url
url = 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2'


# add token to header (see Authentication)
token = 'YOUR_API_TOKEN'
headers = {'X-Auth-Token' : token }

# set parameters

# Parameter	Values	Default Value
# input-files	/path/to/file	
# input-mode	qconly, imputation	imputation
# input-password	user-defined password	auto
# input-files-source	file-upload, sftp, http	default: file-upload
# input-refpanel	apps@topmed-r2@1.0.0	-
# input-phasing	eagle	eagle
# input-population	all, mixed	all

# QC
data = {
  'input-mode': 'qconly',
  'input-refpanel': 'apps@topmed-r2@1.0.0',
  'input-password': 'YOUR_PASSWORD',
  'input-population': 'mixed',
  'input-phasing': 'eagle',
  'input-build': 'hg38'
}


##################
# submit new job

path = '/lustre/scratch123/hgi/projects/ibdgwas/IIBDGC/imputation_ready/'
ancestry = ['eur','noneur']

# EUROPEAN SAMPLES

vcfs = []
for chr in range(1,24):
  vcfs.append(path + study + '/' + study + '_' + ancestry[0] + '_chr' + str(chr) + '_2022.vcf.gz')


files = {
  ('input-files', open(vcfs[0], 'rb')),
  ('input-files', open(vcfs[1], 'rb')),
  ('input-files', open(vcfs[2], 'rb')),
  ('input-files', open(vcfs[3], 'rb')),
  ('input-files', open(vcfs[4], 'rb')),
  ('input-files', open(vcfs[5], 'rb')),
  ('input-files', open(vcfs[6], 'rb')),
  ('input-files', open(vcfs[7], 'rb')),
  ('input-files', open(vcfs[8], 'rb')),
  ('input-files', open(vcfs[9], 'rb')),
  ('input-files', open(vcfs[10], 'rb')),
  ('input-files', open(vcfs[11], 'rb')),
  ('input-files', open(vcfs[12], 'rb')),
  ('input-files', open(vcfs[13], 'rb')),
  ('input-files', open(vcfs[14], 'rb')),
  ('input-files', open(vcfs[15], 'rb')),
  ('input-files', open(vcfs[16], 'rb')),
  ('input-files', open(vcfs[17], 'rb')),
  ('input-files', open(vcfs[18], 'rb')),
  ('input-files', open(vcfs[19], 'rb')),
  ('input-files', open(vcfs[20], 'rb')),
  ('input-files', open(vcfs[21], 'rb')),
  ('input-files', open(vcfs[22], 'rb'))
}


r = requests.post(url + "/jobs/submit/imputationserver@1.6.6", files=files, headers=headers, data=data)
if r.status_code != 200:
  raise Exception('POST /jobs/submit/imputationserver@1.6.6 {}'.format(r.status_code))

# print message
print(r.json()['message'])
print("QC " + ' ' + study + ' job_' + ancestry[0] + '=' + r.json()['id'])

exit()

  

