import csv
import os

# convert CSV to TXT
def convert_csv_to_txt(csv_file, txt_file):
    with open(csv_file, 'r') as csv_in, open(txt_file, 'w') as txt_out:
        reader = csv.reader(csv_in)
        for row in reader:
            txt_out.write('\t'.join(row) + '\n')

csv_folder = 'waveform_data_csv'
txt_folder = 'waveform_data_txt'

if not os.path.exists(txt_folder):
    os.makedirs(txt_folder)

# CSV file to TXT
for filename in os.listdir(csv_folder):
    if filename.endswith('.CSV'):
        csv_file = os.path.join(csv_folder, filename)
        txt_file = os.path.join(txt_folder, os.path.splitext(filename)[0] + '.txt')
        convert_csv_to_txt(csv_file, txt_file)
