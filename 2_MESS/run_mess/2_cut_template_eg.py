""" Cut template data
"""
import os, shutil

# i/o paths
mess_dir = '/home/zhouyj/software/2_MESS'
data_dir = '/data/Example_data'
out_root = 'output/Example_templates'
temp_pha = 'input/eg_pal.temp'
cut_method = 'intense' # 'intense' or 'long'

shutil.copyfile('config_eg.py', os.path.join(mess_dir, 'config.py'))
os.system("python {}/cut_template_{}.py \
    --data_dir={} --temp_pha={} --out_root={}"\
    .format(mess_dir, cut_method, data_dir, temp_pha, out_root))

