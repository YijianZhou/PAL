import glob

def merge(fnames, out):
  for fname in fnames:
    f=open(fname); lines=f.readlines(); f.close()
    for line in lines:
      out.write(line)
  out.close()
  return

net='zsy'
ctlg_files = './output/%s/catalog_20*'%net
pha_files  = './output/%s/phase_20*'%net
ctlgs = sorted(glob.glob(ctlg_files))
phas  = sorted(glob.glob(pha_files))
ctlg_out = open('./output/catalog_xj_%s.dat'%net,'w')
pha_out  = open('./output/phase_xj_%s.dat'%net,'w')
merge(ctlgs, ctlg_out)
merge(phas, pha_out)

