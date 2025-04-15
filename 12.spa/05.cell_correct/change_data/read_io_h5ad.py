import sys
sys.path.append('/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/liuhui/script')
input_file = sys.argv[1]
output_file = sys.argv[2]
print('input file:',input_file)
print('output file:',output_file)
import read_io as gem
adata = gem.read_gem(input_file,label_column="CellID", n_workers=10)
adata.write_h5ad(output_file)
