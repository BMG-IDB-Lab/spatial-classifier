from downloader import download_zhuang_abca
from reference_subsetting import create_meta_and_adata

files = download_zhuang_abca(output_dir="data")

adata_ref = create_meta_and_adata(RoI, IF_WRITE)
