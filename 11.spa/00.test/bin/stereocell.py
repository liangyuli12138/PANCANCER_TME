import os
import datetime
import argparse
import glob
import subprocess


###### Version and Date
PROG_VERSION = 'StereoCell Python Wrapper: v1.1.0'
PROG_DATE = '2023-08-25'

###### Usage
usage="""StereoCell Wrapper v1.1.0"""


def parser_arguments():
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument("--version", action="version", version=PROG_VERSION)
    parser.add_argument("-i", "--stich_images",  dest="stich_images", type=str, required=True, help="Input stich_images file.")
    parser.add_argument("-o", "--output",  dest="output", type=str, required=True, help="Result output dir.")
    parser.add_argument("-c", "--chipNo", dest="chipNo", type=str, default='', help="chip number.")
    parser.add_argument("-g", "--gem", dest="gem", type=str, default='', help="gem file (*.gef, *.gem, *.gem.gz)")
    parser.add_argument("-p", "--stereocell_path", dest="stereocell_path", type=str, default='/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/11.spa/00.test/bin/', help="path to stereocell_1_1_0 directory")

    args = parser.parse_args()
    return (args)


def _install_whl(stereocell_path):
    # 安装stio whl
    print('install ....')
    commands = ['pip', 'install', f'{stereocell_path}/stereocell_1_1_0/whl/stio-0.1.0-py3-none-any.whl', '--no-dependencies']
    subprocess.run([str(i) for i in commands])
    # !pip install stereocell_1_1_0/whl/stio-0.1.0-py3-none-any.whl --no-dependencies
    # 安装cellbin whl
    commands = ['pip', 'install', f'{stereocell_path}/stereocell_1_1_0/whl/cell_bin-1.1.0-py3-none-any.whl', '--no-dependencies']
    subprocess.run([str(i) for i in commands])
    # !pip install stereocell_1_1_0/whl/cell_bin-1.1.0-py3-none-any.whl --no-dependencies


def _main_qc(output_path, stich_images_path, chipNo, stereocell_path, stainType='ssDNA'):

    # 用户需要检查
    qc_input_json = {
        'output': output_path,  # qc 输出位置 '/data/work/previous/01.DEAF/25.cell_bin/cellbin_stereoCell/v1.1.1/qc_output/FP200000571BR_C5_1'
        'image': stich_images_path,  # qc 输入图片路径 '/data/work/previous/01.DEAF/25.cell_bin/cell_bin_v3/FP200000571BR_C5_1/registration/7_result/FP200000571BR_C5_1_stitched.tif'
        'chipNo': chipNo,  # 芯片号 'FP200000571BR_C5_1'
        'stainType': stainType,  # 染色类型 HE
        'config': f'{stereocell_path}/stereocell_1_1_0/py/qc_config.json',  # qc配置文件
        'magnification': 10,  # 显微镜物镜倍率，目前只支持10倍镜
        'zoo_dir': f'{stereocell_path}/stereocell_1_1_0/weights',  # 权重路径
        'fov_size': (2000, 2000),  # 大图裁剪尺寸
        'debug_mode': False, 
    }
    print(qc_input_json)
    main_qc(
    stain_type=qc_input_json['stainType'], 
    src_data=qc_input_json['image'], 
    chip_name=qc_input_json['chipNo'], 
    save_dir=qc_input_json['output'],
    config_path=qc_input_json['config'],
    zoo_dir=qc_input_json['zoo_dir'],
    debug_mode=qc_input_json['debug_mode'], 
    fov_size=qc_input_json['fov_size'],
    magnification=qc_input_json['magnification']
    )
    
    return qc_input_json


def _run_pipeline(qc_input_json, output_path, gem_path, stereocell_path):
    import glob
    import os

    # 找到并使用最新的qc结果
    search_dir = qc_input_json['output']
    files = list(filter(os.path.isfile, glob.glob(search_dir + "/*.ipr")))
    files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
    the_most_recent_ipr = files[0]
    print(f"using {the_most_recent_ipr}")
    
    pipeline_input_json = {
    'output': output_path,  # qc 输出位置
    'config': f'{stereocell_path}/stereocell_1_1_0/py/pipeline_config.json',  # qc配置文件 '/data/work/previous/01.DEAF/25.cell_bin/cellbin_stereoCell/v1.1.1/pipeline_output/FP200000571BR_C5_1'
    'gene': gem_path,  # 矩阵文件 '/data/work/previous/01.DEAF/25.cell_bin/data/raw_data/raw_gem/FP200000571BR_C5/FP200000571BR_C5.gem.gz'
    'cell_type': 'CELL'  
    }
    print(pipeline_input_json)
    # pipeline
    from stereocell_1_1_0.py.pipeline import run_pipeline

    run_pipeline(
        input_=qc_input_json['image'], 
        output=pipeline_input_json['output'], 
        ipr_path=the_most_recent_ipr, 
        matrix=pipeline_input_json['gene'], 
        cell_type=pipeline_input_json['cell_type'],
        config=pipeline_input_json['config'], 
        zoo_dir=qc_input_json['zoo_dir'],
    )


def _log(info):
    current_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(' '.join([current_time, info]))    

    
if __name__ == '__main__':
    print(usage)
    args = parser_arguments()
    
    args.output = os.path.join(args.output, args.chipNo)
    if not os.path.exists(args.output):os.makedirs(args.output)

    _log('Install whl...')
    _install_whl(stereocell_path=args.stereocell_path)
 
    from stereocell_1_1_0.py.qc import main_qc

    _log('Main qc...')
    qc_input_json = _main_qc(output_path=args.output, stich_images_path=args.stich_images, chipNo=args.chipNo, stereocell_path=args.stereocell_path)

    _log('Run pipeline...')
    _run_pipeline(qc_input_json=qc_input_json, output_path=args.output, gem_path=args.gem, stereocell_path=args.stereocell_path)
