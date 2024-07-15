import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

# Update version from VERSION file into module
with open('VERSION', 'r') as fversion:
    version = fversion.readline().rstrip()
with open('_version.py', 'wt') as fversion:
    fversion.write('__version__ = "'+version+'"')

name = 'dypylib'
setuptools.setup(
    name=name, # Replace with your own username
    version=version,
    author='You Duan',
    author_email='duanyou@outlook.com',
    description='genomics research',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/yduanBioinfo/dypylib',
    install_requires = ['case-insensitive-dictionary'],
    #packages=setuptools.find_packages(),
    package_dir={name: '.'},
    packages=[name] + ['.'.join((name, x)) for x in setuptools.find_packages()],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
    #data_files = [
    #    ('',['VERSION'])
    #],
    scripts=[
        'scripts/subset_fasta.py',
        'scripts/subset_fastq.py',
        'scripts/annot_from_file.py',
        'scripts/filter.py',
        'bio/seq/plotsts.py',
        'bio/seq/ncbi_chr.py',
        'bio/seq/gtf2table.py',
        'bio/seq/countgtf.py',
        'bio/seq/gtf2stat.py',
        'bio/seq/format/merge_overlap_gene.py',
    ],
    python_requires='>=3.6',
)
